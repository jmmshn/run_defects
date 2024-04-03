"""Builders related to MP specific tasks."""

from __future__ import annotations

from itertools import chain
from typing import TYPE_CHECKING

from maggma.builders import Builder
from maggma.stores import MemoryStore, Store
from monty.json import MontyDecoder
from pymatgen.analysis.structure_matcher import StructureMatcher

if TYPE_CHECKING:
    from collections.abc import Generator

mdecoder = MontyDecoder().process_decoded


class MPIDTagger(Builder):
    """Tag TASK documents with the MPID."""

    def __init__(
        self,
        job_store: Store,
        mp_materials: Store,
        query: dict | None = None,
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5,
        **kwargs,
    ) -> None:
        """Initialize the MPIDTagger.

        Args:
            job_store: Store of task documents.
            mp_materials: Store of materials documents.
            query: Query to limit the tasks to tag.
            ltol: Lattice parameter tolerance.
            stol: Site tolerance.
            angle_tol: Angle tolerance.
            kwargs: Other kwargs to pass to the parent class.
        """
        self.job_store = job_store
        self.mp_materials = mp_materials
        self.query = query or {}
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol
        super().__init__(
            sources=[job_store, mp_materials], targets=MemoryStore(), **kwargs
        )

    def get_items(self) -> Generator:
        """Get all tasks with a matching formula and symmetry keys."""
        self.logger.info("MPIDTagger: Getting items")
        # return groups of ([materials_doc,], [tasks_doc,])
        # group by formula_pretty
        self.job_store.distinct("output.formula_pretty", self.query)
        q = {
            "output.formula_pretty": {"$exists": True},
            "output.calcs_reversed.0.task_type": {"$exists": True},
            "metadata.material_id": {"$exists": False},
        }
        q.update(self.query)
        agg_pipe = self.job_store.docs_store._collection.aggregate(
            [
                {"$match": q},
                {
                    "$group": {
                        "_id": {
                            "formula_pretty": "$output.formula_pretty",
                            "symbol": "$output.symmetry.symbol",
                            "point_group": "$output.symmetry.point_group",
                            "number": "$output.symmetry.number",
                        },
                        "docs": {
                            "$push": {
                                "_id": "$_id",
                                "uuid": "$uuid",
                                "formula_pretty": "$output.formula_pretty",
                                "symbol": "$output.symmetry.symbol",
                                "point_group": "$output.symmetry.point_group",
                                "number": "$output.symmetry.number",
                                "structure": "$output.structure",
                            }
                        },
                    }
                },
            ]
        )
        for cur in agg_pipe:
            _id = cur["_id"]
            job_docs = cur["docs"]
            mp_q = {
                "formula_pretty": _id["formula_pretty"],
                "symmetry.symbol": _id["symbol"],
                "symmetry.number": _id["number"],
                "symmetry.point_group": _id["point_group"],
            }
            mp_docs = list(
                self.mp_materials.query(mp_q, properties=["material_id", "structure"])
            )
            yield mp_docs, job_docs

    def process_item(self, item: tuple) -> list[dict]:
        """Group the documents by structure similarity."""
        self.logger.info(
            f"Processing {len(item[1])} tasks for {item[1][0]['formula_pretty']}"
        )
        mp_docs, job_docs = item
        sm = StructureMatcher(ltol=self.ltol, stol=self.stol, angle_tol=self.angle_tol)
        for job_doc in job_docs:
            j_structure = mdecoder(job_doc["structure"])
            for mp_doc in mp_docs:
                mp_structure = mdecoder(mp_doc["structure"])
                if sm.fit(j_structure, mp_structure) or sm.fit(
                    mp_structure, j_structure
                ):
                    job_doc["mpid"] = mp_doc["material_id"]
        return job_docs

    def update_targets(self, items: list[dict]) -> None:
        """Update each task document with the matching MPID."""
        self.logger.info(f"Updating targets {sum(len(i) for i in items)} documents")
        for doc in chain.from_iterable(items):
            mpid = doc.get("mpid", None)
            self.job_store._collection.update_one(
                {"_id": doc["_id"]},
                {"$set": {"metadata.material_id": mpid}},
            )
