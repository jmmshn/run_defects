"""Collect all the usable bulk locpot data in one place."""

# %%
from __future__ import annotations

import collections
import itertools
from datetime import datetime
from typing import TYPE_CHECKING

from maggma.builders import Builder
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.core import Structure

if TYPE_CHECKING:
    from collections.abc import Generator

    from jobflow import JobStore
    from maggma.stores import Store


class LocpotBuilder(Builder):
    """Collect all the usable bulk locpot data in one place."""

    def __init__(
        self,
        jobstore: JobStore,
        locpot_store: Store,
        query: dict = None,
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5,
        **kwargs,
    ) -> None:
        """Init.

        Args:
            jobstore: The jobstore to get the locpot data from.
            locpot_store: The store to put the locpot data in.
            query: The query to use to get the locpot data from jobstore.
            ltol: The length tolerance for the structure matcher.
            stol: The site tolerance for the structure matcher.
            angle_tol: The angle tolerance for the structure matcher.
            kwargs: Other kwargs to pass to the parent class.
        """
        self.jobstore = jobstore
        self.locpot_store = locpot_store
        self.query = query or {}
        self.ltol, self.stol, self.angle_tol = ltol, stol, angle_tol
        super().__init__(sources=[self.jobstore], targets=[self.locpot_store], **kwargs)

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        j_query = {}
        j_query = {"output.vasp_objects.locpot.@class": "Chgcar", **self.query}
        valid_formulas = self.jobstore.distinct(
            "output.formula_pretty", criteria=j_query
        )
        properties = [
            "output.vasp_objects.locpot",
            "output.entry",
            "output.calcs_reversed.run_type",
            "output.calcs_reversed.task_type",
            "output.calcs_reversed.calc_type",
            "output.structure",
            "output.input",
        ]
        for formula in valid_formulas:
            group = list(
                self.jobstore.query(
                    {**j_query, "output.formula_pretty": formula}, properties=properties
                )
            )
            yield {"formula": formula, "group": group}

    def process_item(self, item: dict) -> list:
        """Process the item."""
        structure_matcher = StructureMatcher(
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
            comparator=ElementComparator(),
            primitive_cell=True,
        )
        structures = []
        for d in item["group"]:
            s_ = Structure.from_dict(d["output"]["structure"])
            s_.__jobdoc = d
            structures.append(s_)
        res = []
        for g in structure_matcher.group_structures(structures):
            output_doc = {"formula_pretty": item["formula"], "runs": {}}
            run_type2runs = collections.defaultdict(list)
            jdocs = [s.__jobdoc for s in g]
            uuids = [d_["uuid"] for d_ in jdocs]
            for doc in jdocs:
                doc["output"]["uuid"] = doc["uuid"]
                doc["output"]["entry"]["entry_id"] = doc["uuid"]
                run_type2runs[doc["output"]["calcs_reversed"][0]["run_type"]].append(
                    doc["output"]
                )
            output_doc["runs"] = dict(run_type2runs)
            output_doc["task_id"] = min(uuids)
            res.append(output_doc)
        return res

    def update_targets(self, items: dict | list) -> None:
        """Update the target store."""
        items = list(filter(None, itertools.chain.from_iterable(items)))
        for doc in items:
            doc[self.locpot_store.last_updated_field] = datetime.utcnow()
        self.locpot_store.update(items)
