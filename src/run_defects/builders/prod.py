"""Utilities for building the production database."""

import itertools
import logging
from collections.abc import Generator

from jobflow import JobStore
from maggma.core.builder import Builder
from maggma.stores import Store
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher

from run_defects.submit import get_defects
from run_defects.utils import mdecode

azlogger = logging.getLogger("azure.core.pipeline.policies.http_logging_policy")
azlogger.setLevel(logging.WARNING)


class AllDefectsBuilder(Builder):
    """Check the finished jobstore and create a library of defects to compute."""

    def __init__(
        self,
        jobstore: JobStore,
        all_defects: Store,
        query: dict = None,
        **kwargs,
    ) -> None:
        """Initialize the BlessedBulkBuilder."""
        self.jobstore = jobstore
        self.all_defects = all_defects
        self.query = query or {}
        super().__init__(sources=[jobstore], targets=[all_defects], **kwargs)

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        job_query = {
            "output.vasp_objects.locpot.@class": "Chgcar",
            "output.vasp_objects.chgcar.@class": "Chgcar",
            "output.calcs_reversed.0.run_type": {"$in": ["GGA+U", "GGA"]},
            "output.nelements": {"$gt": 1},
            **self.query,
        }
        self.logger.info(f"QUERY: {job_query}")
        agg_pipe = self.jobstore._collection.aggregate(
            [
                {
                    "$match": job_query,
                },
                {
                    "$group": {
                        "_id": {
                            "material_id": "$metadata.material_id",
                            "formula_pretty": "$output.formula_pretty",
                            "run_type": "$output.calcs_reversed.0.run_type",
                            "chemsys": "$output.chemsys",
                        },
                        "uuid2entry": {
                            "$push": {
                                "uuid": "$uuid",
                                "entry": "$output.entry",
                                "chgcar_blob_uuid": "$output"
                                ".vasp_objects.chgcar.blob_uuid",
                            }
                        },
                    }
                },
            ]
        )
        agg_results = list(agg_pipe)
        for group in agg_results:
            best_entry = None
            chgcar_blob_uuid = None
            for u2e in mdecode(group["uuid2entry"]):
                ent_ = u2e["entry"]
                if (
                    best_entry is None
                    or ent_.energy_per_atom < best_entry.energy_per_atom
                ):
                    ent_.entry_id = u2e["uuid"]
                    best_entry = ent_
                    chgcar_blob_uuid = u2e["chgcar_blob_uuid"]
            yield {
                "chgcar": self.jobstore.additional_stores["data"].query_one(
                    criteria={"blob_uuid": chgcar_blob_uuid}
                ),
                "entry": best_entry,
                **group["_id"],
            }

    def process_item(self, item: dict) -> Generator[dict, None, None]:
        """Process an item."""
        entry = item["entry"]
        chgcar = mdecode(item["chgcar"]["data"])
        defects_dict = dict()
        for defect, ii in get_defects(chgcar=chgcar, max_iter=3):
            charge_states = defect.get_charge_states()
            defects_dict[f"{defect.name}:{ii}"] = dict(
                defect=defect.as_dict(),
                completed_charge_states={
                    "GGA_GGA+U": {i: None for i in range(-6, 7)},
                    "HSE06": {i: None for i in range(-6, 7)},
                },
                auto_charge_states=charge_states,
                defect_name=defect.name,
            )

        for defect_key_, dd_ in defects_dict.items():
            yield {
                "material_id": item["material_id"],
                "bulk_formula": item["formula_pretty"],
                "run_type": item["run_type"],
                "chemsys": item["chemsys"],
                "bulk_entry": entry.as_dict(),
                "bulk_uuid": entry.entry_id,
                "task_id": f"{item['material_id']}:{defect_key_}",
                **dd_,
            }

    def update_targets(self, items: list) -> None:
        """Update the target store."""
        items = list(filter(None, itertools.chain.from_iterable(items)))
        self.logger.info(f"Updating {len(items)} documents")
        self.all_defects.update(docs=items)


class TagFinishedDefect(Builder):
    """Tag the finished defect as completed."""

    def __init__(
        self,
        defect_entries: JobStore,
        all_defects: Store,
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5,
        query: dict = None,
        **kwargs,
    ) -> None:
        """Initialize the BlessedBulkBuilder."""
        self.defect_entries = defect_entries
        self.all_defects = all_defects
        self.ltol, self.stol, self.angle_tol = ltol, stol, angle_tol
        self.query = query or {}
        super().__init__(sources=[defect_entries, all_defects, all_defects], **kwargs)

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        all_defects_query = self.query
        formula_and_defect_name_pipe = self.all_defects._collection.aggregate(
            [
                {
                    "$match": all_defects_query,
                },
                {
                    "$group": {
                        "_id": {
                            "bulk_formula": "$bulk_formula",
                            "defect_name": "$defect_name",
                        },
                        "defect_docs": {
                            "$push": {
                                "material_id": "$material_id",
                                "bulk_uuid": "$bulk_uuid",
                                "task_id": "$task_id",
                                "defect": "$defect",
                            }
                        },
                    }
                },
            ]
        )
        formula_and_defect_name_results = list(formula_and_defect_name_pipe)
        self.logger.info(f"Found {len(formula_and_defect_name_results)} groups")
        for group in formula_and_defect_name_results:
            bulk_formula = group["_id"]["bulk_formula"]
            defect_name = group["_id"]["defect_name"]
            defect_docs = mdecode(group["defect_docs"])
            for de_doc in self.defect_entries.query(
                criteria={"bulk_formula": bulk_formula, "defect_name": defect_name}
            ):
                yield {"de_doc": de_doc, "defect_docs": defect_docs}

    def process_item(self, item: dict) -> Generator[list, None, None]:
        """Process an item."""
        sm = StructureMatcher(
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
            comparator=ElementComparator(),
        )
        de_doc = item["de_doc"]
        qq_ = de_doc["defect_entry"]["charge_state"]
        defect_docs = item["defect_docs"]
        run_type = de_doc["defect_run_type"]
        if run_type == "GGA" or run_type == "GGA+U":
            run_type = "GGA_GGA+U"
        dent = mdecode(de_doc["defect_entry"])
        dent_defect_structure = dent.defect.defect_structure
        dent_bulk_structure = dent.defect.structure
        matched_ids = []
        for defect_doc in defect_docs:
            defect = defect_doc["defect"]
            defect_structure = defect.defect_structure
            bulk_structure = defect.structure
            if sm.fit(dent_defect_structure, defect_structure) and sm.fit(
                dent_bulk_structure, bulk_structure
            ):
                yield [
                    {"task_id": defect_doc["task_id"]},
                    {
                        "$set": {
                            f"completed_charge_states.{run_type}.{qq_}": de_doc[
                                "task_id"
                            ]
                        }
                    },
                ]
                matched_ids.append(defect_doc["task_id"])
                if len(matched_ids) > 1:
                    self.logger.warning(
                        f"Multiple matching defect structures, {matched_ids}"
                    )

    def update_targets(self, items: list) -> None:
        """Update the target store."""
        updates = list(filter(None, itertools.chain.from_iterable(items)))
        self.logger.info(f"Updating {len(items)} documents")
        for up in updates:
            self.all_defects._collection.update_one(*up)
