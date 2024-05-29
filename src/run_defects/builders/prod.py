"""Utilities for building the production database."""

import itertools
import logging
from collections.abc import Generator

from jobflow import JobStore
from maggma.core.builder import Builder
from maggma.stores import Store

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
        finished_mpids = self.all_defects.distinct("material_id")
        job_query = {
            "output.vasp_objects.locpot.@class": "Chgcar",
            "output.vasp_objects.chgcar.@class": "Chgcar",
            "output.calcs_reversed.0.run_type": {"$in": ["GGA+U", "GGA"]},
            "output.nelements": {"$gt": 1},
            "metadata.material_id": {"$exists": 1, "$nin": finished_mpids},
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
                completed_charge_states={i: False for i in range(-4, 5)},
                auto_charge_states=charge_states,
            )
        for defect_key_, dd_ in defects_dict.items():
            yield {
                "material_id": item["material_id"],
                "bulk_formula": item["formula_pretty"],
                "run_type": item["run_type"],
                "chemsys": item["chemsys"],
                "bulk_entry": entry.as_dict(),
                "defect": dd_["defect"],
                "bulk_uuid": entry.entry_id,
                "task_id": f"{entry.entry_id}:{defect_key_}",
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
        jobstore: JobStore,
        all_defects: Store,
        **kwargs,
    ) -> None:
        """Initialize the BlessedBulkBuilder."""
        self.jobstore = jobstore
        self.all_defects = all_defects
        super().__init__(
            sources=[jobstore, all_defects], targets=[all_defects], **kwargs
        )

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
