"""Collect all the usable bulk locpot data in one place."""

# %%
from __future__ import annotations

import itertools
from datetime import datetime
from typing import TYPE_CHECKING, Any

from icecream import ic
from maggma.builders import Builder
from maggma.builders.map_builder import MapBuilder
from monty.json import MontyDecoder
from pymatgen.analysis.defects.supercells import get_closest_sc_mat
from pymatgen.analysis.defects.thermo import DefectEntry
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.entries.computed_entries import ComputedStructureEntry

if TYPE_CHECKING:
    from collections.abc import Generator

    from jobflow import JobStore
    from maggma.stores import Store

DEFECT_JOB_QUERY = {
    "output.additional_json.info.defect_name": {"$exists": True},
}

DEFECT_JOB_PROPERTIES = [
    "output.structure",
    "output.vasp_objects.locpot",
    "output.entry",
    "output.calcs_reversed.run_type",
    "output.additional_json.info",
]


class DefectEntryBuilder(Builder):
    """Collect the Defect Entry Data."""

    def __init__(
        self,
        jobstore: JobStore,
        locpot_store: Store,
        defect_entry_store: Store,
        query: dict = None,
        stol: float = 0.3,
        ltol: float = 0.2,
        angle_tol: float = 5,
        match_sc_mat: bool = False,
        **kwargs,
    ) -> None:
        """Init.

        Args:
            jobstore: The jobstore to get the locpot data from.
            locpot_store: The store to put the locpot data in.
            defect_entry_store: The store to put the defect entry data in.
            query: The query to use to get the locpot data from jobstore.
            ltol: The length tolerance for the structure matcher.
            stol: The site tolerance for the structure matcher.
            angle_tol: The angle tolerance for the structure matcher.
            primitive_cell: Whether to use the primitive cell for the structure matcher.
            match_sc_mat: Whether to compare only with bulk SC structures.
            kwargs: Other kwargs to pass to the parent class.
        """
        self.jobstore = jobstore
        self.locpot_store = locpot_store
        self.defect_entry_store = defect_entry_store
        self.query = query or {}
        self.get_bulk_sc = match_sc_mat
        self._primitive_cell = not self.get_bulk_sc
        (
            self.ltol,
            self.stol,
            self.angle_tol,
        ) = (
            ltol,
            stol,
            angle_tol,
        )
        super().__init__(
            sources=[self.jobstore, self.locpot_store],
            targets=[self.defect_entry_store],
            **kwargs,
        )

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        j_query = {**DEFECT_JOB_QUERY, **self.query}

        defect_run_combos = set()
        for d in self.jobstore.query(
            j_query,
            properties=[
                "output.additional_json.info.defect_name",
                "output.additional_json.info.bulk_formula",
                "output.calcs_reversed.run_type",
            ],
        ):
            run_type_ = d["output"]["calcs_reversed"][0]["run_type"]
            bulk_formula_ = d["output"]["additional_json"]["info"]["bulk_formula"]
            defect_name_ = d["output"]["additional_json"]["info"]["defect_name"]
            defect_run_combos.add((bulk_formula_, defect_name_, run_type_))

        for formula, defect_name, run_type in defect_run_combos:
            combo_query = {
                "output.additional_json.info.bulk_formula": formula,
                "output.additional_json.info.defect_name": defect_name,
                "output.calcs_reversed.run_type": run_type,
            }
            ic(combo_query)
            job_docs = list(
                self.jobstore.query(combo_query, properties=DEFECT_JOB_PROPERTIES)
            )
            locpot_query = {
                "formula_pretty": formula,
                f"runs.{run_type}": {"$exists": True},
            }
            locpot_docs = list(self.locpot_store.query(locpot_query))
            yield {
                "formula": formula,
                "defect_name": defect_name,
                "run_type": run_type,
                "job_docs": job_docs,
                "locpot_docs": locpot_docs,
            }

    def process_item(self, item: dict) -> Generator:
        """Process the item."""
        structure_matcher = StructureMatcher(
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
            comparator=ElementComparator(),
            primitive_cell=self._primitive_cell,
        )

        for job_doc in item["job_docs"]:
            defect_obj = mdecode(job_doc["output"]["additional_json"]["info"]["defect"])
            charge_state = job_doc["output"]["additional_json"]["info"]["charge_state"]
            defect_run_type = job_doc["output"]["calcs_reversed"][0]["run_type"]
            defect_sc_struct = mdecode(job_doc["output"]["structure"])
            ref_struct = defect_obj.structure
            if self.get_bulk_sc:
                sc_mat = get_closest_sc_mat(
                    uc_struct=defect_obj.structure, sc_struct=defect_sc_struct
                )
                ref_struct = ref_struct * sc_mat

            defect_entry = DefectEntry(
                defect=defect_obj,
                charge_state=charge_state,
                sc_entry=_get_sc_entry(job_doc),
            )
            candidate_bulk_docs = []

            for locpot_doc in item["locpot_docs"]:
                for run_doc in locpot_doc["runs"][defect_run_type]:
                    bulk_struct = mdecode(run_doc["structure"])
                    if structure_matcher.fit(bulk_struct, ref_struct):
                        candidate_bulk_docs.append(run_doc)

            yield {
                "defect_entry": defect_entry.as_dict(),
                "candidate_bulk_docs": candidate_bulk_docs,
                "defect_locpot": job_doc["output"]["vasp_objects"]["locpot"],
                "task_id": defect_entry.entry_id,
            }

    def update_targets(self, items: dict | list) -> None:
        """Update the target store."""
        items = list(filter(None, itertools.chain.from_iterable(items)))
        for doc in items:
            doc[self.defect_entry_store.last_updated_field] = datetime.utcnow()
        self.defect_entry_store.update(items)


class FreysoldtBuilder(MapBuilder):
    """Perform the Freysoldt Correction."""

    def __init__(
        self,
        defect_entry_store: Store,
        corrected_defect_entry_store: Store,
        **kwargs,
    ) -> None:
        """Init."""
        self.defect_entry_store = defect_entry_store
        self.corrected_defect_entry_store = corrected_defect_entry_store
        super().__init__(
            source=self.defect_entry_store,
            target=self.corrected_defect_entry_store,
            **kwargs,
        )


def mdecode(data: dict) -> Any:
    """Decode the data."""
    return MontyDecoder().process_decoded(data)


def _get_sc_entry(doc: dict) -> ComputedStructureEntry:
    """Get the ComputedStructureEntry."""
    entry_dict = doc["output"]["entry"]
    structure = mdecode(doc["output"]["structure"])
    entry_id = doc["uuid"]
    return ComputedStructureEntry.from_dict(
        {**entry_dict, "structure": structure, "entry_id": entry_id}
    )


# %%
