"""Collect all the usable bulk locpot data in one place."""

# %%
from __future__ import annotations

import collections
import datetime
import itertools
import logging
from functools import lru_cache
from typing import TYPE_CHECKING, Any

from icecream import ic
from maggma.builders import Builder
from maggma.builders.map_builder import MapBuilder
from maggma.core import Store
from monty.json import MontyDecoder
from pymatgen.analysis.defects.supercells import get_closest_sc_mat
from pymatgen.analysis.defects.thermo import DefectEntry
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.core import IStructure, Structure
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.ext.matproj import MPRester

azlogger = logging.getLogger("azure.core.pipeline.policies.http_logging_policy")
azlogger.setLevel(logging.WARNING)


def _utc() -> datetime.datetime:
    """Get the current time in UTC."""
    return datetime.datetime.now(datetime.UTC)


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


@lru_cache(maxsize=200)
def get_dielectric_data(istruct: IStructure) -> dict:
    """Get the dielectric data from MP."""
    ss = Structure.from_sites(istruct)
    with MPRester() as mp:
        mp_id = mp.find_structure(ss.remove_oxidation_states())
        (data,) = mp.materials.dielectric.search(mp_id)
    return data.model_dump()


class DielectricBuilder(MapBuilder):
    """Grab and store the dielectric data from MP."""

    def __init__(
        self, defect_entry_store: Store, dielectric_store: Store, **kwargs
    ) -> None:
        """Init."""
        self.defect_entry_store = defect_entry_store
        self.dielectric_store = dielectric_store
        super().__init__(source=defect_entry_store, target=dielectric_store, **kwargs)

    def unary_function(self, item: dict) -> dict:
        """Get the dielectric data."""
        defect_entry = DefectEntry.from_dict(item["defect_entry"])
        istruct = IStructure.from_sites(defect_entry.defect.structure)
        dielectric_data = get_dielectric_data(istruct)
        return {"task_id": item["task_id"], "dielectric_data": dielectric_data}


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
        self.match_sc_mat = match_sc_mat
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
            bulk_formula = item["formula"]
            for locpot_doc in item["locpot_docs"]:
                for run_doc in locpot_doc["runs"][defect_run_type]:
                    bulk_struct = mdecode(run_doc["structure"])
                    if structure_matcher.fit(bulk_struct, ref_struct):
                        candidate_bulk_docs.append(run_doc)

            elements = set(defect_entry.sc_entry.composition.elements)
            chemsys = "-".join(sorted([el.symbol for el in elements]))

            yield {
                "defect_entry": defect_entry.as_dict(),
                "candidate_bulk_docs": candidate_bulk_docs,
                "defect_locpot": job_doc["output"]["vasp_objects"]["locpot"],
                "task_id": defect_entry.entry_id,
                "defect_run_type": defect_run_type,
                "defect_chemsys": chemsys,
                "defect_name": defect_obj.name,
                "bulk_formula": bulk_formula,
                "defect": defect_obj.as_dict(),
            }

    def update_targets(self, items: dict | list) -> None:
        """Update the target store."""
        items = list(filter(None, itertools.chain.from_iterable(items)))
        for doc in items:
            doc[self.defect_entry_store.last_updated_field] = _utc()
        self.defect_entry_store.update(items)


class PDBuilder(Builder):
    """Grab the elemental entries and make the PD."""

    def __init__(
        self,
        jobstore: JobStore,
        defect_entry_store: Store,
        pd_store: Store,
        thermo_type: str = "GGA_GGA+U",
        **kwargs,
    ) -> None:
        """Init."""
        self.jobstore = jobstore
        self.defect_entry_store = defect_entry_store
        self.pd_store = pd_store
        self.thermo_type = thermo_type
        if self.pd_store.key != "chemsys_and_type":
            raise ValueError("PD Store key must be 'chemsys_and_type'")
        super().__init__(
            sources=[self.jobstore, self.defect_entry_store],
            targets=[self.pd_store],
            **kwargs,
        )

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        all_chemsys = self.defect_entry_store.distinct("defect_chemsys")
        for chemsys in all_chemsys:
            with MPRester() as mp:
                yield {
                    "chemsys": chemsys,
                    "phase_diagram": mp.thermo.get_phase_diagram_from_chemsys(
                        "Ga-N", "GGA_GGA+U"
                    ).as_dict(),
                }

    def process_item(self, item: dict) -> dict:
        """Process the item."""
        chemsys = item["chemsys"]
        return {
            "phase_diagram": item["phase_diagram"],
            "chemsys": chemsys,
            "thermo_type": self.thermo_type,
            "chemsys_and_type": f"{chemsys}:{self.thermo_type}",
        }

    def update_targets(self, items: dict | list) -> None:
        """Update the target store."""
        items = list(filter(None, items))
        for doc in items:
            doc[self.pd_store.last_updated_field] = _utc()
        self.pd_store.update(items)


class ElementBuilder(Builder):
    """Build the elemental entries."""

    def __init__(
        self, jobstore: JobStore, element_store: Store, query: dict = None, **kwargs
    ) -> None:
        """Init."""
        self.query = query or {}
        self.jobstore = jobstore
        self.element_store = element_store
        if self.element_store.key != "chemsys_and_type":
            raise ValueError("PD Store key must be 'chemsys_and_type'")
        super().__init__(
            sources=[self.jobstore], targets=[self.element_store], **kwargs
        )

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        unique_atoms = self.jobstore.distinct(
            "output.formula_pretty", criteria={"output.nelements": 1}
        )
        run_type2chemsys = collections.defaultdict(set)
        for chemsys in unique_atoms:
            for d in self.jobstore.query(
                {"output.formula_pretty": chemsys},
                properties=["output.calcs_reversed.run_type"],
            ):
                run_type2chemsys[d["output"]["calcs_reversed"][0]["run_type"]].add(
                    chemsys
                )
        for run_type, chemsys_set in run_type2chemsys.items():
            elements = list(set("-".join(chemsys_set).split("-")))
            js_query = {
                "output.calcs_reversed.run_type": run_type,
                "output.formula_pretty": {"$in": elements},
                **self.query,
            }
            js_properties = [
                "output.formula_pretty",
                "output.entry",
                "output.structure",
                "output.calcs_reversed.run_type",
            ]
            yield {
                "run_type": run_type,
                "job_docs": list(
                    self.jobstore.query(js_query, properties=js_properties)
                ),
            }

    def process_item(self, item: Any) -> Any:
        """Process the item."""
        job_docs = item["job_docs"]
        run_type = item["run_type"]
        for g, docs in itertools.groupby(
            job_docs, key=lambda d: d["output"]["formula_pretty"]
        ):
            entries = [*map(_get_sc_entry, docs)]
            ic(len(entries))
            stable_entry = min(entries, key=lambda ent: ent.energy_per_atom)
            yield {
                "entries": [ent_.as_dict() for ent_ in entries],
                "stable_entry": stable_entry.as_dict(),
                "chemsys_and_type": f"{g}:{run_type}",
            }

    def update_targets(self, items: dict | list) -> None:
        """Update the target store."""
        items = list(filter(None, itertools.chain.from_iterable(items)))
        for doc in items:
            doc[self.element_store.last_updated_field] = _utc()
        self.element_store.update(items)


class FreysoldtBuilder(MapBuilder):
    """Perform the Freysoldt Correction."""

    def __init__(
        self,
        defect_entry_store: Store,
        dielectric_store: Store,
        corrected_defect_entry_store: Store,
        jobstore: JobStore,
        **kwargs,
    ) -> None:
        """Init."""
        self.defect_entry_store = defect_entry_store
        self.jobstore = jobstore
        self.dielectric_store = dielectric_store
        self.corrected_defect_entry_store = corrected_defect_entry_store
        super().__init__(
            source=self.defect_entry_store,
            target=self.corrected_defect_entry_store,
            **kwargs,
        )
        self.sources.extend([self.jobstore, self.dielectric_store])

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        # call the parent class method
        for doc in super().get_items():
            doc_ = self._replace_blob(doc, dry_run=False)
            doc_["dielectric_data"] = self.dielectric_store.query_one(
                {"task_id": doc["task_id"]}
            )["dielectric_data"]
            yield doc_

    def unary_function(self, item: dict) -> dict:
        """Perform the Freysoldt Correction."""
        defect_locpot = mdecode(item["defect_locpot"])
        bulk_doc = min(item["candidate_bulk_docs"], key=_get_energy)
        bulk_locpot = mdecode(bulk_doc["vasp_objects"]["locpot"])
        defect_entry = DefectEntry.from_dict(item.pop("defect_entry"))
        dielectric = item["dielectric_data"]["e_total"]
        correction_results = defect_entry.get_freysoldt_correction(
            defect_locpot=defect_locpot,
            bulk_locpot=bulk_locpot,
            dielectric=dielectric,
        )
        defect_entry.bulk_entry = ComputedEntry.from_dict(bulk_doc["entry"])
        return {
            "freysoldt_data": correction_results.as_dict(),
            "defect_entry": defect_entry.as_dict(),
            "defect_run_type": item["defect_run_type"],
            "defect_chemsys": item["defect_chemsys"],
            "bulk_formula": item["bulk_formula"],
            "bulk_structure": bulk_doc["structure"],
            "defect_name": item["defect_name"],
            "defect": item["defect"],
            "task_id": item["task_id"],
        }

    def _replace_blob(self, doc: dict | list | Any, dry_run: bool = True) -> dict:
        """Replace the blob search docs with the data."""
        if isinstance(doc, dict):
            d_out = {}
            if "blob_uuid" in doc and "store" in doc:
                return _read_blob(
                    doc["blob_uuid"],
                    doc["store"],
                    dry_run,
                    self.jobstore.to_json(),
                )
            for k, v in doc.items():
                d_out[k] = self._replace_blob(v, dry_run=dry_run)
            return d_out
        if isinstance(doc, list):
            return [self._replace_blob(d, dry_run=dry_run) for d in doc]  # type: ignore[return-value]
        return doc


class FormationEnergyBuilder(Builder):
    """Build the formation energy data."""

    def __inti__(
        self,
        corrected_defect_entry_store: Store,
        element_store: Store,
        pd_store: Store,
        formation_energy_store: Store,
        **kwargs,
    ) -> None:
        """Init."""
        self.corrected_defect_entry_store = corrected_defect_entry_store
        self.element_store = element_store
        self.pd_store = pd_store
        self.formation_energy_store = formation_energy_store
        super().__init__(
            sources=[
                self.corrected_defect_entry_store,
                self.element_store,
                self.pd_store,
            ],
            targets=[self.formation_energy_store],
            **kwargs,
        )

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        bulk_formulas = self.corrected_defect_entry_store.distinct("bulk_formula")
        for bulk_formula in bulk_formulas:
            elements = set()
            # defect_names = self.corrected_defect_entry_store.distinct()
            for doc in self.corrected_defect_entry_store.query(
                {"bulk_formula": bulk_formula}
            ):
                # run_type = doc["defect_run_type"]
                defect_chemsys = doc["defect_chemsys"]
                elements |= set(defect_chemsys.split("-"))
                yield doc


@lru_cache(maxsize=200)
def _read_blob(blob_uuid: str, store: str, dry_run: bool, jobstore_str: str) -> dict:
    jobstore = MontyDecoder().decode(jobstore_str)
    blob_store = jobstore.additional_stores[store]
    blob_or_index_store = blob_store.index if dry_run else blob_store
    with blob_or_index_store as s:
        blob_data = s.query_one({"blob_uuid": blob_uuid})
    return blob_data if dry_run else blob_data["data"]


def _get_energy(bulk_doc: dict) -> float:
    """Get the effective energy for each bulk doc."""
    return ComputedEntry.from_dict(bulk_doc["entry"]).energy_per_atom


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
