"""Collect all the usable bulk locpot data in one place."""

from __future__ import annotations

import collections
import datetime
import itertools
import logging
from functools import lru_cache
from typing import TYPE_CHECKING, Any

from maggma.builders import Builder
from maggma.builders.map_builder import MapBuilder
from maggma.core import Store
from monty.json import MontyDecoder
from pymatgen.analysis.defects.supercells import get_closest_sc_mat
from pymatgen.analysis.defects.thermo import (
    DefectEntry,
    FormationEnergyDiagram,
    group_defect_entries,
)
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.entries.computed_entries import Composition, ComputedEntry
from pymatgen.ext.matproj import MPRester

from run_defects.utils import jdoc_to_entry, mdecode

azlogger = logging.getLogger("azure.core.pipeline.policies.http_logging_policy")
azlogger.setLevel(logging.WARNING)
pmglogger = logging.getLogger("pymatgen.core.structure")
pmglogger.setLevel(logging.ERROR)


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
    "metadata",
]


class DielectricBuilder(MapBuilder):
    """Grab and store the dielectric data from MP."""

    def __init__(
        self,
        defect_entry_store: Store,
        jobstore: JobStore,
        dielectric_store: Store,
        **kwargs,
    ) -> None:
        """Init."""
        self.defect_entry_store = defect_entry_store
        self.dielectric_store = dielectric_store
        self.jobstore = jobstore
        super().__init__(source=defect_entry_store, target=dielectric_store, **kwargs)
        self.sources.append(jobstore)

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        for doc in super().get_items():
            bulk_uuid = doc["candidate_bulk_docs"][0]["uuid"]
            mp_doc_ = self.jobstore.query_one(
                {"uuid": bulk_uuid}, properties=["metadata"]
            )
            if mp_doc_ is None:
                continue
            mp_id = mp_doc_["metadata"]["material_id"]
            de_data = _get_dielectric_data(mp_id)
            if de_data is None:
                continue
            yield {
                "task_id": doc["task_id"],
                "mp_id": mp_id,
                "dielectric_data": _get_dielectric_data(mp_id),
            }

    def unary_function(self, item: dict) -> dict:
        """Get the dielectric data."""
        return item


class DefectEntryBuilder(Builder):
    """Collect the Defect Entry Data."""

    def __init__(
        self,
        jobstore: JobStore,
        bulk_store: Store,
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
            jobstore: The jobstore.
            bulk_store: The store to put the bulk data in.
            defect_entry_store: The store to put the defect entry data in.
            query: The query to use to get the defect ent data from jobstore.
            ltol: The length tolerance for the structure matcher.
            stol: The site tolerance for the structure matcher.
            angle_tol: The angle tolerance for the structure matcher.
            match_sc_mat: Whether to compare only with bulk SC structures.
            kwargs: Other kwargs to pass to the parent class.
        """
        self.jobstore = jobstore
        self.bulk_store = bulk_store
        self.defect_entry_store = defect_entry_store
        self.query = query or {}
        self.match_sc_mat = match_sc_mat
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
            sources=[self.jobstore, self.bulk_store],
            targets=[self.defect_entry_store],
            **kwargs,
        )

    @property
    def _primitive_cell(self) -> bool:
        return not self.match_sc_mat

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        j_query = {**DEFECT_JOB_QUERY, **self.query}
        # get groups of uuids with the same bulk formula and defect name
        agg_pipe = self.jobstore.docs_store._collection.aggregate(
            [
                {"$match": j_query},
                {
                    "$group": {
                        "_id": {
                            "name": "$output.additional_json.info.defect_name",
                            "bulk_formula": "$output.additional_json.info.bulk_formula",
                        },
                        "uuids": {"$addToSet": "$uuid"},
                    }
                },
            ]
        )
        finished_task_ids = set(self.defect_entry_store.distinct("task_id"))
        agg_results = list(agg_pipe)
        self.total = len(agg_results)
        self.logger.info(
            f"Found {self.total} groups (defect_name + bulk_formula) to process."
        )
        for group in agg_results:
            # `_id` contains keys `name` and `bulk`
            # `uuids` contains a list of uuids
            if finished_task_ids.issuperset(group["uuids"]):
                continue

            missing_uuids = list(set(group["uuids"]) - finished_task_ids)
            self.logger.info(
                f"Processing {len(missing_uuids)}) missing uuids for group: "
                "{group['_id']}"
            )
            formula = group["_id"]["bulk_formula"]
            defect_name = group["_id"]["name"]
            job_docs = list(
                self.jobstore.query(
                    {"uuid": {"$in": missing_uuids}},
                    properties=DEFECT_JOB_PROPERTIES,
                )
            )
            bulk_query = {
                "formula_pretty": formula,
            }
            bulk_docs = list(self.bulk_store.query(bulk_query))
            yield {
                "formula": formula,
                "defect_name": defect_name,
                "job_docs": job_docs,
                "bulk_docs": bulk_docs,
            }

        # import ipdb; ipdb.set_trace()
        # for d in self.jobstore.query(
        #     j_query,
        #     properties=[
        #         "output.additional_json.info.defect_name",
        #         "output.additional_json.info.bulk_formula",
        #         "output.calcs_reversed.run_type",
        #     ],
        # ):
        #     run_type_ = d["output"]["calcs_reversed"][0]["run_type"]
        #     bulk_formula_ = d["output"]["additional_json"]["info"]["bulk_formula"]
        #     defect_name_ = d["output"]["additional_json"]["info"]["defect_name"]
        #     defect_run_combos.add((bulk_formula_, defect_name_, run_type_))
        # finished_task_ids = self.defect_entry_store.distinct("task_id")
        # for formula, defect_name, run_type in defect_run_combos:
        #     combo_query = {
        #         "output.additional_json.info.bulk_formula": formula,
        #         "output.additional_json.info.defect_name": defect_name,
        #         "output.calcs_reversed.run_type": run_type,
        #     }
        #     job_docs = list(
        #         self.jobstore.query(
        #             combo_query | {"uuid": {"$nin": finished_task_ids}},
        #             properties=DEFECT_JOB_PROPERTIES,
        #         )
        #     )
        #     bulk_query = {
        #         "formula_pretty": formula,
        #         f"runs.{run_type}": {"$exists": True},
        #     }
        #     bulk_docs = list(self.bulk_store.query(bulk_query))
        #     yield {
        #         "formula": formula,
        #         "defect_name": defect_name,
        #         "run_type": run_type,
        #         "job_docs": job_docs,
        #         "locpot_docs": bulk_docs,
        #     }

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
            if self.match_sc_mat:
                sc_mat = get_closest_sc_mat(
                    uc_struct=defect_obj.structure, sc_struct=defect_sc_struct
                )
                ref_struct = ref_struct * sc_mat

            defect_entry = DefectEntry(
                defect=defect_obj,
                charge_state=charge_state,
                sc_entry=jdoc_to_entry(job_doc, inc_structure=True),
            )

            if int(defect_entry.charge_state) != int(
                defect_entry.sc_entry.structure._charge
            ):
                raise ValueError(
                    f"Defect OBJ charge state ({int(defect_entry.charge_state)}) "
                    "does not match structure charge "
                    f"({int(defect_entry.sc_entry.structure._charge)})\n"
                    f"UUID: {job_doc['uuid']}\n"
                )
            candidate_bulk_docs = []
            bulk_formula = item["formula"]
            for bulk_doc in item["bulk_docs"]:
                for run_doc in bulk_doc["runs"][defect_run_type]:
                    bulk_struct = mdecode(run_doc["structure"])
                    if structure_matcher.fit(bulk_struct, ref_struct):
                        candidate_bulk_docs.append(run_doc)
            bulk_uuids = [d["uuid"] for d in candidate_bulk_docs]
            elements = set(defect_entry.sc_entry.composition.elements)
            chemsys = "-".join(sorted([el.symbol for el in elements]))
            if candidate_bulk_docs:
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
                    "bulk_uuids": bulk_uuids,
                }
            else:
                self.logger.warning(
                    f"No bulk structures found for {defect_obj.name} "
                    "in {bulk_formula}.\n"
                    f"Defect Run uuid: {job_doc['uuid']}.\n"
                    f"Defect Run metadata: {job_doc['metadata']}.\n"
                    f"defect_run_type: {defect_run_type}.\n"
                )

    def update_targets(self, items: dict | list) -> None:
        """Update the target store."""
        items = list(filter(None, itertools.chain.from_iterable(items)))
        self.logger.info(f"Updating {len(items)} defect entries.")
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
        all_chemsys_and_type = self.pd_store.distinct("chemsys_and_type")
        for chemsys in all_chemsys:
            if f"{chemsys}:{self.thermo_type}" in all_chemsys_and_type:
                self.logger.info(f"Skipping {chemsys} for {self.thermo_type}")
                continue
            with MPRester() as mp:
                try:
                    pd_entries = mp.get_entries_in_chemsys(
                        chemsys,
                        additional_criteria={
                            "energy_above_hull": (0.0, 0.1),
                            "thermo_types": [self.thermo_type],
                        },
                    )
                    pd_mp = PhaseDiagram(pd_entries)
                except Exception as e:
                    self.logger.warning(f"Error getting PD for {chemsys}: {e}")
                    continue
                yield {
                    "chemsys": chemsys,
                    "phase_diagram": pd_mp.as_dict(),
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
        self.logger.info(f"Updating {len(items)} phase diagrams.")
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
            entries = [*map(jdoc_to_entry, docs)]
            stable_entry = min(entries, key=lambda ent: ent.energy_per_atom)
            comp_ = stable_entry.composition
            chemsys = "-".join(sorted([el.symbol for el in comp_.elements]))
            yield {
                "entries": [ent_.as_dict() for ent_ in entries],
                "stable_entry": stable_entry.as_dict(),
                "chemsys_and_type": f"{g}:{run_type}",
                "run_type": run_type,
                "chemsys": chemsys,
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
            dielectric_doc = self.dielectric_store.query_one(
                {"task_id": doc["task_id"]}
            )
            if dielectric_doc is None:
                self.logger.warning(
                    f"Dielectric data not found for {doc['task_id']}. Skipping."
                )
                continue
            doc_["dielectric_data"] = dielectric_doc["dielectric_data"]
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
        sc_mat = get_closest_sc_mat(
            uc_struct=bulk_locpot.structure, sc_struct=defect_locpot.structure
        )
        bulk_sc_composition_ = (bulk_locpot.structure * sc_mat).composition
        bulk_doc_composition_ = Composition(bulk_doc["entry"]["composition"])

        bulk_sc_fu, bulk_sc_factor = (
            bulk_sc_composition_.get_reduced_composition_and_factor()
        )
        bulk_doc_fu, bulk_doc_factor = (
            bulk_doc_composition_.get_reduced_composition_and_factor()
        )
        if bulk_sc_fu != bulk_doc_fu:
            raise ValueError(f"bulk_sc_fu: {bulk_sc_fu}, bulk_doc_fu: {bulk_doc_fu}")

        self.logger.info(
            f"SC COMPOSITION: {bulk_sc_composition_}, FACTOR: {bulk_sc_factor}"
            f"UC COMPOSITION: {bulk_doc_composition_}, FACTOR: {bulk_doc_factor}"
        )
        bulk_entry_dict = bulk_doc["entry"]
        bulk_entry_dict["composition"] = bulk_sc_composition_
        bulk_entry_dict["energy"] = bulk_doc["entry"]["energy"] * (
            bulk_sc_factor / bulk_doc_factor
        )
        defect_entry.bulk_entry = ComputedEntry.from_dict(bulk_entry_dict)

        return {
            "freysoldt_data": correction_results.as_dict(),
            "defect_entry": defect_entry.as_dict(),
            "defect_run_type": item["defect_run_type"],
            "defect_chemsys": item["defect_chemsys"],
            "bulk_formula": item["bulk_formula"],
            "bulk_structure": bulk_doc["structure"],
            "bulk_entry": bulk_doc["entry"],
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

    def __init__(
        self,
        corrected_defect_entry_store: Store,
        element_store: Store,
        pd_store: Store,
        bandgap_store: Store,
        formation_energy_store: Store,
        thermo_type: str = "GGA_GGA+U",
        query: dict = None,
        **kwargs,
    ) -> None:
        """Init."""
        self.corrected_defect_entry_store = corrected_defect_entry_store
        self.element_store = element_store
        self.pd_store = pd_store
        self.bandgap_store = bandgap_store
        self.formation_energy_store = formation_energy_store
        self.query = query or {}
        self.thermo_type = thermo_type
        super().__init__(
            sources=[
                self.corrected_defect_entry_store,
                self.element_store,
                self.bandgap_store,
                self.pd_store,
            ],
            targets=[self.formation_energy_store],
            **kwargs,
        )

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        agg_pipe = self.corrected_defect_entry_store._collection.aggregate(
            [
                {
                    "$match": {
                        "defect_name": {"$exists": 1},
                        "bulk_formula": {"$exists": 1},
                        "defect_run_type": {"$exists": 1},
                        **self.query,
                    }
                },
                {
                    "$group": {
                        "_id": {
                            "defect_run_type": "$defect_run_type",
                            "bulk_formula": "$bulk_formula",
                            "defect_name": "$defect_name",
                        },
                        "task_ids": {"$addToSet": "$task_id"},
                    }
                },
            ]
        )
        agg_results = list(agg_pipe)
        self.total = len(agg_results)
        self.logger.info(
            f"Found {self.total} groups (defect_run_type + bulk_formula) to process."
        )
        for group in agg_results:
            bulk_formula = group["_id"]["bulk_formula"]
            defect_name = group["_id"]["defect_name"]
            defect_uuids = group["task_ids"]
            run_type = group["_id"]["defect_run_type"]
            completed_uuids = self.formation_energy_store.distinct(
                "task_id", {"task_id": {"$in": defect_uuids}}
            )
            if set(defect_uuids) == set(completed_uuids):
                continue
            self.logger.info(
                f"Getting ({len(defect_uuids)}) uuids for group: {group['_id']}"
            )
            defect_entry_docs = list(
                self.corrected_defect_entry_store.query(
                    {"task_id": {"$in": defect_uuids}}
                )
            )
            defect_chemsys = defect_entry_docs[0]["defect_chemsys"]
            elements = list(set(defect_chemsys.split("-")))
            element_docs = list(
                self.element_store.query(
                    {"chemsys": {"$in": elements}, "run_type": run_type}
                )
            )
            pd_doc = self.pd_store.query_one(
                {"chemsys": defect_chemsys, "thermo_type": self.thermo_type}
            )
            if pd_doc is None:
                self.logger.error(
                    f"Phase diagram data not found for {defect_chemsys}. Skipping."
                )
                continue
            bandgap_docs = list(
                self.bandgap_store.query({"formula_pretty": bulk_formula})
            )
            yield {
                "bulk_formula": bulk_formula,
                "defect_name": defect_name,
                "defect_chemsys": defect_chemsys,
                "defect_entry_docs": defect_entry_docs,
                "element_docs": element_docs,
                "band_gap_docs": bandgap_docs,
                "pd_doc": pd_doc,
                "run_type": run_type,
            }

    def process_item(self, item: dict) -> Generator:
        """Process the item."""
        pd_entries = mdecode(item["pd_doc"]["phase_diagram"]["all_entries"])
        elements = mdecode([d["stable_entry"] for d in item["element_docs"]])
        defect_entries = mdecode([d["defect_entry"] for d in item["defect_entry_docs"]])
        bandgap_data = {
            tuple(d["uuids"]): d["band_gaps"][item["run_type"]]
            for d in item["band_gap_docs"]
        }

        for _defect_name, gdents_ in group_defect_entries(defect_entries):
            bulk_uuid = gdents_[0].bulk_entry.entry_id
            bg_key = next(k for k in bandgap_data if bulk_uuid in k)
            vbm, cbm = bandgap_data[bg_key]["vbm"], bandgap_data[bg_key]["cbm"]
            if vbm is None or cbm is None:
                self.logger.warning(
                    f"Band gap data not found for {bulk_uuid}. Skipping."
                )
                continue
            pd_ = PhaseDiagram(pd_entries)
            self.logger.info(
                f"Elements: {[e.composition.reduced_formula for e in elements]},"
                " VBM: {vbm}, CBM: {cbm}"
            )
            fed = FormationEnergyDiagram.with_atomic_entries(
                defect_entries=gdents_,
                atomic_entries=elements,
                phase_diagram=pd_,
                vbm=vbm,
                band_gap=cbm - vbm,
            )
            defect_uuids = [d.entry_id for d in gdents_]
            yield {
                "bulk_uuid": bulk_uuid,
                "defect_uuids": defect_uuids,
                "task_id": min(defect_uuids),
                "vbm": vbm,
                "cbm": cbm,
                "bulk_formula": item["bulk_formula"],
                "defect_chemsys": item["defect_chemsys"],
                "defect_name": _defect_name,
                "fed": fed.as_dict(),
                "run_type": item["run_type"],
            }

    def update_targets(self, items: dict | list) -> None:
        """Update the target store."""
        items = list(filter(None, itertools.chain.from_iterable(items)))
        for doc in items:
            doc[self.element_store.last_updated_field] = _utc()
        self.formation_energy_store.update(items)


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


# def jdoc_to_entry(doc: dict) -> ComputedStructureEntry:
#     """Get the ComputedStructureEntry."""
#     entry_dict = doc["output"]["entry"]
#     structure = mdecode(doc["output"]["structure"])
#     entry_id = doc["uuid"]
#     return ComputedStructureEntry.from_dict(
#         {**entry_dict, "structure": structure, "entry_id": entry_id}
#     )


@lru_cache(maxsize=200)
def _get_dielectric_data(mp_id: str) -> dict:
    """Get the dielectric data from MP."""
    with MPRester() as mp:
        res = mp.materials.dielectric.search(mp_id)
    if not res:
        return None
    return res[0].model_dump()
