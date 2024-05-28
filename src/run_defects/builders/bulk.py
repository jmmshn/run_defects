"""Collect all the usable bulk properties data in one place."""

# %%
from __future__ import annotations

import collections
import itertools
import logging
from datetime import datetime
from typing import TYPE_CHECKING

from icecream import ic
from maggma.builders import Builder
from maggma.builders.map_builder import MapBuilder
from maggma.core.store import Store
from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.core import Structure

azlogger = logging.getLogger("azure.core.pipeline.policies.http_logging_policy")
azlogger.setLevel(logging.WARNING)

if TYPE_CHECKING:
    from collections.abc import Generator

    from jobflow import JobStore
    from maggma.stores import Store


class BulkBuilder(Builder):
    """Collect all the usable bulk locpot data in one place."""

    def __init__(
        self,
        jobstore: JobStore,
        bulk_store: Store,
        query: dict = None,
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5,
        **kwargs,
    ) -> None:
        """Init.

        Args:
            jobstore: The jobstore to get the locpot data from.
            bulk_store: The store to put the locpot data in.
            query: The query to use to get the locpot data from jobstore.
            ltol: The length tolerance for the structure matcher.
            stol: The site tolerance for the structure matcher.
            angle_tol: The angle tolerance for the structure matcher.
            kwargs: Other kwargs to pass to the parent class.
        """
        self.jobstore = jobstore
        self.bulk_store = bulk_store
        self.query = query or {}
        self.ltol, self.stol, self.angle_tol = ltol, stol, angle_tol
        super().__init__(sources=[self.jobstore], targets=[self.bulk_store], **kwargs)

    def get_items(self) -> Generator[dict, None, None]:
        """Get the items to process."""
        j_query = {}
        j_query = {"output.vasp_objects.locpot.@class": "Chgcar", **self.query}
        valid_formulas = self.jobstore.distinct(
            "output.formula_pretty", criteria=j_query
        )
        ic(j_query)
        ic(valid_formulas)
        properties = [
            "output.vasp_objects.locpot",
            "output.entry",
            "output.calcs_reversed.run_type",
            "output.calcs_reversed.task_type",
            "output.calcs_reversed.calc_type",
            "output.structure",
            "output.input",
            "output.calcs_reversed.output",
            "output.calcs_reversed.input",
        ]
        self.total = len(valid_formulas)
        for formula in valid_formulas:
            uuids_in = self.jobstore.distinct(
                "uuid", {**j_query, "output.formula_pretty": formula}
            )
            uuids_completed = self.bulk_store.distinct(
                "task_id", {**j_query, "output.formula_pretty": formula}
            )
            if set(uuids_in) == set(uuids_completed):
                continue
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
            all_calcs = list(itertools.chain.from_iterable(run_type2runs.values()))
            oldest_doc = min(all_calcs, key=lambda x: x["uuid"])
            oldest_struct = oldest_doc["structure"]
            output_doc["runs"] = dict(run_type2runs)
            output_doc["task_id"] = min(uuids)
            output_doc["uuids"] = uuids
            output_doc["structure"] = oldest_struct
            res.append(output_doc)
        return res

    def update_targets(self, items: dict | list) -> None:
        """Update the target store."""
        items = list(filter(None, itertools.chain.from_iterable(items)))
        for doc in items:
            doc[self.bulk_store.last_updated_field] = datetime.utcnow()
        self.bulk_store.update(items)


class BandGapBuilder(MapBuilder):
    """Get the band gap data from the bulk properties data."""

    def __init__(
        self, bulk_store: Store, bandgap_store: Store, query: dict = None, **kwargs
    ) -> None:
        """Init."""
        self.bulk_store = bulk_store
        self.bandgap_store = bandgap_store
        super().__init__(
            source=self.bulk_store, target=self.bandgap_store, query=query, **kwargs
        )

    def unary_function(self, item: dict) -> dict:
        """Get the band gap data from the selected document."""
        runs = item.pop("runs")
        band_gaps = dict()
        uuids = []
        for run_type, run_list in runs.items():
            best_run = min(run_list, key=_get_min_ranking_val)
            band_gaps[run_type] = {
                "uuid": best_run["uuid"],
                "vbm": best_run["calcs_reversed"][0]["output"]["vbm"],
                "cbm": best_run["calcs_reversed"][0]["output"]["cbm"],
            }
            uuids.extend([run["uuid"] for run in run_list])
        item["band_gaps"] = band_gaps
        item["uuids"] = uuids
        return item


def _get_min_ranking_val(bulk_doc: dict) -> float:
    """Quantify the quality of a bulk calculation."""
    return -bulk_doc["calcs_reversed"][0]["input"]["nkpoints"]
