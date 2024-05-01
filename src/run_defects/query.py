"""Query module for the run_defects package."""

from __future__ import annotations

from typing import TYPE_CHECKING

from run_defects.utils import JOB_STORE, jdoc_to_entry

if TYPE_CHECKING:
    from collections.abc import Generator

    from jobflow.core.store import JobStore
    from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry

RUNTYPE = {"hse06": "HSE06", "gga": {"$in": ["GGA", "GGA+U"]}}
TASKTYPE = {
    "relax": "Structure Optimization",
    "static": "Static",
    "deformation": "Deformation",
}


def get_structure_with_volumetric_data(
    query: dict = None, jobstore: JobStore = None
) -> Generator:
    """Get the structure with volumetric data."""
    query = query or {}
    js_query = {
        "output.vasp_objects.locpot.@class": {"$exists": True},
        "output.vasp_objects.chgcar.@class": {"$exists": True},
        **query,
    }
    jobstore = jobstore or JOB_STORE
    properties = [
        "output.structure",
        "output.vasp_objects.locpot",
        "output.vasp_objects.chgcar",
        "metadata",
    ]
    yield from jobstore.query(js_query, properties=properties)


def get_outputs(
    formula: str,
    run_type: str = None,
    task_type: str = None,
    jobstore: JobStore = None,
    query: dict = None,
) -> Generator[dict, None, None]:
    """Get the output for a formula and run type."""
    query_ = query or {}
    query_["output.formula_pretty"] = formula
    if run_type:
        query_["output.calcs_reversed.0.run_type"] = RUNTYPE[run_type]
    if task_type:
        query_["output.calcs_reversed.0.task_type"] = TASKTYPE[task_type]

    jobstore = jobstore or JOB_STORE
    with jobstore as store:
        yield from store.query(query_)


def get_entries(
    formula: str,
    run_type: str = None,
    task_type: str = None,
    jobstore: JobStore = None,
    inc_structures: bool = False,
) -> list[ComputedEntry | ComputedStructureEntry]:
    """Get the output for a formula and run type."""
    return [
        jdoc_to_entry(doc_, inc_structure=inc_structures)
        for doc_ in get_outputs(
            formula=formula, run_type=run_type, task_type=task_type, jobstore=jobstore
        )
    ]


def get_min_energy_entry(
    formula: str,
    run_type: str = None,
    task_type: str = None,
    jobstore: JobStore = None,
    inc_structures: bool = False,
) -> ComputedEntry | ComputedStructureEntry:
    """Get the minimum energy entry for a formula and run type."""
    entries = get_entries(
        formula=formula,
        run_type=run_type,
        task_type=task_type,
        jobstore=jobstore,
        inc_structures=inc_structures,
    )
    return min(entries, key=lambda x: x.energy_per_atom)
