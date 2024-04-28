"""Query module for the run_defects package."""

from collections.abc import Generator

from jobflow.core.store import JobStore

from run_defects.utils import JOB_STORE

RUNTYPE = {"hse06": "HSE06", "gga": {"$in": ["GGA", "GGA+U"]}}
TASKTYPE = {"relax": "Structure Optimization", "static": "Static"}


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


def get_output(
    formula: str,
    # run_type: str = None,
    task_type: str = None,
    jobstore: JobStore = None,
) -> Generator:
    """Get the output for a formula and run type."""
    query_ = {"output.formula_pretty": formula}
    # if run_type:
    #     query_["calcs_reversed.0.run_type"] = RUNTYPE[run_type]
    if task_type:
        query_["calcs_reversed.0.task_type"] = TASKTYPE[task_type]

    return (query_, jobstore)

    jobstore = jobstore or JOB_STORE

    with jobstore as store:
        yield from store.query(query_)
