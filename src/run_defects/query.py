"""Query module for the run_defects package."""

from collections.abc import Generator

from jobflow.core.store import JobStore

from run_defects.utils import JOB_STORE


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


# def _has_dielectric(struct):
#     with MPRester() as mpr:
#         data = mpr.get_dielectric_data(struct)
