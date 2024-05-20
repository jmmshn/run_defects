"""Utility functions for the defects module."""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING, Any

import jobflow
from fireworks import LaunchPad
from monty.json import MontyDecoder
from monty.serialization import loadfn
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry

if TYPE_CHECKING:
    from pymatgen.core import Structure

logger = logging.getLogger(__name__)


class ComboMaker(jobflow.Maker):
    """Combine multiple makers into a single maker."""

    def __init__(self, makers: jobflow.Maker) -> None:
        """Init."""
        self.makers = makers

    def make(self, structure: Structure) -> jobflow.Job:
        """Make the job."""
        prv_struct = structure
        prev_dir = None
        jobs = []
        for maker in self.makers:
            job = maker.make(prv_struct, prev_dir=prev_dir)
            jobs.append(job)
            prv_struct = job.output.structure
            prev_dir = job.output.dir_name
        return jobflow.Flow(jobs, output=jobs[-1].output)


def mdecode(data: dict | list) -> Any:
    """Decode the data."""
    return MontyDecoder().process_decoded(data)


def rget(d: dict, keys: str) -> Any:
    """Parse keys separated by dots.

    Args:
        d: dictionary
        keys: keys separated by dots

    Return:
        value of the key
    """
    if "." in keys:
        key, rest = keys.split(".", 1)
        if key.isdigit():
            key_ = int(key)
            return rget(d[key_], rest)
        if isinstance(d, list):
            key = 0
            rest = keys
        return rget(d[key], rest)
    return d[keys]


defect_lpad_path = os.environ.get("DEFECT_LAUNCHPAD", None)
if defect_lpad_path is None:
    logger.info("DEFECT_LAUNCHPAD not set. Using default LaunchPad.")
    LPAD = LaunchPad.auto_load()
else:
    logger.info(f"DEFECT_LAUNCHPAD set to {defect_lpad_path}")
    LPAD = LaunchPad.from_file(defect_lpad_path)

defect_jobstore_path = os.environ.get("DEFECT_JOBSTORE", None)
if defect_jobstore_path is None:
    logger.info("DEFECT_JOBSTORE not set. Using default jobstore.")
    JOB_STORE = jobflow.SETTINGS.JOB_STORE
else:
    logger.info(f"DEFECT_JOBSTORE set to {defect_jobstore_path}")
    JOB_STORE = loadfn(defect_jobstore_path)


def jdoc_to_entry(
    doc: dict, inc_structure: bool = False
) -> ComputedStructureEntry | ComputedEntry:
    """Get the entry from a job doc.

    Args:
        doc: job doc
        inc_structure: include the structure

    Returns:
        ComputedStructureEntry | ComputedEntry
    """
    entry_dict = doc["output"]["entry"]
    entry_id = doc["uuid"]
    if inc_structure:
        structure = mdecode(doc["output"]["structure"])
        return ComputedStructureEntry.from_dict(
            {**entry_dict, "structure": structure, "entry_id": entry_id}
        )
    return ComputedEntry.from_dict({**entry_dict, "entry_id": entry_id})


def update_metadata(flow: jobflow.Flow | jobflow.Job, metadata_updates: dict) -> None:
    """Update the metadata for a flow.

    Args:
        flow: The flow to update.
        metadata_updates: The updates to make.
    """
    for job in flow.jobs:
        if isinstance(job, jobflow.Flow):
            update_metadata(job, metadata_updates)
        else:
            job.metadata.update(metadata_updates)
