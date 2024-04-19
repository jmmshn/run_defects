"""Utility functions for the defects module."""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING, Any

import jobflow
from fireworks import LaunchPad
from monty.json import MontyDecoder
from monty.serialization import loadfn

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
        jobs = []
        for maker in self.makers:
            job = maker.make(prv_struct)
            jobs.append(job)
            prv_struct = job.output.structure
        return jobflow.Flow(jobs)


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
