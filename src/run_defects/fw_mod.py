"""Modify atomate2 fireworks in database."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from fireworks import LaunchPad
from monty.json import MontyDecoder, jsanitize
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Incar, Poscar

from run_defects.utils import rget

if TYPE_CHECKING:
    from collections.abc import Generator


logger = logging.getLogger(__name__)

LPAD_AUTO_LOAD = LaunchPad.auto_load()


def _get_fw_id_by_state(
    query: dict | None = None, state: str = "FIZZLED", lpad: LaunchPad = None
) -> list:
    """Find the fizzed relaxation jobs.

    The atomate2 created fireworks for long running relaxation jobs should be simple.
    Basically `spec._tasks.0` should be the relaxation job.

    Args:
        query: query to use.  The basic query is already included.
        state: state to search for, default FIZZLED
        lpad: LaunchPad object


    Returns:
        list: list of fizzled fireworks ids
    """
    if lpad is None:
        lpad = LPAD_AUTO_LOAD
    fw_query = {
        "spec._tasks.0.job.function.@callable": "BaseVaspMaker.make",
        "state": state,
        "$or": [
            {"spec._tasks.0.job.function_args.0.@class": "Structure"},
            {"spec._tasks.0.job.function_args.0.@class": "OutputReference"},
            {"spec._tasks.0.job.function_args.0.@class": "Poscar"},
        ],
    }
    fw_query.update(query or {})
    return lpad.get_fw_ids(query=fw_query)


def _get_last_launch_dir(fw_dict: dict) -> str:
    """Get the last launch directory for a fizzled job.

    Args:
        fw_dict: fireworks id

    Returns:
        str: last launch directory
    """
    launches = fw_dict.get("launches")
    if launches is None:
        launches = fw_dict.get("archived_launches")
    if launches is None or len(launches) == 0:
        raise ValueError(f"Could not find launches for {fw_dict['fw_id']}")
    return launches[-1]["launch_dir"]


def _read_contcar(dir_name: str) -> Poscar:
    """Read the CONTCAR currently in the directory.

    Args:
            dir_name: directory to read CONTCAR from

    Returns:
            Structure: structure from CONTCAR
    """
    # check for CONTCAR or CONTCAR.gz
    sub_dir = Path(dir_name)
    try:
        contcar_file = next(sub_dir.glob("CONTCAR*"))
        contcar = Poscar.from_file(contcar_file)
    except Exception:
        contcar = Poscar.from_file(sub_dir / "POSCAR")
    return contcar


def _read_incar(dir_name: str) -> tuple[Incar, Incar]:
    """Read the INCAR currently in the directory.

    Args:
            dir_name: directory to read INCAR from

    Returns:
            incar_orig: original INCAR
            incar: current INCAR
    """
    try:
        incar_orig = Incar.from_file(Path(dir_name) / "INCAR.orig")
        incar = Incar.from_file(Path(dir_name) / "INCAR")
    except Exception as e:
        raise RuntimeError(f"Could not read INCAR from {dir_name}") from e
        return None, None
    return incar_orig, incar


def _recursive_update(
    orig: dict | list, update_dict: dict, path: str = "_tasks"
) -> Generator:
    """Update a dictionary recursively.

    Check if orig[key] is in update_dict, if it is then update it.
    Else, call _recursive_update(orig[key], update_dict)
    """
    if isinstance(orig, dict):
        for key, v in orig.items():
            if key in update_dict:
                if not isinstance(v, dict):
                    raise RuntimeError("Can only update dict.")
                for k_, v_ in update_dict[key].items():
                    yield f"{path}.{key}.{k_}", v_
            new_path = f"{path}.{key}"
            if "archived_launches" not in path:
                yield from _recursive_update(orig[key], update_dict, path=new_path)

    elif isinstance(orig, list):
        for i, v in enumerate(orig):
            new_path = f"{path}.{i}"
            if "archived_launches" not in path:
                yield from _recursive_update(v, update_dict, path=new_path)


def _update_fizzled_firework(
    fw_id: int,
    update_dict: dict,
    dry_run: bool = True,
    check_fizzled: bool = True,
    lpad: LaunchPad = None,
) -> None:
    """Update a fizzled firework.

    Args:
        fw_id: fireworks id.
        update_dict: dictionary of updates.
        dry_run: dry run, default True, only print update dict.
        check_fizzled: check if the firework is fizzled, default True.
        lpad: LaunchPad object, default will use auto_load.
    """
    # make sure the firework is fizzled
    if lpad is None:
        lpad = LPAD_AUTO_LOAD
    fw_dict = lpad.get_fw_dict_by_id(fw_id)
    if check_fizzled and fw_dict["state"] != "FIZZLED":
        raise RuntimeError("Firework must be fizzled.")

    updates_dict = dict(_recursive_update(fw_dict["spec"]["_tasks"], update_dict))
    if dry_run:
        logger.debug(f"Would update {fw_id} with {updates_dict}")
    else:
        lpad.update_spec([fw_id], updates_dict)


def update_incar_settings(
    fw_id: int,
    user_incar_setting: dict,
    check_fizzled: bool = True,
    dry_run: bool = True,
) -> None:
    """Update a fizzled firework.

    Args:
        fw_id: fireworks id.
        user_incar_setting: dictionary of updates. Ex: {"EDIFFG": -0.005}
        check_fizzled: check if the firework is fizzled, default True.
        dry_run: dry run, default True, only print update dict.
    """
    _update_fizzled_firework(
        fw_id=fw_id,
        update_dict={"user_incar_settings": user_incar_setting},
        dry_run=dry_run,
        check_fizzled=check_fizzled,
    )


def update_struct(fw_id: int, structure: Structure, lpad: LaunchPad = None) -> None:
    """Update the spec structure for a fizzled job.

    Args:
        fw_id: fireworks id.
        structure: structure to update to.
        lpad: LaunchPad object, default will use auto_load.
    """
    if lpad is None:
        lpad = LPAD_AUTO_LOAD
    update_dict = {
        "_tasks.0.job.function_args.0": jsanitize(structure.as_dict()),
    }
    lpad.update_spec([fw_id], update_dict)


def _get_incar_diff(incar_old: Incar, incar_new: Incar) -> dict:
    """Get the INCAR differences between two INCARs."""
    incar_diff = {}
    for k, v in incar_new.items():
        if incar_old.get(k, None) != v:
            incar_diff[k] = v
    return incar_diff


def get_charged_structure(fw_id: int, lpad: LaunchPad = None) -> Structure:
    """Read the CONTCAR if available and decorate with the charge state.

    Args:
        fw_id: fireworks id
        lpad: LaunchPad object, default will use auto_load.
    """
    if lpad is None:
        lpad = LPAD_AUTO_LOAD
    fw_dict = lpad.get_fw_dict_by_id(fw_id)
    dir_name = _get_last_launch_dir(fw_dict)
    if dir_name is None:
        raise RuntimeError(f"Could not find a launch directory for {fw_id}")
    logger.debug(f"Get the charged structure for {fw_id} from {dir_name}")
    contcar = _read_contcar(dir_name)
    structure: Structure = contcar.structure
    defect_charge = None
    try:
        defect_charge_path = (
            "spec._tasks.0.job.function.@bound"
            ".write_additional_data.info:json.charge_state"
        )
        defect_charge = rget(fw_dict, defect_charge_path)
    except KeyError:
        pass

    if defect_charge is None:
        old_struct = fw_dict["spec"]["_tasks"][0]["job"]["function_args"][0]
        old_struct = MontyDecoder().process_decoded(
            fw_dict["spec"]["_tasks"][0]["job"]["function_args"][0]
        )
        defect_charge = old_struct._charge
        if not isinstance(old_struct, Structure):
            raise TypeError("Check that the arg is parsed correctly")
    structure._charge = defect_charge
    return structure
