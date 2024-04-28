"""Restart fizzled FireWorks with updated INCAR settings."""

from __future__ import annotations

import json
import logging

import click
from fireworks import LaunchPad

from run_defects.fw_mod import (
    _get_fw_id_by_state,
    _get_incar_diff,
    _get_last_launch_dir,
    _read_incar,
    get_charged_structure,
    update_incar_settings,
    update_struct,
)

LPAD = LaunchPad.auto_load()
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


# # %%
# ready_ids = _get_fw_id_by_state(state="READY")
# for fw_id in ready_ids:
#     last_launch_dir = _get_last_launch_dir(LPAD.get_fw_dict_by_id(fw_id))
#     update_incar_settings(fw_id, {"LREAL": "Auto"}, check_fizzled=False, dry_run=True)

# #%%
# fw_id = _get_fizzled()[0]
# fw_dict = LPAD.get_fw_dict_by_id(fw_id)
# last_launch_dir = _get_last_launch_dir(fw_dict)
# contcar = _read_contcar(last_launch_dir)
# incar_orig, incar = _read_incar(last_launch_dir)
# incar_update_dict = _get_incar_diff(incar_orig, incar)

# # get the incar values that are different
# incar_diff = _get_incar_diff(incar_orig, incar)
# if swap_ibrion:
#     _swap_val = {1: 2, 2: 1}
#     incar_diff["IBRION"] = _swap_val[incar["IBRION"]]

# update_fizzled_incar_settings(fw_id, incar_diff, dry_run=False)
# cstruct = get_charged_structure(fw_id)
# update_struct(fw_id, cstruct)


# %%
@click.command()
@click.option("-q", "--query", default={})
@click.option("-d", "--dry-run", is_flag=True, show_default=True, default=False)
@click.option("-i", "--fw-id")
@click.option(
    "-b",
    "--ibrion-swap",
    is_flag=True,
    show_default=True,
    default=False,
    help="Swap IBRION value.",
)
def main(
    query: str, dry_run: bool, fw_id: int | None = None, ibrion_swap: bool = False
) -> None:
    """Run the CLI."""
    q = json.loads(query)
    fizzled = [int(fw_id)] if fw_id is not None else _get_fw_id_by_state(q)
    for fw_id in fizzled:
        structure = get_charged_structure(fw_id)
        last_launch_dir = _get_last_launch_dir(LPAD.get_fw_dict_by_id(fw_id))
        incar_orig, incar = _read_incar(last_launch_dir)
        incar_update_dict = _get_incar_diff(incar_orig, incar)
        logger.info(f"Launch Dir: {last_launch_dir}")
        logger.info(f"INCAR Update: {incar_update_dict}")
        logger.info(f"Parsed Structure: {structure.formula}")

        if not dry_run:
            update_struct(fw_id, structure)
            if ibrion_swap:
                _swap_val = {1: 2, 2: 1, 3: 1}
                incar_update_dict["IBRION"] = _swap_val[incar["IBRION"]]

            update_incar_settings(fw_id, incar_update_dict, dry_run=False)
            LPAD.rerun_fw(fw_id)


# %%
if __name__ == "__main__":
    main()
