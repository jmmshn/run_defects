"""Submit the workflows to firework."""

from __future__ import annotations

from typing import TYPE_CHECKING

import jobflow
from atomate2.vasp.flows.defect import FormationEnergyMaker
from atomate2.vasp.flows.mp import MPGGADoubleRelaxMaker, MPGGAStaticMaker
from atomate2.vasp.jobs.mp import MPGGARelaxMaker
from atomate2.vasp.powerups import (
    update_user_incar_settings,
    update_user_kpoints_settings,
)
from atomate2.vasp.sets.defect import SPECIAL_KPOINT
from atomate2.vasp.sets.mp import MPGGARelaxSetGenerator
from jobflow.managers.fireworks import flow_to_workflow

if TYPE_CHECKING:
    from fireworks import LaunchPad
    from pymatgen.core import Structure

INCAR_UPDATES = {
    "EDIFF": 1e-5,
    "EDIFFG": -0.05,
    "NELMIN": 5,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "NELM": 300,
    "ENCUT": 520,
    "ALGO": "Normal",
    "POTIM": 0.25,
    "LREAL": "Auto",
}

BULK_RELAX_UC = update_user_incar_settings(
    MPGGADoubleRelaxMaker(), incar_updates=INCAR_UPDATES
)

BULK_STATIC_UC = MPGGAStaticMaker(
    task_document_kwargs={"store_volumetric_data": ["locpot", "chgcar"]},
)

BULK_STATIC_UC = update_user_incar_settings(
    BULK_STATIC_UC, incar_updates=INCAR_UPDATES | {"LVHAR": True}
)

BULK_STATIC_SC = update_user_incar_settings(
    BULK_STATIC_UC, incar_updates={"LREAL": "Auto"}
)
BULK_STATIC_SC = update_user_kpoints_settings(
    BULK_STATIC_SC, kpoints_updates=SPECIAL_KPOINT
)

DEFECT_RELAX_SC = update_user_incar_settings(
    MPGGARelaxMaker(
        input_set_generator=MPGGARelaxSetGenerator(use_structure_charge=True),
        task_document_kwargs={"store_volumetric_data": ["locpot"]},
    ),
    incar_updates=INCAR_UPDATES | {"ISIF": 2, "LVHAR": True},
)
DEFECT_RELAX_SC = update_user_kpoints_settings(
    DEFECT_RELAX_SC, kpoints_updates=SPECIAL_KPOINT
)

F_MAKER_UC = FormationEnergyMaker(
    defect_relax_maker=DEFECT_RELAX_SC,
    uc_bulk=True,
    perturb=0.2,
    collect_defect_entry_data=False,
    relax_radius="auto",
)

F_MAKER_SC = FormationEnergyMaker(
    defect_relax_maker=DEFECT_RELAX_SC,
    uc_bulk=False,
    perturb=0.2,
    collect_defect_entry_data=False,
    relax_radius="auto",
)


def add_tag(flow: jobflow.Flow, tag: str) -> None:
    """Update the metadata of the flow to include a tag."""
    for job in flow.jobs:
        if isinstance(job, jobflow.Flow):
            add_tag(job, tag)
        else:
            tags = set(job.metadata.get("tags", set())) | {
                tag,
            }
            job.metadata.update({"tags": list(tags)})


def submit_flow(flow: jobflow.Flow, tag: str, lpad: LaunchPad) -> None:
    """Submit the flow to firework.

    Args:
        flow (Flow): flow to submit to firework
        tag (str): tag to add to the flow
        lpad (LaunchPad): firework launchpad
    """
    add_tag(flow, tag)
    lpad.add_wf(flow_to_workflow(flow))


def get_bulk_flow(
    structure: Structure,
    relax_maker: jobflow.Maker | None = None,
    static_maker: jobflow.Maker | None = None,
    remove_symmetry: bool = True,
) -> jobflow.Flow:
    """Get the bulk workflow for a structure.

    (double relax, static)

    Args:
        structure (Structure): structure to compute
        relax_maker (jobflow.Maker, optional): maker for the relax job.
            Defaults to None.
        static_maker (jobflow.Maker, optional): maker for the static job.
            Defaults to None.
        remove_symmetry (bool, optional): remove symmetry from the structure.
            Defaults to True.

    Returns:
        jobflow.Flow: bulk flow
    """
    if relax_maker is None:
        relax_maker = BULK_RELAX_UC
    if static_maker is None:
        static_maker = BULK_STATIC_UC

    struct_ = structure.copy()
    struct_.remove_oxidation_states()
    struct_.remove_site_property("magmom")
    relax_job = relax_maker.make(struct_)
    jobs = [relax_job]
    formula = struct_.composition.reduced_formula
    static_job = static_maker.make(relax_job.output.structure)
    jobs.append(static_job)
    flow = jobflow.Flow(jobs, name=f"{formula} bulk")
    if remove_symmetry:
        magmom_vals = {el: 0.6 for el in map(str, struct_.elements)}
        flow = update_user_incar_settings(flow, incar_updates={"MAGMOM": magmom_vals})
    return flow
