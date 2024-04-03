"""Submit the workflows to firework."""

from __future__ import annotations

import itertools
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
from pymatgen.analysis.defects.generators import (
    ChargeInterstitialGenerator,
    SubstitutionGenerator,
    VacancyGenerator,
)

from run_defects.utils import LPAD, rget

VGEN = VacancyGenerator()
IGEN = ChargeInterstitialGenerator(max_insertions=3)
SGEN = SubstitutionGenerator()

if TYPE_CHECKING:
    from collections.abc import Generator

    from fireworks import LaunchPad
    from pymatgen.analysis.defects.core import Defect
    from pymatgen.core import Structure
    from pymatgen.io.vasp.outputs import VolumetricData

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
}

BULK_RELAX_UC = update_user_incar_settings(
    MPGGADoubleRelaxMaker(), incar_updates=INCAR_UPDATES
)

BULK_STATIC_UC = MPGGAStaticMaker(
    task_document_kwargs={"store_volumetric_data": ["locpot", "chgcar"]},
)

BULK_STATIC_UC = update_user_incar_settings(
    BULK_STATIC_UC, incar_updates=INCAR_UPDATES | {"LVHAR": True, "LREAL": False}
)

BULK_STATIC_SC = update_user_kpoints_settings(
    BULK_STATIC_UC, kpoints_updates=SPECIAL_KPOINT
)

DEFECT_RELAX_SC = update_user_incar_settings(
    MPGGARelaxMaker(
        input_set_generator=MPGGARelaxSetGenerator(use_structure_charge=True),
        task_document_kwargs={"store_volumetric_data": ["locpot"]},
    ),
    incar_updates=INCAR_UPDATES | {"ISIF": 2, "LVHAR": True, "LREAL": "Auto"},
)

DEFECT_RELAX_SC = update_user_kpoints_settings(
    DEFECT_RELAX_SC, kpoints_updates=SPECIAL_KPOINT
)

DEFECT_STATIC_SC = MPGGAStaticMaker(
    task_document_kwargs={"store_volumetric_data": ["locpot"]},
)

DEFECT_STATIC_SC = update_user_incar_settings(
    DEFECT_STATIC_SC, incar_updates=INCAR_UPDATES | {"LVHAR": True, "LREAL": False}
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


def _get_elements(struct: Structure) -> tuple[tuple[str], ...]:
    """Get the elements inputs for vacancy and interstitial defects."""
    s_atoms = sorted(atom.symbol for atom in struct.elements)
    return tuple((aa_,) for aa_ in s_atoms)


def _get_subs(struct: Structure) -> tuple[dict[str, str], ...]:
    """Get the substitutions."""
    s_atoms = sorted(atom.symbol for atom in struct.elements)
    return tuple(
        itertools.chain.from_iterable(
            ({a: b}, {b: a}) for a, b in itertools.combinations(s_atoms, 2)
        )
    )


def get_defects(
    chgcar: VolumetricData, max_iter: int = 3
) -> Generator[tuple[Defect, int], None, None]:
    """Generate the defects for a chgcar."""
    tup_el = _get_elements(chgcar.structure)
    tup_sub = _get_subs(chgcar.structure)
    for sub_d in tup_sub[:max_iter]:
        for ii, defect in enumerate(SGEN.generate(chgcar.structure, sub_d)):
            yield defect, ii

    for el in tup_el[:max_iter]:
        for ii, defect in enumerate(VGEN.generate(chgcar.structure, el)):
            yield defect, ii

    for el in tup_el[:max_iter]:
        for ii, defect in enumerate(
            IGEN.generate(
                chgcar,
                el,
            )
        ):
            yield defect, ii


def get_submitted_defect_run_ids() -> set:
    """Get the defect run ids that have been submitted."""
    submitted_dicts = list(
        LPAD.fireworks.find(
            {
                "spec._tasks.0.job.metadata.defect_run_id": {"$exists": 1},
            },
            projection={"spec._tasks.job.metadata.defect_run_id": True},
        )
    )
    return {
        rget(dd_, "spec._tasks.0.job.metadata.defect_run_id") for dd_ in submitted_dicts
    }
