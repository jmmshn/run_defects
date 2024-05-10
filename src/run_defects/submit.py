"""Submit the workflows to firework."""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import jobflow
from atomate2.common.jobs.defect import check_charge_state
from atomate2.vasp.flows.defect import FormationEnergyMaker
from atomate2.vasp.flows.mp import MPGGADoubleRelaxMaker, MPGGAStaticMaker
from atomate2.vasp.jobs.mp import MPGGARelaxMaker
from atomate2.vasp.powerups import (
    update_user_incar_settings,
    update_user_kpoints_settings,
)
from atomate2.vasp.sets.defect import SPECIAL_KPOINT
from atomate2.vasp.sets.mp import MPGGARelaxSetGenerator, MPGGAStaticSetGenerator
from jobflow.managers.fireworks import flow_to_workflow
from pymatgen.analysis.defects.generators import (
    ChargeInterstitialGenerator,
    SubstitutionGenerator,
    VacancyGenerator,
)
from pymatgen.core import Structure
from tqdm import tqdm

from run_defects.utils import JOB_STORE, LPAD, ComboMaker, rget

VGEN = VacancyGenerator()
IGEN = ChargeInterstitialGenerator(max_insertions=3)
SGEN = SubstitutionGenerator()

if TYPE_CHECKING:
    from collections.abc import Generator

    from fireworks import LaunchPad
    from jobflow import Maker
    from maggma.stores import Store
    from pymatgen.analysis.defects.core import Defect
    from pymatgen.io.vasp.outputs import VolumetricData

INCAR_UPDATES = {
    "EDIFF": 1e-5,
    "EDIFFG": -0.05,
    "NELMIN": 5,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "NELM": 100,
    "ENCUT": 520,
    "ALGO": "Normal",
    "POTIM": 0.25,
    "LORBIT": 11,
    "NCORE": 4,
}

HSE_INCAR_UPDATES = {
    "ALGO": "Damped",
    "TIME": 0.5,
    "HFSCREEN": 0.2,
    "GGA": "PE",
    "LHFCALC": True,
    "PRECFOCK": "Normal",
}

BULK_RELAX_UC = update_user_incar_settings(
    MPGGADoubleRelaxMaker(), incar_updates=INCAR_UPDATES | {"KPAR": 2}
)

BULK_RELAX_SC = update_user_kpoints_settings(
    BULK_RELAX_UC, kpoints_updates=SPECIAL_KPOINT
)

BULK_STATIC_UC = MPGGAStaticMaker(
    task_document_kwargs={"store_volumetric_data": ["locpot", "chgcar"]},
)

BULK_STATIC_UC = update_user_incar_settings(
    BULK_STATIC_UC,
    incar_updates=INCAR_UPDATES | {"LVHAR": True, "LREAL": False, "KPAR": 4},
)

BULK_STATIC_UC_HSE = update_user_incar_settings(
    BULK_STATIC_UC,
    incar_updates=HSE_INCAR_UPDATES
    | {"KPAR": 4, "NCORE": 1, "LREAL": False, "LMAXMIX": 6},
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
DEFECT_RELAX_SC.name = "defect relax"

DEFECT_STATIC_SC = MPGGAStaticMaker(
    input_set_generator=MPGGAStaticSetGenerator(use_structure_charge=True),
    task_document_kwargs={"store_volumetric_data": ["locpot"]},
    copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
)
DEFECT_STATIC_SC.name = "defect static"

DEFECT_STATIC_SC = update_user_incar_settings(
    DEFECT_STATIC_SC,
    incar_updates=INCAR_UPDATES
    | {"LVHAR": True, "LREAL": False, "LMAXMIX": 6, "LWAVE": True},
)

DEFECT_STATIC_SC = update_user_kpoints_settings(
    DEFECT_STATIC_SC, kpoints_updates=SPECIAL_KPOINT
)

DEFECT_STATIC_SC_HSE = update_user_incar_settings(
    DEFECT_STATIC_SC, incar_updates=HSE_INCAR_UPDATES | {"NCORE": 1}
)
DEFECT_STATIC_SC_HSE.name = "defect static HSE"

F_MAKER_UC = FormationEnergyMaker(
    defect_relax_maker=DEFECT_RELAX_SC,
    uc_bulk=True,
    perturb=0.2,
    collect_defect_entry_data=False,
    relax_radius="auto",
)

F_MAKER_SC = FormationEnergyMaker(
    defect_relax_maker=DEFECT_RELAX_SC,
    uc_bulk=True,
    perturb=0.2,
    collect_defect_entry_data=False,
    relax_radius="auto",
)

COMBO_MAKER = ComboMaker([DEFECT_RELAX_SC, DEFECT_STATIC_SC, DEFECT_STATIC_SC_HSE])
F_MAKER_COMBO = FormationEnergyMaker(
    defect_relax_maker=COMBO_MAKER,
    uc_bulk=True,
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
    static_maker: jobflow.Maker,
    relax_maker: jobflow.Maker | None = None,
    additional_maker: jobflow.Maker | None = None,
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
        additional_maker (jobflow.Maker, optional): additional maker for the flow.

    Returns:
        jobflow.Flow: bulk flow
    """
    struct_ = structure.copy()
    formula = struct_.composition.reduced_formula
    if remove_symmetry:
        struct_.remove_oxidation_states()
        struct_.remove_site_property("magmom")

    jobs = []
    if relax_maker:
        relax_job = relax_maker.make(struct_)
        jobs.append(relax_job)
        relax_job.name = f"{formula} static"
        struct_out = relax_job.output.structure
    else:
        struct_out = struct_

    static_job = static_maker.make(struct_out)
    static_job.name = f"{formula} static"
    jobs.append(static_job)
    if additional_maker:
        add_job = additional_maker.make(struct_out)
        jobs.append(add_job)

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
    """Get the defect run ids that have been submitted to the queue."""
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


def _get_completed_defects_calcs(
    query: dict = None, jobstore: Store = None, **kwargs
) -> Generator[dict, None, None]:
    """Query for completed defect calculations."""
    query = query or {}
    js_query = {
        "output.vasp_objects.locpot.@class": {"$exists": True},
        "output.additional_json.info.defect_name": {"$exists": True},
        **query,
    }
    properties = ["output.additional_json", "output.structure"]
    jobstore = jobstore or JOB_STORE
    ndocs = jobstore.count(js_query)
    for job_doc in tqdm(
        jobstore.query(criteria=js_query, properties=properties, **kwargs), total=ndocs
    ):
        defect_info = job_doc["output"]["additional_json"]["info"]
        structure = Structure.from_dict(job_doc["output"]["structure"])
        charge_state = defect_info["charge_state"]
        structure_charge = structure._charge
        if int(charge_state) != int(structure_charge):
            raise ValueError(
                f"Defect OBJ charge state ({int(charge_state)}) "
                "does not match structure charge ({int(structure_charge)})"
            )
        yield {"info": defect_info, "structure": structure}


def _get_defect_flow_from_info(
    maker: Maker, add_info: dict, structure: Structure
) -> jobflow.Flow:
    """Create the defect flow for a structure."""
    structure = structure.copy()
    charge_state = add_info["charge_state"]
    structure._charge = charge_state
    job_ = maker.make(structure)
    if not maker.input_set_generator.use_structure_charge:
        raise ValueError("Maker must use structure charge")
    job_.update_maker_kwargs(
        {"_set": {"write_additional_data->info:json": add_info}}, dict_mod=True
    )
    job_.name = f"{add_info['bulk_formula']} {add_info['defect_name']}"
    check_job = check_charge_state(charge_state, job_.output.structure)
    return jobflow.Flow([job_, check_job], name=job_.name)


def get_defect_sc_job_from_completed(
    maker: Maker, query: dict = None, jobstore: Store = None, **kwargs
) -> Generator[dict, None, None]:
    """Query for completed GGA calculations and rerun them.

    Args:
        maker (Maker): maker for the defect calculations
        query (dict, optional): query for the defect calculations. Defaults to None.
        jobstore (Store, optional): jobstore to query. Defaults to None.
        **kwargs: additional kwargs for the query.
    """
    for dd_ in _get_completed_defects_calcs(query=query, jobstore=jobstore, **kwargs):
        yield _get_defect_flow_from_info(maker, dd_["info"], dd_["structure"])


def chunker(seq: list, size: int) -> Generator:
    """Chunk a sequence into chunks of size."""
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))
