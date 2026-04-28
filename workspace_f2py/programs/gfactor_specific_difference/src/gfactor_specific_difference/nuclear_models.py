from __future__ import annotations

import math

from .cases import IonCase


# These constants intentionally mirror the values in
# workspace_f2py/codes/bspqed/src/bspqed_physconsts.f90.
HBARC_MEV_FM = 197.3269788
ELECTRON_REST_ENERGY_MEV = 0.5109989461
ELECTRON_REDUCED_COMPTON_FM = HBARC_MEV_FM / ELECTRON_REST_ENERGY_MEV


def electron_reduced_compton_fm() -> float:
    return ELECTRON_REDUCED_COMPTON_FM


def rms_fm_to_spherical_radius_fm(rms_fm: float) -> float:
    return math.sqrt(5.0 / 3.0) * rms_fm


def rms_fm_to_spherical_radius_relativistic(rms_fm: float) -> float:
    return rms_fm_to_spherical_radius_fm(rms_fm) / electron_reduced_compton_fm()


def convert_r_fm_to_relativistic(r_fm: float) -> float:
    return r_fm / electron_reduced_compton_fm()


def summarize_nuclear_model(case: IonCase) -> dict[str, object]:
    return {
        "charge_model": case.nuclear_model.charge_model,
        "weak_model": case.nuclear_model.weak_model,
        "charge_radius_convention": case.nuclear_model.charge_radius_convention,
        "rms_fm": case.rms_fm,
        "sphere_radius_fm": rms_fm_to_spherical_radius_fm(case.rms_fm),
        "sphere_radius_relativistic": rms_fm_to_spherical_radius_relativistic(case.rms_fm),
        "weak_radius_relativistic": convert_r_fm_to_relativistic(case.rms_fm),
        "notes": case.nuclear_model.notes,
        "source_path": "workspace_f2py/codes/bspqed/src/bspqed_physconsts.f90",
    }
