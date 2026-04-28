from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from ._paths import ensure_python_api_on_path
from .cases import IonCase


ensure_python_api_on_path()


@dataclass
class WeakContributionResult:
    case: IonCase
    status: str
    selected_value: float | None
    weak_potential_summary: dict[str, float]
    source_paths: tuple[str, ...]
    notes: tuple[str, ...]


def _workspace_api():
    import workspace_f2py_api

    return workspace_f2py_api


def _build_rad_f2py_if_needed() -> None:
    api = _workspace_api()
    try:
        api.load_rad_kernel_suite_module()
    except ModuleNotFoundError:
        api.build_rad_f2py()


@lru_cache(maxsize=None)
def _cached_weak_summary(znucl: int, nneutrons: int, r_fm: float) -> tuple[float, float, float, float, float, int]:
    api = _workspace_api()
    _build_rad_f2py_if_needed()
    return tuple(api.compute_rad_weak_potential_summary(znucl, nneutrons, r_fm))


def summarize_weak_potential(case: IonCase) -> WeakContributionResult:
    r_rel, rho0, qweak, pc_center, pv_center, ierr = _cached_weak_summary(case.znucl, case.nneutrons, case.rms_fm)
    status = "ready" if ierr == 0 else "error"
    notes = (
        "This summary reuses the copied weak-potential normalization from workspace_f2py/programs/rad_f2py/src/rad_kernel_f2py.f90.",
    )
    if ierr != 0:
        notes = notes + ("The weak-potential summary returned a nonzero ierr flag.",)
    return WeakContributionResult(
        case=case,
        status=status,
        selected_value=None,
        weak_potential_summary={
            "r_rel": float(r_rel),
            "rho0": float(rho0),
            "qweak": float(qweak),
            "pc_center": float(pc_center),
            "pv_center": float(pv_center),
            "ierr": float(ierr),
        },
        source_paths=(
            "workspace_f2py/programs/rad_f2py/src/rad_kernel_f2py.f90",
            "workspace_f2py/codes/bspqed/src/bspqed_potentials.f90",
        ),
        notes=notes,
    )


def estimate_hlike_weak_contribution(case: IonCase) -> WeakContributionResult:
    summary = summarize_weak_potential(case)
    return WeakContributionResult(
        case=case,
        status="provisional_not_implemented",
        selected_value=None,
        weak_potential_summary=summary.weak_potential_summary,
        source_paths=summary.source_paths + ("Reading/Thesisvf2.pdf",),
        notes=summary.notes
        + (
            "No validated legacy weak-g driver exists in the repository yet.",
            "The weak-g contribution must be assembled from reused DKB orbitals and weak-potential machinery in a later milestone.",
        ),
    )
