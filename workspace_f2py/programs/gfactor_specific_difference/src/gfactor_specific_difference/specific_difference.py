from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .cases import IonCase
from .dirac_bspline_interface import run_hlike_dirac, run_lilike_dirac
from .gfactor_evaluation import GFactorComponents
from .uncertainty_propagation import centered_finite_difference, secant_slope
from .weak_contribution import WeakContributionResult


@dataclass
class FNSResponseResult:
    case: IonCase
    system: str
    rms_values_fm: tuple[float, ...]
    g_values: tuple[float, ...]
    derivative_value: float
    reference_rms_fm: float
    source_paths: tuple[str, ...]
    notes: tuple[str, ...]


@dataclass
class XiResult:
    li_case: IonCase
    h_case: IonCase
    value: float
    li_derivative: float
    h_derivative: float
    source_paths: tuple[str, ...]
    notes: tuple[str, ...]


@dataclass
class SpecificDifferenceResult:
    li_case: IonCase
    h_case: IonCase
    xi: XiResult
    li_value: float
    h_value: float
    delta_xi_g: float
    weak_residual_estimate: float | None
    weak_residual_mode: str
    source_paths: tuple[str, ...]
    notes: tuple[str, ...]


def _normalize_system_label(system: str) -> str:
    return system.strip().lower().replace("-", "").replace(" ", "")


def estimate_fns_response(case: IonCase, system: str, rms_steps_fm: Sequence[float]) -> FNSResponseResult:
    rms_values = tuple(float(value) for value in rms_steps_fm)
    if len(rms_values) < 2:
        raise ValueError("At least two RMS-radius points are required to estimate the FNS response.")

    normalized_system = _normalize_system_label(system)
    if normalized_system == "hlike" and not case.is_hlike:
        raise ValueError("The provided case is not H-like.")
    if normalized_system == "lilike" and not case.is_lilike:
        raise ValueError("The provided case is not Li-like.")

    g_values: list[float] = []
    source_paths: list[str] = []
    notes: list[str] = []
    for rms_fm in rms_values:
        shifted_case = case.with_rms_fm(rms_fm)
        if normalized_system == "hlike":
            components = run_hlike_dirac(shifted_case)
        elif normalized_system == "lilike":
            components = run_lilike_dirac(shifted_case)
        else:
            raise ValueError(f"Unsupported system label {system!r}. Use 'hlike' or 'lilike'.")
        g_values.append(components.selected_value)
        source_paths.extend(components.source_paths)
        notes.extend(components.notes)

    if len(rms_values) == 3 and abs(rms_values[1] - case.rms_fm) <= 1.0e-12:
        derivative = centered_finite_difference(rms_values[0], rms_values[2], g_values[0], g_values[2])
        deriv_note = "Centered finite difference using the outer scan points."
    else:
        derivative = secant_slope(rms_values, tuple(g_values))
        deriv_note = "Secant slope across the provided RMS-radius range."

    return FNSResponseResult(
        case=case,
        system=system,
        rms_values_fm=rms_values,
        g_values=tuple(g_values),
        derivative_value=derivative,
        reference_rms_fm=case.rms_fm,
        source_paths=tuple(dict.fromkeys(source_paths)),
        notes=tuple(dict.fromkeys(notes + [deriv_note])),
    )


def determine_xi(
    li_case: IonCase,
    h_case: IonCase,
    fns_li: FNSResponseResult,
    fns_h: FNSResponseResult,
) -> XiResult:
    if li_case.znucl != h_case.znucl or li_case.nneutrons != h_case.nneutrons:
        raise ValueError("The H-like and Li-like cases must refer to the same nucleus when determining xi.")
    if fns_h.derivative_value == 0.0:
        raise ZeroDivisionError("The H-like FNS derivative is zero, so xi cannot be determined.")

    xi_value = fns_li.derivative_value / fns_h.derivative_value
    return XiResult(
        li_case=li_case,
        h_case=h_case,
        value=xi_value,
        li_derivative=fns_li.derivative_value,
        h_derivative=fns_h.derivative_value,
        source_paths=tuple(dict.fromkeys(fns_li.source_paths + fns_h.source_paths)),
        notes=(
            "Initial xi is determined from the Dirac+FNS response only.",
            "The cancellation condition is imposed on the leading numerical FNS derivative.",
        ),
    )


def build_specific_difference(
    li_g: GFactorComponents,
    h_g: GFactorComponents,
    xi: XiResult,
    weak_li: WeakContributionResult | None,
    weak_h: WeakContributionResult | None,
) -> SpecificDifferenceResult:
    delta_xi_g = li_g.selected_value - xi.value * h_g.selected_value
    weak_residual_estimate: float | None = None
    weak_mode = "not_available"
    notes = [
        "Specific difference currently combines the available copied H-like and Li-like baseline inputs.",
    ]

    if weak_li is not None and weak_h is not None and weak_li.selected_value is not None and weak_h.selected_value is not None:
        weak_residual_estimate = weak_li.selected_value - xi.value * weak_h.selected_value
        weak_mode = "complete"
    elif weak_h is not None and weak_h.selected_value is not None:
        weak_residual_estimate = -xi.value * weak_h.selected_value
        weak_mode = "hlike_only_partial"
        notes.append("Only the H-like weak contribution is currently available; the residual is partial.")
    else:
        notes.append("Weak residual estimate is not yet available from a validated weak-g backend.")

    source_paths = tuple(
        dict.fromkeys(
            li_g.source_paths
            + h_g.source_paths
            + xi.source_paths
            + (() if weak_li is None else weak_li.source_paths)
            + (() if weak_h is None else weak_h.source_paths)
        )
    )
    return SpecificDifferenceResult(
        li_case=li_g.case,
        h_case=h_g.case,
        xi=xi,
        li_value=li_g.selected_value,
        h_value=h_g.selected_value,
        delta_xi_g=delta_xi_g,
        weak_residual_estimate=weak_residual_estimate,
        weak_residual_mode=weak_mode,
        source_paths=source_paths,
        notes=tuple(notes),
    )
