from __future__ import annotations

import json
from pathlib import Path

from .gfactor_evaluation import GFactorComponents, summarize_components
from .specific_difference import FNSResponseResult, SpecificDifferenceResult, XiResult
from .weak_contribution import WeakContributionResult


def summarize_fns_response(result: FNSResponseResult) -> dict[str, object]:
    return {
        "label": result.case.label,
        "system": result.system,
        "reference_rms_fm": result.reference_rms_fm,
        "rms_values_fm": list(result.rms_values_fm),
        "g_values": list(result.g_values),
        "derivative_value": result.derivative_value,
        "source_paths": list(result.source_paths),
        "notes": list(result.notes),
    }


def summarize_xi(result: XiResult) -> dict[str, object]:
    return {
        "li_label": result.li_case.label,
        "h_label": result.h_case.label,
        "xi": result.value,
        "li_derivative": result.li_derivative,
        "h_derivative": result.h_derivative,
        "source_paths": list(result.source_paths),
        "notes": list(result.notes),
    }


def summarize_weak_result(result: WeakContributionResult) -> dict[str, object]:
    return {
        "label": result.case.label,
        "status": result.status,
        "selected_value": result.selected_value,
        "weak_potential_summary": result.weak_potential_summary,
        "source_paths": list(result.source_paths),
        "notes": list(result.notes),
    }


def summarize_specific_difference(result: SpecificDifferenceResult) -> dict[str, object]:
    return {
        "li_label": result.li_case.label,
        "h_label": result.h_case.label,
        "xi": result.xi.value,
        "li_value": result.li_value,
        "h_value": result.h_value,
        "delta_xi_g": result.delta_xi_g,
        "weak_residual_estimate": result.weak_residual_estimate,
        "weak_residual_mode": result.weak_residual_mode,
        "source_paths": list(result.source_paths),
        "notes": list(result.notes),
    }


def write_json_report(path: str | Path, payload: dict[str, object]) -> None:
    path_obj = Path(path)
    path_obj.parent.mkdir(parents=True, exist_ok=True)
    path_obj.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def summarize_run(
    h_components: GFactorComponents,
    li_components: GFactorComponents,
    fns_h: FNSResponseResult,
    fns_li: FNSResponseResult,
    xi: XiResult,
    specific_difference: SpecificDifferenceResult,
    weak_h: WeakContributionResult | None = None,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "hydrogen_like": summarize_components(h_components),
        "lithium_like": summarize_components(li_components),
        "fns_hydrogen_like": summarize_fns_response(fns_h),
        "fns_lithium_like": summarize_fns_response(fns_li),
        "xi": summarize_xi(xi),
        "specific_difference": summarize_specific_difference(specific_difference),
    }
    if weak_h is not None:
        payload["weak_hydrogen_like"] = summarize_weak_result(weak_h)
    return payload
