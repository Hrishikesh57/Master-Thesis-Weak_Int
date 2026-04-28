#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from First_order import FirstOrderError, import_rad_kernel_suite
from Mass_derivative import MassDerivativeError, relative_difference
from Mass_derivative_direct_dkb import (
    DEFAULT_FERMI_A_FM,
    DirectDkbMassDerivativeError,
    OUTPUT_DIR,
    fermi_c_from_rms,
    load_sampled_reference,
    summarize_ion,
)


JSON_OUTPUT = OUTPUT_DIR / "fermi_route_convergence.json"
CSV_OUTPUT = OUTPUT_DIR / "fermi_route_convergence.csv"
DEFAULT_NPOINTS_VALUES = [600, 800, 1000, 1200, 1500]
NEUTRON_SKIN_VALUES_FM = [0.0, 0.10]
REFERENCE_NPOINTS_FOR_DECISION = 1000

ION_SET = [
    {
        "label": "Pb-208",
        "znucl": 82,
        "nneutrons": 126,
        "rms_fm": 5.5012,
        "proton_rms_fm": 5.5012,
    },
    {
        "label": "U-238",
        "znucl": 92,
        "nneutrons": 146,
        "rms_fm": 5.8571,
        "proton_rms_fm": 5.8571,
    },
]


class FermiRouteConvergenceError(RuntimeError):
    """User-facing final-route convergence error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="npoints convergence for Fermi Coulomb plus p/n Fermi weak-density direct-DKB route."
    )
    parser.add_argument(
        "--npoints",
        type=int,
        nargs="+",
        default=DEFAULT_NPOINTS_VALUES,
        help="B-spline grid point counts to sweep. Defaults to 600 800 1000 1200 1500.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for JSON and CSV convergence outputs.",
    )
    return parser.parse_args()


def validate_npoints_values(values: list[int]) -> list[int]:
    unique_values = sorted(set(values))
    for value in unique_values:
        if value < 20:
            raise FermiRouteConvergenceError(f"npoints={value} is too small for the DKB knot construction.")
        if value % 10 != 0:
            raise FermiRouteConvergenceError(
                f"npoints={value} is not divisible by 10; the DKB knot construction uses npoints/10."
            )
    return unique_values


def ion_cases() -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for ion in ION_SET:
        for skin_fm in NEUTRON_SKIN_VALUES_FM:
            neutron_rms_fm = ion["proton_rms_fm"] + skin_fm
            fermi_c_from_rms(ion["proton_rms_fm"], DEFAULT_FERMI_A_FM)
            fermi_c_from_rms(neutron_rms_fm, DEFAULT_FERMI_A_FM)
            cases.append(
                {
                    **ion,
                    "neutron_rms_fm": neutron_rms_fm,
                    "proton_a_fm": DEFAULT_FERMI_A_FM,
                    "neutron_a_fm": DEFAULT_FERMI_A_FM,
                    "neutron_skin_fm": skin_fm,
                }
            )
    return cases


def compact_record(record: dict[str, Any], elapsed_seconds: float) -> dict[str, Any]:
    mu_zero = next(step for step in record["steps"] if abs(step["h"] - 5.0e-4) < 1.0e-15)["mu_zero"]
    return {
        "isotope_label": record["isotope_label"],
        "Z": record["Z"],
        "N": record["N"],
        "neutron_skin_fm": record["neutron_skin_fm"],
        "proton_rms_fm": record["physical_proton_rms_fm"],
        "neutron_rms_fm": record["physical_neutron_rms_fm"],
        "proton_a_fm": record["physical_proton_a_fm"],
        "neutron_a_fm": record["physical_neutron_a_fm"],
        "proton_c_fm": record["physical_proton_c_fm"],
        "neutron_c_fm": record["physical_neutron_c_fm"],
        "npoints": record["npoints_requested"],
        "order": mu_zero["order"],
        "nknots": mu_zero["nknots"],
        "nbspl": mu_zero["nbspl"],
        "nspectrum": mu_zero["nspectrum"],
        "rcavity": mu_zero["rcavity"],
        "state_index": mu_zero["state_index"],
        "energy_total": mu_zero["energy_total"],
        "binding_energy": mu_zero["binding_energy"],
        "Delta_E_W_mu1": record["delta_epsilon_mu1"],
        "D_W_h_5e-4": record["derivative_dw_reference"],
        "delta_g_W_h_5e-4": record["delta_g_w_reference"],
        "D_W_h_1e-3": record["derivative_dw_by_h"]["1.000000e-03"],
        "delta_g_W_h_1e-3": record["delta_g_w_by_h"]["1.000000e-03"],
        "step_size_relative_difference_D_W": record["step_size_relative_difference_dw"],
        "step_size_relative_difference_delta_g_W": record["step_size_relative_difference_delta_g"],
        "proton_normalization_check_mu1": record["proton_normalization_check_mu1"],
        "neutron_normalization_check_mu1": record["neutron_normalization_check_mu1"],
        "fermi_coulomb_v_cp_mu1": record["fermi_coulomb_v_cp_mu1"],
        "fermi_coulomb_v_cp_plus_5ap_mu1": record["fermi_coulomb_v_cp_plus_5ap_mu1"],
        "elapsed_seconds": elapsed_seconds,
    }


def evaluate_case(
    module: Any,
    sampled_reference: dict[str, dict[str, Any]],
    ion: dict[str, Any],
    npoints_values: list[int],
) -> dict[str, Any]:
    settings: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for npoints in npoints_values:
        start = time.monotonic()
        print(
            f"Running {ion['label']} skin={ion['neutron_skin_fm']:.2f} fm, npoints={npoints}...",
            flush=True,
        )
        try:
            record = summarize_ion(module, ion, sampled_reference, npoints, "pn-fermi", "fermi")
        except Exception as exc:
            failures.append(
                {
                    "npoints": npoints,
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                }
            )
            print(f"  failed: {type(exc).__name__}: {exc}", flush=True)
            continue
        elapsed_seconds = time.monotonic() - start
        compact = compact_record(record, elapsed_seconds)
        settings.append(compact)
        print(
            f"  delta_g_W={compact['delta_g_W_h_5e-4']:.16e}, "
            f"step_rel={compact['step_size_relative_difference_delta_g_W']:.3e}, "
            f"elapsed={elapsed_seconds:.1f}s",
            flush=True,
        )

    if not settings:
        raise FermiRouteConvergenceError(
            f"All npoints values failed for {ion['label']} skin={ion['neutron_skin_fm']:.2f} fm."
        )

    reference = max(settings, key=lambda item: item["npoints"])
    for setting in settings:
        setting["reference_npoints"] = reference["npoints"]
        setting["rel_diff_Delta_E_vs_highest_npoints"] = relative_difference(
            reference["Delta_E_W_mu1"], setting["Delta_E_W_mu1"]
        )
        setting["rel_diff_D_W_vs_highest_npoints"] = relative_difference(
            reference["D_W_h_5e-4"], setting["D_W_h_5e-4"]
        )
        setting["rel_diff_delta_g_W_vs_highest_npoints"] = relative_difference(
            reference["delta_g_W_h_5e-4"], setting["delta_g_W_h_5e-4"]
        )

    decision = convergence_decision(settings, reference["npoints"])
    return {
        "isotope_label": ion["label"],
        "Z": ion["znucl"],
        "N": ion["nneutrons"],
        "proton_rms_fm": ion["proton_rms_fm"],
        "neutron_rms_fm": ion["neutron_rms_fm"],
        "neutron_skin_fm": ion["neutron_skin_fm"],
        "proton_a_fm": ion["proton_a_fm"],
        "neutron_a_fm": ion["neutron_a_fm"],
        "highest_available_npoints": reference["npoints"],
        "decision": decision,
        "settings": settings,
        "failures": failures,
    }


def convergence_decision(settings: list[dict[str, Any]], highest_npoints: int) -> dict[str, Any]:
    target = next((item for item in settings if item["npoints"] == REFERENCE_NPOINTS_FOR_DECISION), None)
    if target is None:
        return {
            "reference_npoints_for_decision": REFERENCE_NPOINTS_FOR_DECISION,
            "highest_available_npoints": highest_npoints,
            "rel_diff_delta_g_W_at_1000": None,
            "verdict": "npoints_1000_unavailable",
        }

    rel_diff = target["rel_diff_delta_g_W_vs_highest_npoints"]
    if rel_diff < 1.0e-6:
        verdict = "npoints_1000_remains_frozen_for_final_internal_route"
    elif rel_diff < 1.0e-5:
        verdict = "npoints_1000_acceptable_internal_consider_1200_or_1500_for_final_tables"
    else:
        verdict = "recommend_increasing_default_npoints"
    return {
        "reference_npoints_for_decision": REFERENCE_NPOINTS_FOR_DECISION,
        "highest_available_npoints": highest_npoints,
        "rel_diff_delta_g_W_at_1000": rel_diff,
        "verdict": verdict,
    }


def overall_decision(case_records: list[dict[str, Any]]) -> dict[str, Any]:
    rel_values = [
        case["decision"]["rel_diff_delta_g_W_at_1000"]
        for case in case_records
        if case["decision"]["rel_diff_delta_g_W_at_1000"] is not None
    ]
    if len(rel_values) != len(case_records):
        verdict = "convergence_remains_unresolved_and_needs_debugging"
    else:
        max_rel = max(rel_values)
        if max_rel < 1.0e-6:
            verdict = "npoints_1000_remains_frozen_for_final_internal_route"
        elif max_rel < 1.0e-5:
            verdict = "npoints_1000_acceptable_internal_consider_1200_or_1500_for_final_tables"
        else:
            verdict = "npoints_should_be_increased_for_final_route_tables"
    return {
        "max_rel_diff_delta_g_W_at_1000": max(rel_values) if rel_values else None,
        "verdict": verdict,
    }


def write_json(
    case_records: list[dict[str, Any]],
    output_path: Path,
    import_route: str,
    module: Any,
    npoints_values: list[int],
) -> None:
    payload = {
        "benchmark_name": "fermi_route_convergence",
        "physics_status": "final_internal_route_convergence_only_no_production_delta_g_claim",
        "route": {
            "energy_route": "direct-DKB/intlw",
            "coulomb_model": "Fermi proton charge distribution",
            "weak_model": "separate proton/neutron Fermi weak density",
            "state": "H-like 1s",
            "order": 9,
            "rcavity": "2.5/alpha",
            "quadrature": "nodes20/weights20",
            "h_reference": 5.0e-4,
            "h_stability_check": 1.0e-3,
            "external_mu2_weak_normalization_scaling": "retained",
            "weak_convention": "unchanged from Sections 31-37",
        },
        "npoints_values_requested": npoints_values,
        "npoints_decision_rule": {
            "keep_1000": "rel_diff_delta_g_W at npoints=1000 below 1e-6 versus highest available npoints",
            "consider_larger_for_final_tables": "between 1e-6 and 1e-5",
            "increase_default": "above 1e-5",
        },
        "import_route": import_route,
        "module_file": getattr(module, "__file__", "n/a"),
        "overall_decision": overall_decision(case_records),
        "cases": case_records,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(case_records: list[dict[str, Any]], output_path: Path) -> None:
    fieldnames = [
        "isotope_label",
        "Z",
        "N",
        "neutron_skin_fm",
        "proton_rms_fm",
        "neutron_rms_fm",
        "proton_a_fm",
        "neutron_a_fm",
        "npoints",
        "reference_npoints",
        "order",
        "nknots",
        "nbspl",
        "nspectrum",
        "rcavity",
        "state_index",
        "energy_total",
        "binding_energy",
        "Delta_E_W_mu1",
        "D_W_h_5e-4",
        "delta_g_W_h_5e-4",
        "D_W_h_1e-3",
        "delta_g_W_h_1e-3",
        "step_size_relative_difference_delta_g_W",
        "rel_diff_Delta_E_vs_highest_npoints",
        "rel_diff_D_W_vs_highest_npoints",
        "rel_diff_delta_g_W_vs_highest_npoints",
        "proton_normalization_check_mu1",
        "neutron_normalization_check_mu1",
        "elapsed_seconds",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in case_records:
            for setting in case["settings"]:
                writer.writerow(
                    {
                        "isotope_label": case["isotope_label"],
                        "Z": case["Z"],
                        "N": case["N"],
                        "neutron_skin_fm": f"{case['neutron_skin_fm']:.10f}",
                        "proton_rms_fm": f"{case['proton_rms_fm']:.10f}",
                        "neutron_rms_fm": f"{case['neutron_rms_fm']:.10f}",
                        "proton_a_fm": f"{case['proton_a_fm']:.16e}",
                        "neutron_a_fm": f"{case['neutron_a_fm']:.16e}",
                        "npoints": setting["npoints"],
                        "reference_npoints": setting["reference_npoints"],
                        "order": setting["order"],
                        "nknots": setting["nknots"],
                        "nbspl": setting["nbspl"],
                        "nspectrum": setting["nspectrum"],
                        "rcavity": f"{setting['rcavity']:.16e}",
                        "state_index": setting["state_index"],
                        "energy_total": f"{setting['energy_total']:.16e}",
                        "binding_energy": f"{setting['binding_energy']:.16e}",
                        "Delta_E_W_mu1": f"{setting['Delta_E_W_mu1']:.16e}",
                        "D_W_h_5e-4": f"{setting['D_W_h_5e-4']:.16e}",
                        "delta_g_W_h_5e-4": f"{setting['delta_g_W_h_5e-4']:.16e}",
                        "D_W_h_1e-3": f"{setting['D_W_h_1e-3']:.16e}",
                        "delta_g_W_h_1e-3": f"{setting['delta_g_W_h_1e-3']:.16e}",
                        "step_size_relative_difference_delta_g_W": (
                            f"{setting['step_size_relative_difference_delta_g_W']:.16e}"
                        ),
                        "rel_diff_Delta_E_vs_highest_npoints": (
                            f"{setting['rel_diff_Delta_E_vs_highest_npoints']:.16e}"
                        ),
                        "rel_diff_D_W_vs_highest_npoints": (
                            f"{setting['rel_diff_D_W_vs_highest_npoints']:.16e}"
                        ),
                        "rel_diff_delta_g_W_vs_highest_npoints": (
                            f"{setting['rel_diff_delta_g_W_vs_highest_npoints']:.16e}"
                        ),
                        "proton_normalization_check_mu1": (
                            f"{setting['proton_normalization_check_mu1']:.16e}"
                        ),
                        "neutron_normalization_check_mu1": (
                            f"{setting['neutron_normalization_check_mu1']:.16e}"
                        ),
                        "elapsed_seconds": f"{setting['elapsed_seconds']:.6f}",
                    }
                )


def main() -> int:
    args = parse_args()
    npoints_values = validate_npoints_values(args.npoints)
    output_dir = args.output_dir.expanduser()
    json_output = output_dir / JSON_OUTPUT.name
    csv_output = output_dir / CSV_OUTPUT.name

    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise FermiRouteConvergenceError("Imported module does not expose `rad_kernel_f2py`.")
    if not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_pn_fermi_fermi_coulomb"):
        raise FermiRouteConvergenceError(
            "Imported f2py module does not expose `direct_dkb_weak_shift_1s_pn_fermi_fermi_coulomb`; rebuild rad_f2py."
        )

    print("Final-route convergence: Fermi Coulomb plus p/n Fermi weak density")
    print(f"Import route: {import_route}")
    print(f"Module file: {getattr(module, '__file__', 'n/a')}")
    print(f"npoints sweep: {', '.join(str(value) for value in npoints_values)}")

    sampled_reference = load_sampled_reference()
    case_records = [evaluate_case(module, sampled_reference, ion, npoints_values) for ion in ion_cases()]

    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(case_records, json_output, import_route, module, npoints_values)
    write_csv(case_records, csv_output)

    decision = overall_decision(case_records)
    print("")
    print(f"Overall decision: {decision['verdict']}")
    print(f"Max rel_diff_delta_g_W at npoints=1000: {decision['max_rel_diff_delta_g_W_at_1000']:.16e}")
    print(f"JSON output: {json_output}")
    print(f"CSV output: {csv_output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FermiRouteConvergenceError,
        DirectDkbMassDerivativeError,
        FirstOrderError,
        MassDerivativeError,
        subprocess.SubprocessError,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
