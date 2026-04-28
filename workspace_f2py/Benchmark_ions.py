#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

from First_order import FORTRAN_ROUTINE_NAME, FirstOrderError, import_rad_kernel_suite
from Mass_derivative import (
    MassDerivativeError,
    StepResult,
    ensure_runtime_dependencies,
    evaluate_step_size,
    preserve_standard_wavefunction,
    relative_difference,
    resolve_physical_radii,
    standard_wavefunction_path,
    validate_h_values,
)


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = WORKSPACE_F2PY_ROOT / "outputs" / "weak_g_benchmark"
JSON_OUTPUT = OUTPUT_DIR / "two_ion_benchmark.json"
CSV_OUTPUT = OUTPUT_DIR / "two_ion_benchmark.csv"
WEAK_CONVENTION_LABEL = "internal_common_radius_Cve_Cvp_1_minus_4sin2_Qimpl_CvpZ_minus_N"
DEFAULT_H_VALUES = [1.0e-3, 5.0e-4]
ION_SET = [
    {"label": "Pb-208", "znucl": 82, "nneutrons": 126, "rms_fm": 5.5012},
    {"label": "U-238", "znucl": 92, "nneutrons": 146, "rms_fm": 5.8571},
]


class BenchmarkError(RuntimeError):
    """User-facing two-ion benchmark error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the two-ion sampled-rad weak-g benchmark using the existing "
            "Mass_derivative.py prototype workflow."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for JSON and CSV benchmark outputs.",
    )
    return parser.parse_args()


def mass_derivative_args(ion: dict[str, Any]) -> argparse.Namespace:
    return argparse.Namespace(
        znucl=ion["znucl"],
        nneutrons=ion["nneutrons"],
        rms_fm=ion["rms_fm"],
        radius_fm=None,
        h_values=list(DEFAULT_H_VALUES),
        wavefunction_mode="regenerate",
        baseline_wavefunction=None,
        radius_scaling_mode="fixed-physical",
        allow_unnormalized_diagnostic=False,
    )


def path_to_string(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: path_to_string(item) for key, item in value.items()}
    if isinstance(value, list):
        return [path_to_string(item) for item in value]
    return value


def step_to_record(step_result: StepResult) -> dict[str, Any]:
    return path_to_string(asdict(step_result))


def stability_verdict(relative_difference_value: float) -> str:
    if relative_difference_value < 1.0e-6:
        return "stable_to_better_than_1e-6"
    if relative_difference_value < 1.0e-4:
        return "stable_to_better_than_1e-4"
    if relative_difference_value < 1.0e-2:
        return "stable_to_percent_level"
    return "not_stable_at_percent_level"


def summarize_steps(step_results: list[StepResult]) -> dict[str, Any]:
    if not step_results:
        raise BenchmarkError("No step results were produced.")

    by_h = {f"{step.h:.6e}": step for step in step_results}
    smallest_h = min(step_results, key=lambda item: item.h)
    largest_h = max(step_results, key=lambda item: item.h)
    derivative_rel_diff = relative_difference(smallest_h.derivative_dw, largest_h.derivative_dw)
    g_rel_diff = relative_difference(smallest_h.delta_g_w, largest_h.delta_g_w)

    return {
        "delta_epsilon_mu1": smallest_h.mu_zero.delta_epsilon,
        "delta_epsilon_mu1_h": smallest_h.h,
        "derivative_dw_by_h": {h_key: step.derivative_dw for h_key, step in by_h.items()},
        "delta_g_w_by_h": {h_key: step.delta_g_w for h_key, step in by_h.items()},
        "reference_h": smallest_h.h,
        "derivative_dw_reference": smallest_h.derivative_dw,
        "delta_g_w_reference": smallest_h.delta_g_w,
        "step_size_relative_difference_dw": derivative_rel_diff,
        "step_size_relative_difference_delta_g": g_rel_diff,
        "stability_verdict": stability_verdict(derivative_rel_diff),
        "steps": [step_to_record(step) for step in sorted(step_results, key=lambda item: item.h)],
    }


def run_ion(module: Any, ion: dict[str, Any]) -> dict[str, Any]:
    args = mass_derivative_args(ion)
    h_values = validate_h_values(args.h_values)
    physical_rms_fm, physical_radius_fm, radius_input_mode = resolve_physical_radii(args)
    standard_path = standard_wavefunction_path(args.znucl)
    step_results: list[StepResult] = []

    with preserve_standard_wavefunction(standard_path):
        for h_value in h_values:
            step_results.append(
                evaluate_step_size(
                    module,
                    args,
                    h_value,
                    physical_rms_fm,
                    physical_radius_fm,
                    args.radius_scaling_mode,
                    None,
                )
            )

    summary = summarize_steps(step_results)
    return {
        "isotope_label": ion["label"],
        "Z": args.znucl,
        "N": args.nneutrons,
        "rms_fm": physical_rms_fm,
        "common_sphere_radius_fm": physical_radius_fm,
        "radius_input_mode": radius_input_mode,
        "h_values": h_values,
        "wavefunction_mode": args.wavefunction_mode,
        "radius_scaling_mode": args.radius_scaling_mode,
        "weak_convention_label": WEAK_CONVENTION_LABEL,
        **summary,
    }


def write_json(records: list[dict[str, Any]], output_path: Path, import_route: str, module: Any) -> None:
    payload = {
        "benchmark_name": "two_ion_sampled_rad_common_radius",
        "physics_status": "internal_common_radius_prototype_only",
        "weak_convention_label": WEAK_CONVENTION_LABEL,
        "weak_convention": {
            "Cv_e": "1 - 4 sin^2(theta_W)",
            "Cv_p": "1 - 4 sin^2(theta_W)",
            "Q_impl": "Cv_p * Z - N",
            "rho0": "hard_sphere_density * Q_impl * G_F / (2 sqrt(2))",
            "V_W_pc": "Cv_e * rho",
        },
        "workflow": {
            "driver": str(Path(__file__).resolve()),
            "mass_derivative_module": str((WORKSPACE_F2PY_ROOT / "Mass_derivative.py").resolve()),
            "fortran_routine": FORTRAN_ROUTINE_NAME,
            "import_route": import_route,
            "module_file": getattr(module, "__file__", "n/a"),
        },
        "ions": records,
    }
    output_path.write_text(json.dumps(path_to_string(payload), indent=2, sort_keys=True) + "\n")


def write_csv(records: list[dict[str, Any]], output_path: Path) -> None:
    fieldnames = [
        "isotope_label",
        "Z",
        "N",
        "rms_fm",
        "common_sphere_radius_fm",
        "delta_epsilon_mu1",
        "derivative_dw_h_1e-3",
        "delta_g_w_h_1e-3",
        "derivative_dw_h_5e-4",
        "delta_g_w_h_5e-4",
        "step_size_relative_difference_dw",
        "stability_verdict",
        "wavefunction_mode",
        "radius_scaling_mode",
        "weak_convention_label",
    ]

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "isotope_label": record["isotope_label"],
                    "Z": record["Z"],
                    "N": record["N"],
                    "rms_fm": f"{record['rms_fm']:.10f}",
                    "common_sphere_radius_fm": f"{record['common_sphere_radius_fm']:.10f}",
                    "delta_epsilon_mu1": f"{record['delta_epsilon_mu1']:.16e}",
                    "derivative_dw_h_1e-3": f"{record['derivative_dw_by_h']['1.000000e-03']:.16e}",
                    "delta_g_w_h_1e-3": f"{record['delta_g_w_by_h']['1.000000e-03']:.16e}",
                    "derivative_dw_h_5e-4": f"{record['derivative_dw_by_h']['5.000000e-04']:.16e}",
                    "delta_g_w_h_5e-4": f"{record['delta_g_w_by_h']['5.000000e-04']:.16e}",
                    "step_size_relative_difference_dw": f"{record['step_size_relative_difference_dw']:.16e}",
                    "stability_verdict": record["stability_verdict"],
                    "wavefunction_mode": record["wavefunction_mode"],
                    "radius_scaling_mode": record["radius_scaling_mode"],
                    "weak_convention_label": record["weak_convention_label"],
                }
            )


def print_record(record: dict[str, Any]) -> None:
    print("")
    print(f"{record['isotope_label']} sampled-rad common-radius benchmark")
    print(f"Z = {record['Z']}, N = {record['N']}, rms_fm = {record['rms_fm']:.10f}")
    print(f"Common sphere radius R_fm = {record['common_sphere_radius_fm']:.10f}")
    print(f"Delta_epsilon_W(mu=1) = {record['delta_epsilon_mu1']:.16e}")
    for h_key in sorted(record["derivative_dw_by_h"]):
        print(
            f"h = {h_key}: D_W = {record['derivative_dw_by_h'][h_key]:.16e}, "
            f"delta_g_W = {record['delta_g_w_by_h'][h_key]:.16e}"
        )
    print(f"Relative h-difference in D_W: {record['step_size_relative_difference_dw']:.16e}")
    print(f"Stability verdict: {record['stability_verdict']}")


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.expanduser()
    json_output = output_dir / "two_ion_benchmark.json"
    csv_output = output_dir / "two_ion_benchmark.csv"

    ensure_runtime_dependencies()
    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise BenchmarkError("Imported module does not expose `rad_kernel_f2py`.")

    print("Two-ion sampled-rad common-radius benchmark")
    print("Physics status: internal prototype only; no production delta_g_W claim")
    print(f"Weak convention label: {WEAK_CONVENTION_LABEL}")
    print(f"Import route: {import_route}")

    records = [run_ion(module, ion) for ion in ION_SET]
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(records, json_output, import_route, module)
    write_csv(records, csv_output)

    for record in records:
        print_record(record)

    print("")
    print(f"JSON output: {json_output}")
    print(f"CSV output: {csv_output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (BenchmarkError, FirstOrderError, MassDerivativeError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
