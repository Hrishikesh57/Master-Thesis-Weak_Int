#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

from First_order import FirstOrderError, import_rad_kernel_suite
from Mass_derivative import MassDerivativeError, relative_difference
from Mass_derivative_direct_dkb import (
    DirectDkbMassDerivativeError,
    FORTRAN_ROUTINE_NAME,
    FORTRAN_ROUTINE_NPOINTS_NAME,
    ION_SET,
    OUTPUT_DIR,
    WEAK_CONVENTION_LABEL,
    evaluate_step,
)


JSON_OUTPUT = OUTPUT_DIR / "direct_dkb_convergence.json"
CSV_OUTPUT = OUTPUT_DIR / "direct_dkb_convergence.csv"
DEFAULT_NPOINTS_VALUES = [600, 800, 1000, 1200, 1500]
DEFAULT_H = 5.0e-4
SECTION_29_REFERENCES = {
    "Pb-208": {
        "delta_epsilon_mu1": -3.8780004284383617e-12,
        "derivative_dw": -9.9691587855362991e-12,
        "delta_g_w": -1.3292211714048398e-11,
    },
    "U-238": {
        "delta_epsilon_mu1": -9.5696639911281104e-12,
        "derivative_dw": -2.3349791620499020e-11,
        "delta_g_w": -3.1133055493998691e-11,
    },
}


class DirectDkbConvergenceError(RuntimeError):
    """User-facing direct-DKB convergence error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convergence sweep for the direct DKB/intlw weak-energy mass-derivative route."
    )
    parser.add_argument(
        "--npoints",
        type=int,
        nargs="+",
        default=DEFAULT_NPOINTS_VALUES,
        help="B-spline grid point counts to sweep. Values must be positive multiples of 10.",
    )
    parser.add_argument(
        "--h",
        type=float,
        default=DEFAULT_H,
        help="Finite-difference step for the mass derivative. Default: 5e-4.",
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
            raise DirectDkbConvergenceError(f"npoints={value} is too small for the DKB knot construction.")
        if value % 10 != 0:
            raise DirectDkbConvergenceError(
                f"npoints={value} is not divisible by 10; this sweep keeps ninside=npoints/10 exact."
            )
    return unique_values


def state_check(mu_zero: dict[str, Any]) -> str:
    energy = mu_zero["energy_total"]
    if 0.0 < energy < 1.0 and mu_zero["binding_energy"] < 0.0:
        return "positive_bound_1s_slot"
    return "check_state_index_or_spectrum"


def evaluate_ion(module: Any, ion: dict[str, Any], npoints_values: list[int], h_value: float) -> dict[str, Any]:
    settings: list[dict[str, Any]] = []
    for npoints in npoints_values:
        cache: dict[str, dict[str, Any]] = {}
        step = evaluate_step(module, ion, h_value, npoints, cache)
        mu_zero = step["mu_zero"]
        settings.append(
            {
                "npoints": npoints,
                "ninside": npoints // 10,
                "nknots": mu_zero["nknots"],
                "nbspl": mu_zero["nbspl"],
                "nspectrum": mu_zero["nspectrum"],
                "order": mu_zero["order"],
                "rcavity": mu_zero["rcavity"],
                "state_index": mu_zero["state_index"],
                "energy_total": mu_zero["energy_total"],
                "binding_energy": mu_zero["binding_energy"],
                "state_check": state_check(mu_zero),
                "delta_epsilon_mu1": mu_zero["delta_epsilon"],
                "raw_direct_dkb_delta_epsilon_mu1": mu_zero["raw_direct_dkb_delta_epsilon"],
                "derivative_dw": step["derivative_dw"],
                "delta_g_w": step["delta_g_w"],
                "h": h_value,
                "mu_points": {
                    "mu_minus": step["mu_minus"],
                    "mu_zero": step["mu_zero"],
                    "mu_plus": step["mu_plus"],
                },
            }
        )

    reference = max(settings, key=lambda item: item["npoints"])
    for setting in settings:
        setting["abs_diff_delta_e_vs_highest_npoints"] = (
            setting["delta_epsilon_mu1"] - reference["delta_epsilon_mu1"]
        )
        setting["rel_diff_delta_e_vs_highest_npoints"] = relative_difference(
            setting["delta_epsilon_mu1"], reference["delta_epsilon_mu1"]
        )
        setting["abs_diff_derivative_dw_vs_highest_npoints"] = (
            setting["derivative_dw"] - reference["derivative_dw"]
        )
        setting["rel_diff_derivative_dw_vs_highest_npoints"] = relative_difference(
            setting["derivative_dw"], reference["derivative_dw"]
        )
        setting["abs_diff_delta_g_vs_highest_npoints"] = setting["delta_g_w"] - reference["delta_g_w"]
        setting["rel_diff_delta_g_vs_highest_npoints"] = relative_difference(
            setting["delta_g_w"], reference["delta_g_w"]
        )

        section_29 = SECTION_29_REFERENCES.get(ion["label"])
        if section_29 is not None and setting["npoints"] == 1000:
            setting["section_29_default_check"] = {
                "abs_diff_delta_epsilon_mu1": setting["delta_epsilon_mu1"] - section_29["delta_epsilon_mu1"],
                "rel_diff_delta_epsilon_mu1": relative_difference(
                    setting["delta_epsilon_mu1"], section_29["delta_epsilon_mu1"]
                ),
                "abs_diff_derivative_dw": setting["derivative_dw"] - section_29["derivative_dw"],
                "rel_diff_derivative_dw": relative_difference(setting["derivative_dw"], section_29["derivative_dw"]),
                "abs_diff_delta_g_w": setting["delta_g_w"] - section_29["delta_g_w"],
                "rel_diff_delta_g_w": relative_difference(setting["delta_g_w"], section_29["delta_g_w"]),
            }

    return {
        "isotope_label": ion["label"],
        "Z": ion["znucl"],
        "N": ion["nneutrons"],
        "state": "1s",
        "physical_rms_fm": ion["rms_fm"],
        "physical_common_sphere_radius_fm": math.sqrt(5.0 / 3.0) * ion["rms_fm"],
        "energy_route": "direct-dkb",
        "radius_scaling_mode": "fixed-physical",
        "weak_convention_label": WEAK_CONVENTION_LABEL,
        "reference_npoints": reference["npoints"],
        "settings": settings,
    }


def write_json(records: list[dict[str, Any]], output_path: Path, import_route: str, module: Any, h_value: float) -> None:
    payload = {
        "benchmark_name": "direct_dkb_convergence",
        "physics_status": "internal_common_radius_prototype_only",
        "weak_convention_label": WEAK_CONVENTION_LABEL,
        "workflow": {
            "energy_route": "direct-dkb",
            "fortran_routine_default": FORTRAN_ROUTINE_NAME,
            "fortran_routine_variable_npoints": FORTRAN_ROUTINE_NPOINTS_NAME,
            "import_route": import_route,
            "module_file": getattr(module, "__file__", "n/a"),
            "h": h_value,
            "fixed_parameters": {
                "order": 9,
                "rcavity_rule": "2.5 / alpha",
                "quadrature": "nodes20 / weights20",
                "state": "1s, kappa=-1",
                "state_index_rule": "nspectrum/2 + 1",
                "ninside_rule": "npoints/10",
            },
        },
        "ions": records,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(records: list[dict[str, Any]], output_path: Path) -> None:
    fieldnames = [
        "isotope_label",
        "Z",
        "N",
        "physical_rms_fm",
        "npoints",
        "ninside",
        "order",
        "nknots",
        "nbspl",
        "nspectrum",
        "rcavity",
        "state_index",
        "energy_total",
        "binding_energy",
        "delta_epsilon_mu1",
        "derivative_dw",
        "delta_g_w",
        "rel_diff_delta_e_vs_highest_npoints",
        "rel_diff_delta_g_vs_highest_npoints",
        "state_check",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            for setting in record["settings"]:
                writer.writerow(
                    {
                        "isotope_label": record["isotope_label"],
                        "Z": record["Z"],
                        "N": record["N"],
                        "physical_rms_fm": f"{record['physical_rms_fm']:.10f}",
                        "npoints": setting["npoints"],
                        "ninside": setting["ninside"],
                        "order": setting["order"],
                        "nknots": setting["nknots"],
                        "nbspl": setting["nbspl"],
                        "nspectrum": setting["nspectrum"],
                        "rcavity": f"{setting['rcavity']:.16e}",
                        "state_index": setting["state_index"],
                        "energy_total": f"{setting['energy_total']:.16e}",
                        "binding_energy": f"{setting['binding_energy']:.16e}",
                        "delta_epsilon_mu1": f"{setting['delta_epsilon_mu1']:.16e}",
                        "derivative_dw": f"{setting['derivative_dw']:.16e}",
                        "delta_g_w": f"{setting['delta_g_w']:.16e}",
                        "rel_diff_delta_e_vs_highest_npoints": (
                            f"{setting['rel_diff_delta_e_vs_highest_npoints']:.16e}"
                        ),
                        "rel_diff_delta_g_vs_highest_npoints": (
                            f"{setting['rel_diff_delta_g_vs_highest_npoints']:.16e}"
                        ),
                        "state_check": setting["state_check"],
                    }
                )


def print_record(record: dict[str, Any]) -> None:
    print("")
    print(f"{record['isotope_label']} direct-DKB npoints convergence")
    print(f"Z = {record['Z']}, N = {record['N']}, rms_fm = {record['physical_rms_fm']:.10f}")
    print(f"Reference npoints for relative changes: {record['reference_npoints']}")
    print(
        "npoints   Delta_E_W(mu=1)        D_W(h=5e-4)           delta_g_W"
        "             rel_dE_ref       rel_dg_ref       E_total"
    )
    for setting in record["settings"]:
        print(
            f"{setting['npoints']:7d}  "
            f"{setting['delta_epsilon_mu1']: .16e}  "
            f"{setting['derivative_dw']: .16e}  "
            f"{setting['delta_g_w']: .16e}  "
            f"{setting['rel_diff_delta_e_vs_highest_npoints']: .3e}  "
            f"{setting['rel_diff_delta_g_vs_highest_npoints']: .3e}  "
            f"{setting['energy_total']: .12e}"
        )
        if setting.get("section_29_default_check"):
            check = setting["section_29_default_check"]
            print(
                "          Section 29 npoints=1000 check: "
                f"rel_dE={check['rel_diff_delta_epsilon_mu1']:.3e}, "
                f"rel_dg={check['rel_diff_delta_g_w']:.3e}"
            )


def main() -> int:
    args = parse_args()
    npoints_values = validate_npoints_values(args.npoints)
    output_dir = args.output_dir.expanduser()
    json_output = output_dir / "direct_dkb_convergence.json"
    csv_output = output_dir / "direct_dkb_convergence.csv"

    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise DirectDkbConvergenceError("Imported module does not expose `rad_kernel_f2py`.")
    if not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_npoints"):
        raise DirectDkbConvergenceError("Imported f2py module does not expose `direct_dkb_weak_shift_1s_npoints`; rebuild rad_f2py.")

    print("Direct-DKB weak mass-derivative convergence sweep")
    print("Physics status: internal common-radius prototype only")
    print(f"Import route: {import_route}")
    print(f"Module file: {getattr(module, '__file__', 'n/a')}")
    print(f"Fortran routine: {FORTRAN_ROUTINE_NPOINTS_NAME}")
    print(f"h = {args.h:.6e}")
    print(f"npoints values: {', '.join(str(value) for value in npoints_values)}")

    records = [evaluate_ion(module, ion, npoints_values, args.h) for ion in ION_SET]
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(records, json_output, import_route, module, args.h)
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
    except (
        DirectDkbConvergenceError,
        DirectDkbMassDerivativeError,
        FirstOrderError,
        MassDerivativeError,
        subprocess.SubprocessError,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
