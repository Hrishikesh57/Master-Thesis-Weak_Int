#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from First_order import FirstOrderError, import_rad_kernel_suite
from Mass_derivative import MassDerivativeError
from Mass_derivative_direct_dkb import (
    DEFAULT_FERMI_A_FM,
    DirectDkbMassDerivativeError,
    OUTPUT_DIR,
    fermi_c_from_rms,
    load_sampled_reference,
    summarize_ion,
)


JSON_OUTPUT = OUTPUT_DIR / "final_pb_u_neutron_skin_table.json"
CSV_OUTPUT = OUTPUT_DIR / "final_pb_u_neutron_skin_table.csv"
NPOINTS = 1000
NEUTRON_SKIN_VALUES_FM = [0.0, 0.05, 0.10, 0.20]
SECTION_38_SKIN0_CHECKS = {
    "Pb-208": -1.3328870635281932e-11,
    "U-238": -3.1233277218115705e-11,
}

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


class FinalNeutronSkinTableError(RuntimeError):
    """User-facing final internal table error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Final internal Pb/U neutron-skin table for Fermi Coulomb plus p/n Fermi weak density."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for JSON and CSV final-table outputs.",
    )
    return parser.parse_args()


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
        "ion": record["isotope_label"],
        "Z": record["Z"],
        "N": record["N"],
        "proton_rms_fm": record["physical_proton_rms_fm"],
        "neutron_rms_fm": record["physical_neutron_rms_fm"],
        "neutron_skin_fm": record["neutron_skin_fm"],
        "proton_a_fm": record["physical_proton_a_fm"],
        "neutron_a_fm": record["physical_neutron_a_fm"],
        "proton_c_fm": record["physical_proton_c_fm"],
        "neutron_c_fm": record["physical_neutron_c_fm"],
        "Delta_E_W_mu1": record["delta_epsilon_mu1"],
        "D_W_h_5e-4": record["derivative_dw_reference"],
        "delta_g_W_h_5e-4": record["delta_g_w_reference"],
        "D_W_h_1e-3": record["derivative_dw_by_h"]["1.000000e-03"],
        "delta_g_W_h_1e-3": record["delta_g_w_by_h"]["1.000000e-03"],
        "step_size_relative_difference": record["step_size_relative_difference_delta_g"],
        "npoints": record["npoints_requested"],
        "order": mu_zero["order"],
        "rcavity": mu_zero["rcavity"],
        "quadrature": "nodes20/weights20",
        "state": "H-like 1s",
        "proton_normalization_check_mu1": record["proton_normalization_check_mu1"],
        "neutron_normalization_check_mu1": record["neutron_normalization_check_mu1"],
        "fermi_coulomb_v_cp_mu1": record["fermi_coulomb_v_cp_mu1"],
        "elapsed_seconds": elapsed_seconds,
    }


def add_baseline_changes(rows: list[dict[str, Any]]) -> None:
    baselines: dict[str, float] = {}
    for row in rows:
        if abs(row["neutron_skin_fm"]) < 1.0e-14:
            baselines[row["ion"]] = row["delta_g_W_h_5e-4"]

    for row in rows:
        baseline = baselines.get(row["ion"])
        if baseline in (None, 0.0):
            raise FinalNeutronSkinTableError(f"Missing nonzero skin=0 baseline for {row['ion']}.")
        row["signed_relative_change_vs_skin0"] = (row["delta_g_W_h_5e-4"] - baseline) / baseline
        row["relative_magnitude_change_vs_skin0"] = (
            abs(row["delta_g_W_h_5e-4"]) - abs(baseline)
        ) / abs(baseline)


def check_section_38_reproduction(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    for row in rows:
        if abs(row["neutron_skin_fm"]) >= 1.0e-14:
            continue
        expected = SECTION_38_SKIN0_CHECKS[row["ion"]]
        actual = row["delta_g_W_h_5e-4"]
        abs_diff = actual - expected
        rel_diff = abs(abs_diff) / abs(expected)
        checks.append(
            {
                "ion": row["ion"],
                "section_38_delta_g_W": expected,
                "final_table_delta_g_W": actual,
                "absolute_difference": abs_diff,
                "relative_difference": rel_diff,
            }
        )
        if rel_diff > 1.0e-14:
            raise FinalNeutronSkinTableError(
                f"{row['ion']} skin=0 row does not reproduce Section 38 within roundoff: rel_diff={rel_diff:.3e}."
            )
    return checks


def build_rows(module: Any, sampled_reference: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for ion in ion_cases():
        print(f"Running {ion['label']} skin={ion['neutron_skin_fm']:.2f} fm...", flush=True)
        start = time.monotonic()
        record = summarize_ion(module, ion, sampled_reference, NPOINTS, "pn-fermi", "fermi")
        row = compact_record(record, time.monotonic() - start)
        rows.append(row)
        print(
            f"  delta_g_W={row['delta_g_W_h_5e-4']:.16e}, "
            f"step_rel={row['step_size_relative_difference']:.3e}",
            flush=True,
        )
    add_baseline_changes(rows)
    return rows


def write_json(
    rows: list[dict[str, Any]],
    reproduction_checks: list[dict[str, Any]],
    output_path: Path,
    import_route: str,
    module: Any,
) -> None:
    payload = {
        "benchmark_name": "final_pb_u_neutron_skin_table",
        "physics_status": "final_internal_numerical_result_set_no_experimental_observability_claim",
        "route": {
            "energy_route": "direct-DKB/intlw",
            "coulomb_model": "Fermi proton charge distribution",
            "weak_model": "separate proton/neutron Fermi weak density",
            "state": "H-like 1s",
            "order": 9,
            "npoints": NPOINTS,
            "rcavity": "2.5/alpha",
            "quadrature": "nodes20/weights20",
            "h_reference": 5.0e-4,
            "h_stability_check": 1.0e-3,
            "external_mu2_weak_normalization_scaling": "retained",
            "weak_convention": "frozen internal convention",
        },
        "inputs": {
            "neutron_skin_values_fm": NEUTRON_SKIN_VALUES_FM,
            "default_a_fm": DEFAULT_FERMI_A_FM,
        },
        "import_route": import_route,
        "module_file": getattr(module, "__file__", "n/a"),
        "section_38_reproduction_checks": reproduction_checks,
        "rows": rows,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(rows: list[dict[str, Any]], output_path: Path) -> None:
    fieldnames = [
        "ion",
        "Z",
        "N",
        "proton_rms_fm",
        "neutron_rms_fm",
        "neutron_skin_fm",
        "proton_a_fm",
        "neutron_a_fm",
        "Delta_E_W_mu1",
        "D_W_h_5e-4",
        "delta_g_W_h_5e-4",
        "signed_relative_change_vs_skin0",
        "relative_magnitude_change_vs_skin0",
        "step_size_relative_difference",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "ion": row["ion"],
                    "Z": row["Z"],
                    "N": row["N"],
                    "proton_rms_fm": f"{row['proton_rms_fm']:.10f}",
                    "neutron_rms_fm": f"{row['neutron_rms_fm']:.10f}",
                    "neutron_skin_fm": f"{row['neutron_skin_fm']:.10f}",
                    "proton_a_fm": f"{row['proton_a_fm']:.16e}",
                    "neutron_a_fm": f"{row['neutron_a_fm']:.16e}",
                    "Delta_E_W_mu1": f"{row['Delta_E_W_mu1']:.16e}",
                    "D_W_h_5e-4": f"{row['D_W_h_5e-4']:.16e}",
                    "delta_g_W_h_5e-4": f"{row['delta_g_W_h_5e-4']:.16e}",
                    "signed_relative_change_vs_skin0": f"{row['signed_relative_change_vs_skin0']:.16e}",
                    "relative_magnitude_change_vs_skin0": f"{row['relative_magnitude_change_vs_skin0']:.16e}",
                    "step_size_relative_difference": f"{row['step_size_relative_difference']:.16e}",
                }
            )


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.expanduser()
    json_output = output_dir / JSON_OUTPUT.name
    csv_output = output_dir / CSV_OUTPUT.name

    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise FinalNeutronSkinTableError("Imported module does not expose `rad_kernel_f2py`.")
    if not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_pn_fermi_fermi_coulomb"):
        raise FinalNeutronSkinTableError(
            "Imported f2py module does not expose `direct_dkb_weak_shift_1s_pn_fermi_fermi_coulomb`; rebuild rad_f2py."
        )

    print("Final internal Pb/U p/n Fermi neutron-skin table")
    print(f"Import route: {import_route}")
    print(f"Module file: {getattr(module, '__file__', 'n/a')}")
    print(f"npoints: {NPOINTS}")

    sampled_reference = load_sampled_reference()
    rows = build_rows(module, sampled_reference)
    reproduction_checks = check_section_38_reproduction(rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(rows, reproduction_checks, json_output, import_route, module)
    write_csv(rows, csv_output)

    print("")
    for check in reproduction_checks:
        print(
            f"Section 38 check {check['ion']}: rel_diff={check['relative_difference']:.3e}, "
            f"abs_diff={check['absolute_difference']:.3e}"
        )
    print(f"JSON output: {json_output}")
    print(f"CSV output: {csv_output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        FinalNeutronSkinTableError,
        DirectDkbMassDerivativeError,
        FirstOrderError,
        MassDerivativeError,
        subprocess.SubprocessError,
    ) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
