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

from First_order import FirstOrderError, import_rad_kernel_suite, normalization_trapz, read_sampled_wavefunction
from Mass_derivative import (
    MassDerivativeError,
    preserve_standard_wavefunction,
    run_practise_for_1s,
    standard_wavefunction_path,
)


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = WORKSPACE_F2PY_ROOT / "outputs" / "weak_g_benchmark"
WEAK_CONVENTION_LABEL = "internal_common_radius_Cve_Cvp_1_minus_4sin2_Qimpl_CvpZ_minus_N"
ION_SET = [
    {"label": "Pb-208", "znucl": 82, "nneutrons": 126, "rms_fm": 5.5012},
    {"label": "U-238", "znucl": 92, "nneutrons": 146, "rms_fm": 5.8571},
]


class DirectDkbComparisonError(RuntimeError):
    """User-facing direct-DKB comparison error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare direct coefficient-basis DKB weak-energy integrals against "
            "the existing exported sampled-rad weak-energy path."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for JSON and CSV comparison outputs.",
    )
    return parser.parse_args()


def direct_dkb_record(module: Any, ion: dict[str, Any]) -> dict[str, Any]:
    routine = module.rad_kernel_f2py.direct_dkb_weak_shift_1s
    (
        delta_e_direct,
        radius_fm,
        r_rel,
        order,
        npoints,
        nknots,
        nbspl,
        nspectrum,
        rcavity,
        state_index,
        energy_total,
        binding_energy,
        qweak,
        rho0,
        pc_half_radius,
        ierr,
    ) = routine(ion["znucl"], ion["nneutrons"], ion["rms_fm"])
    if ierr != 0:
        raise DirectDkbComparisonError(f"direct_dkb_weak_shift_1s returned ierr={int(ierr)} for {ion['label']}")

    return {
        "delta_e_direct_dkb": float(delta_e_direct),
        "common_sphere_radius_fm": float(radius_fm),
        "r_rel": float(r_rel),
        "order": int(order),
        "npoints": int(npoints),
        "nknots": int(nknots),
        "nbspl": int(nbspl),
        "nspectrum": int(nspectrum),
        "rcavity": float(rcavity),
        "state_index": int(state_index),
        "energy_total": float(energy_total),
        "binding_energy": float(binding_energy),
        "q_impl": float(qweak),
        "rho0": float(rho0),
        "v_w_pc_half_radius": float(pc_half_radius),
    }


def sampled_record(module: Any, ion: dict[str, Any], radius_fm: float) -> dict[str, Any]:
    standard_path = standard_wavefunction_path(ion["znucl"])
    generated_path = run_practise_for_1s(ion["znucl"], ion["rms_fm"])
    if generated_path != standard_path:
        raise DirectDkbComparisonError(
            f"Unexpected sampled wavefunction path for {ion['label']}: {generated_path}"
        )

    r_array, g_array, f_array = read_sampled_wavefunction(generated_path)
    normalization = normalization_trapz(r_array, g_array, f_array)
    delta_e_sampled, r_rel_sampled, ierr = module.rad_kernel_f2py.first_order_weak_shift(
        ion["znucl"], ion["nneutrons"], radius_fm, r_array, g_array, f_array
    )
    if ierr != 0:
        raise DirectDkbComparisonError(f"first_order_weak_shift returned ierr={int(ierr)} for {ion['label']}")

    return {
        "delta_e_sampled_rad": float(delta_e_sampled),
        "sampled_r_rel": float(r_rel_sampled),
        "sampled_wavefunction_path": str(generated_path),
        "sampled_normalization_trapz": float(normalization),
        "sampled_npoints": int(r_array.size),
        "sampled_r_min": float(r_array[0]),
        "sampled_r_second": float(r_array[1]) if r_array.size > 1 else math.nan,
        "sampled_r_max": float(r_array[-1]),
        "sampled_points_inside_R": int((r_array <= r_rel_sampled).sum()),
    }


def compare_ion(module: Any, ion: dict[str, Any]) -> dict[str, Any]:
    direct = direct_dkb_record(module, ion)
    standard_path = standard_wavefunction_path(ion["znucl"])
    with preserve_standard_wavefunction(standard_path):
        sampled = sampled_record(module, ion, direct["common_sphere_radius_fm"])

    abs_diff = direct["delta_e_direct_dkb"] - sampled["delta_e_sampled_rad"]
    rel_diff_direct = abs(abs_diff) / abs(direct["delta_e_direct_dkb"])
    rel_diff_sampled = abs(abs_diff) / abs(sampled["delta_e_sampled_rad"])

    return {
        "isotope_label": ion["label"],
        "Z": ion["znucl"],
        "N": ion["nneutrons"],
        "state": "1s",
        "rms_fm": ion["rms_fm"],
        "weak_convention_label": WEAK_CONVENTION_LABEL,
        **direct,
        **sampled,
        "absolute_difference_direct_minus_sampled": abs_diff,
        "relative_difference_vs_direct": rel_diff_direct,
        "relative_difference_vs_sampled": rel_diff_sampled,
        "diagnosis": (
            "sampled_rad_under_resolves_nuclear_region"
            if sampled["sampled_points_inside_R"] <= 1
            else "requires_further_diagnosis"
        ),
    }


def write_json(records: list[dict[str, Any]], output_path: Path, import_route: str, module: Any) -> None:
    payload = {
        "comparison_name": "direct_dkb_vs_sampled_weak_energy",
        "physics_status": "representation_cross_check_only",
        "weak_convention_label": WEAK_CONVENTION_LABEL,
        "workflow": {
            "direct_dkb": "solve -> selected coefficients -> bspqed_integrals::intlw",
            "sampled_rad": "practise solve -> exported r,G,F table -> first_order_weak_shift",
            "import_route": import_route,
            "module_file": getattr(module, "__file__", "n/a"),
        },
        "ions": records,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(records: list[dict[str, Any]], output_path: Path) -> None:
    fieldnames = [
        "isotope_label",
        "Z",
        "N",
        "state",
        "rms_fm",
        "common_sphere_radius_fm",
        "r_rel",
        "energy_total",
        "binding_energy",
        "q_impl",
        "rho0",
        "v_w_pc_half_radius",
        "delta_e_direct_dkb",
        "delta_e_sampled_rad",
        "absolute_difference_direct_minus_sampled",
        "relative_difference_vs_direct",
        "relative_difference_vs_sampled",
        "sampled_normalization_trapz",
        "sampled_npoints",
        "sampled_r_min",
        "sampled_r_second",
        "sampled_r_max",
        "sampled_points_inside_R",
        "diagnosis",
        "weak_convention_label",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            row = {key: record[key] for key in fieldnames}
            writer.writerow(row)


def print_record(record: dict[str, Any]) -> None:
    print("")
    print(f"{record['isotope_label']} direct-DKB versus sampled-rad weak-energy check")
    print(f"Z = {record['Z']}, N = {record['N']}, state = {record['state']}, rms_fm = {record['rms_fm']:.10f}")
    print(f"Common weak sphere radius R_fm = {record['common_sphere_radius_fm']:.10f}")
    print(f"R_rel = {record['r_rel']:.16e}, rcavity = {record['rcavity']:.10f}")
    print(f"order = {record['order']}, npoints = {record['npoints']}, nknots = {record['nknots']}")
    print(f"nbspl = {record['nbspl']}, nspectrum = {record['nspectrum']}, state_index = {record['state_index']}")
    print(f"energy_total = {record['energy_total']:.16e}, binding_energy = {record['binding_energy']:.16e}")
    print(f"Q_impl = {record['q_impl']:.16e}, rho0 = {record['rho0']:.16e}")
    print(f"V_W^pc(0.5R) = {record['v_w_pc_half_radius']:.16e}")
    print(f"Delta_E_W^direct_DKB = {record['delta_e_direct_dkb']:.16e}")
    print(f"Delta_E_W^sampled_rad = {record['delta_e_sampled_rad']:.16e}")
    print(f"direct - sampled = {record['absolute_difference_direct_minus_sampled']:.16e}")
    print(f"relative difference vs direct = {record['relative_difference_vs_direct']:.16e}")
    print(
        "sampled grid: "
        f"n={record['sampled_npoints']}, r_min={record['sampled_r_min']:.3e}, "
        f"r_second={record['sampled_r_second']:.3e}, points_inside_R={record['sampled_points_inside_R']}"
    )
    print(f"diagnosis: {record['diagnosis']}")


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.expanduser()
    json_output = output_dir / "direct_dkb_vs_sampled_weak_energy.json"
    csv_output = output_dir / "direct_dkb_vs_sampled_weak_energy.csv"

    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise DirectDkbComparisonError("Imported module does not expose `rad_kernel_f2py`.")
    if not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s"):
        raise DirectDkbComparisonError("Imported f2py module does not expose `direct_dkb_weak_shift_1s`; rebuild rad_f2py.")

    print("Direct DKB versus sampled-rad weak-energy comparison")
    print("Physics status: representation cross-check only; no production delta_g_W claim")
    print(f"Weak convention label: {WEAK_CONVENTION_LABEL}")
    print(f"Import route: {import_route}")

    records = [compare_ion(module, ion) for ion in ION_SET]
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
    except (DirectDkbComparisonError, FirstOrderError, MassDerivativeError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
