#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad

from First_order import import_rad_kernel_suite


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = WORKSPACE_F2PY_ROOT / "outputs" / "weak_g_benchmark"
JSON_OUTPUT = OUTPUT_DIR / "fermi_coulomb_validation.json"
CSV_OUTPUT = OUTPUT_DIR / "fermi_coulomb_validation.csv"

IALPHA = 137.035999139
ALPHA = 1.0 / IALPHA
HBARC_MEV_FM = 197.3269788
ELECTRON_REST_ENERGY_MEV = 0.5109989461
ELECTRON_REDUCED_COMPTON_FM = HBARC_MEV_FM / ELECTRON_REST_ENERGY_MEV

DEFAULT_FERMI_A_FM = 2.3 / (4.0 * math.log(3.0))
IONS = [
    {"isotope": "Pb-208", "Z": 82, "rms_fm": 5.5012},
    {"isotope": "U-238", "Z": 92, "rms_fm": 5.8571},
]


class FermiCoulombValidationError(RuntimeError):
    """User-facing validation error."""


def fermi_c_from_rms(rms_fm: float, a_fm: float) -> float:
    c_sq = (5.0 * rms_fm * rms_fm - 7.0 * (math.pi * a_fm) ** 2) / 3.0
    if c_sq <= 0.0:
        raise FermiCoulombValidationError("Fermi c_fm is not real for the requested rms_fm and a_fm.")
    return math.sqrt(c_sq)


def fermi_shape(r_rel: float, c_rel: float, a_rel: float) -> float:
    arg = (r_rel - c_rel) / a_rel
    if arg > 80.0:
        return math.exp(-arg)
    if arg < -80.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(arg))


def integrate_shape(c_rel: float, a_rel: float, power: int, lower: float, upper: float) -> float:
    def integrand(r_rel: float) -> float:
        return (r_rel**power) * fermi_shape(r_rel, c_rel, a_rel)

    value, error = quad(integrand, lower, upper, epsabs=1.0e-16, epsrel=1.0e-12, limit=300)
    if not math.isfinite(value):
        raise FermiCoulombValidationError("Non-finite numerical Fermi integral.")
    if abs(error) > max(1.0e-14, 1.0e-9 * abs(value)):
        raise FermiCoulombValidationError(
            f"Numerical Fermi integral did not converge tightly enough: value={value:.16e}, error={error:.3e}"
        )
    return float(value)


def numerical_reference(
    znucl: int,
    c_rel: float,
    a_rel: float,
    r_values: np.ndarray,
) -> tuple[float, float, list[float]]:
    rmax = max(20.0 * c_rel, c_rel + 100.0 * a_rel, 120.0 * a_rel)
    shape_norm = 4.0 * math.pi * integrate_shape(c_rel, a_rel, 2, 0.0, rmax)
    rho0 = 1.0 / shape_norm
    normalization_check = rho0 * shape_norm

    potentials: list[float] = []
    for r_rel in r_values:
        if r_rel == 0.0:
            central = 4.0 * math.pi * rho0 * integrate_shape(c_rel, a_rel, 1, 0.0, rmax)
            potentials.append(-float(znucl) * ALPHA * central)
            continue

        q_inside = 4.0 * math.pi * rho0 * integrate_shape(c_rel, a_rel, 2, 0.0, float(r_rel))
        tail_integral = 4.0 * math.pi * rho0 * integrate_shape(c_rel, a_rel, 1, float(r_rel), rmax)
        potentials.append(-float(znucl) * ALPHA * (q_inside / float(r_rel) + tail_integral))

    return rho0, normalization_check, potentials


def radius_grid(c_rel: float, a_rel: float) -> list[tuple[str, float]]:
    return [
        ("0", 0.0),
        ("1e-8", 1.0e-8),
        ("1e-6", 1.0e-6),
        ("0.01c", 0.01 * c_rel),
        ("0.1c", 0.1 * c_rel),
        ("c", c_rel),
        ("c+5a", c_rel + 5.0 * a_rel),
        ("2c", 2.0 * c_rel),
        ("5c", 5.0 * c_rel),
        ("10c", 10.0 * c_rel),
    ]


def analytic_values(module: Any, znucl: int, rms_fm: float, a_fm: float, r_values: np.ndarray) -> dict[str, Any]:
    if not hasattr(module.rad_kernel_f2py, "fermi_coulomb_eval_points"):
        raise FermiCoulombValidationError(
            "Imported f2py module does not expose `fermi_coulomb_eval_points`; rebuild rad_f2py."
        )

    result = module.rad_kernel_f2py.fermi_coulomb_eval_points(znucl, rms_fm, a_fm, r_values)
    c_fm, a_fm_used, c_rel, a_rel, v_values, ierr = result
    if int(ierr) != 0:
        raise FermiCoulombValidationError(f"fermi_coulomb_eval_points returned ierr={int(ierr)}.")
    return {
        "c_fm": float(c_fm),
        "a_fm": float(a_fm_used),
        "c_rel": float(c_rel),
        "a_rel": float(a_rel),
        "v_values": [float(value) for value in v_values],
    }


def tail_ratio(znucl: int, r_rel: float, potential: float) -> float | None:
    if r_rel == 0.0:
        return None
    return potential / (-float(znucl) * ALPHA / r_rel)


def build_validation() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    module, import_source = import_rad_kernel_suite()
    json_records: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []

    for ion in IONS:
        isotope = ion["isotope"]
        znucl = int(ion["Z"])
        rms_fm = float(ion["rms_fm"])
        a_fm = DEFAULT_FERMI_A_FM
        c_fm = fermi_c_from_rms(rms_fm, a_fm)
        c_rel = c_fm / ELECTRON_REDUCED_COMPTON_FM
        a_rel = a_fm / ELECTRON_REDUCED_COMPTON_FM
        labels_and_radii = radius_grid(c_rel, a_rel)
        r_values = np.asarray([r_rel for _, r_rel in labels_and_radii], dtype=np.float64)

        analytic = analytic_values(module, znucl, rms_fm, a_fm, r_values)
        rho0, normalization_check, v_num = numerical_reference(znucl, c_rel, a_rel, r_values)

        if abs(analytic["c_rel"] - c_rel) > 1.0e-15 or abs(analytic["a_rel"] - a_rel) > 1.0e-15:
            raise FermiCoulombValidationError("Python and Fortran Fermi parameter conversions disagree.")

        point_records: list[dict[str, Any]] = []
        for (label, r_rel), v_num_i, v_analytic_i in zip(labels_and_radii, v_num, analytic["v_values"]):
            abs_diff = abs(v_analytic_i - v_num_i)
            rel_diff = abs_diff / abs(v_num_i) if v_num_i != 0.0 else None
            point = {
                "radius_label": label,
                "r_rel": r_rel,
                "V_num": v_num_i,
                "V_analytic": v_analytic_i,
                "absolute_difference": abs_diff,
                "relative_difference": rel_diff,
                "tail_ratio_num": tail_ratio(znucl, r_rel, v_num_i),
                "tail_ratio_analytic": tail_ratio(znucl, r_rel, v_analytic_i),
            }
            point_records.append(point)
            csv_rows.append(
                {
                    "isotope": isotope,
                    "Z": znucl,
                    "rms_fm": rms_fm,
                    "a_fm": a_fm,
                    "c_fm": c_fm,
                    "c_rel": c_rel,
                    "a_rel": a_rel,
                    "rho0": rho0,
                    "normalization_check": normalization_check,
                    **point,
                }
            )

        json_records.append(
            {
                "isotope": isotope,
                "Z": znucl,
                "rms_fm": rms_fm,
                "a_fm": a_fm,
                "c_fm": c_fm,
                "c_rel": c_rel,
                "a_rel": a_rel,
                "rho0": rho0,
                "normalization_check": normalization_check,
                "points": point_records,
            }
        )

    payload = {
        "description": "Independent normalized-Fermi-density validation of bspqed analytic Fermi Coulomb potential.",
        "units": "relativistic units, hbar = c = m_e = 1; radii are divided by electron reduced Compton wavelength.",
        "constants": {
            "alpha": ALPHA,
            "electron_reduced_compton_fm": ELECTRON_REDUCED_COMPTON_FM,
            "a_fm": DEFAULT_FERMI_A_FM,
        },
        "analytic_source": "rad_kernel_suite.rad_kernel_f2py.fermi_coulomb_eval_points -> bspqed_potentials::fermi_potential_eval",
        "f2py_import_source": import_source,
        "records": json_records,
    }
    return payload, csv_rows


def write_outputs(payload: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    JSON_OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    fieldnames = [
        "isotope",
        "Z",
        "rms_fm",
        "a_fm",
        "c_fm",
        "c_rel",
        "a_rel",
        "rho0",
        "normalization_check",
        "radius_label",
        "r_rel",
        "V_num",
        "V_analytic",
        "absolute_difference",
        "relative_difference",
        "tail_ratio_num",
        "tail_ratio_analytic",
    ]
    with CSV_OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    payload, rows = build_validation()
    write_outputs(payload, rows)

    print(f"Saved JSON: {JSON_OUTPUT}")
    print(f"Saved CSV: {CSV_OUTPUT}")
    for record in payload["records"]:
        print(
            f"{record['isotope']}: c_fm={record['c_fm']:.10f}, c_rel={record['c_rel']:.16e}, "
            f"rho0={record['rho0']:.16e}, norm={record['normalization_check']:.16e}"
        )
        for point in record["points"]:
            if point["radius_label"] in {"0", "1e-8", "1e-6", "c", "c+5a", "10c"}:
                rel = point["relative_difference"]
                rel_text = "" if rel is None else f", rel={rel:.3e}"
                print(
                    f"  {point['radius_label']:>4s}: V_num={point['V_num']:.16e}, "
                    f"V_analytic={point['V_analytic']:.16e}, "
                    f"abs={point['absolute_difference']:.3e}{rel_text}"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
