#!/usr/bin/env python3
"""Controlled scalar-potential sign check for the weak-g mass derivative.

This diagnostic compares the repository's direct magnetic-response convention
against a mass-derivative prediction for a non-weak central scalar test
potential. It does not use weak constants and does not update production
tables.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from First_order import import_rad_kernel_suite


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "outputs" / "weak_g_benchmark"
JSON_OUT = OUTPUT_DIR / "mass_derivative_sign_check.json"
CSV_OUT = OUTPUT_DIR / "mass_derivative_sign_check.csv"

ION = {
    "label": "U-238",
    "Z": 92,
    "rms_fm": 5.8571,
}

NPOINTS = 200
TEST_R0_FM = 5.0
TEST_LAMBDA_REL_ABS = 1.0e-6
H_VALUES = (1.0e-3, 5.0e-4)


class MassDerivativeSignCheckError(RuntimeError):
    """Raised when the diagnostic wrapper is unavailable or returns an error."""


def relative_difference(reference: float, candidate: float) -> float | None:
    if reference == 0.0:
        return None
    return abs(candidate - reference) / abs(reference)


def call_response(
    module: Any,
    *,
    znucl: int,
    rms_fm: float,
    test_r0_fm: float,
    test_lambda_rel: float,
    npoints: int,
) -> dict[str, Any]:
    routine = module.rad_kernel_f2py.direct_dkb_exponential_test_response_1s
    result = routine(znucl, rms_fm, test_r0_fm, test_lambda_rel, npoints)
    (
        delta_g_direct,
        delta_e_test,
        test_r0_rel,
        nuclear_radius_fm,
        nuclear_radius_rel,
        order,
        npoints_out,
        nknots,
        nbspl,
        nspectrum,
        rcavity,
        state_index,
        energy_total,
        binding_energy,
        ierr,
    ) = result
    if int(ierr) != 0:
        raise MassDerivativeSignCheckError(
            "direct_dkb_exponential_test_response_1s returned "
            f"ierr={int(ierr)} for Z={znucl}, rms_fm={rms_fm}, "
            f"test_r0_fm={test_r0_fm}, lambda={test_lambda_rel}, npoints={npoints}"
        )
    return {
        "delta_g_direct": float(delta_g_direct),
        "delta_e_test": float(delta_e_test),
        "test_r0_rel": float(test_r0_rel),
        "nuclear_radius_fm": float(nuclear_radius_fm),
        "nuclear_radius_rel": float(nuclear_radius_rel),
        "order": int(order),
        "npoints": int(npoints_out),
        "nknots": int(nknots),
        "nbspl": int(nbspl),
        "nspectrum": int(nspectrum),
        "rcavity": float(rcavity),
        "state_index": int(state_index),
        "energy_total": float(energy_total),
        "binding_energy": float(binding_energy),
    }


def evaluate_mu_point(module: Any, lambda_sign: int, mu: float) -> dict[str, Any]:
    # Fixed physical test potential convention:
    #   R0_code(mu) = mu R0_code(1)
    #   lambda_code(mu) = lambda_code(1) / mu
    # so no weak-like mu^2 scaling is applied.
    lambda_reference = lambda_sign * TEST_LAMBDA_REL_ABS
    return call_response(
        module,
        znucl=ION["Z"],
        rms_fm=mu * ION["rms_fm"],
        test_r0_fm=mu * TEST_R0_FM,
        test_lambda_rel=lambda_reference / mu,
        npoints=NPOINTS,
    )


def build_rows(module: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for lambda_sign in (1, -1):
        lambda_reference = lambda_sign * TEST_LAMBDA_REL_ABS
        zero = evaluate_mu_point(module, lambda_sign, 1.0)
        for h_value in H_VALUES:
            minus = evaluate_mu_point(module, lambda_sign, 1.0 - h_value)
            plus = evaluate_mu_point(module, lambda_sign, 1.0 + h_value)
            d_test = (
                (1.0 + h_value) * plus["delta_e_test"]
                - (1.0 - h_value) * minus["delta_e_test"]
            ) / (2.0 * h_value)
            delta_g_md_plus = (4.0 / 3.0) * d_test
            delta_g_md_minus = -delta_g_md_plus
            rel_plus = relative_difference(zero["delta_g_direct"], delta_g_md_plus)
            rel_minus = relative_difference(zero["delta_g_direct"], delta_g_md_minus)
            rows.append(
                {
                    "ion": ION["label"],
                    "Z": ION["Z"],
                    "rms_fm": ION["rms_fm"],
                    "npoints": NPOINTS,
                    "h": h_value,
                    "lambda_sign": lambda_sign,
                    "lambda_rel_mu1": lambda_reference,
                    "test_r0_fm_mu1": TEST_R0_FM,
                    "test_r0_rel_mu1": zero["test_r0_rel"],
                    "Delta_E_test_mu1": zero["delta_e_test"],
                    "D_test": d_test,
                    "delta_g_direct": zero["delta_g_direct"],
                    "delta_g_MD_plus": delta_g_md_plus,
                    "delta_g_MD_minus": delta_g_md_minus,
                    "relative_difference_direct_vs_plus": rel_plus,
                    "relative_difference_direct_vs_minus": rel_minus,
                    "energy_total_mu1": zero["energy_total"],
                    "binding_energy_mu1": zero["binding_energy"],
                    "state_index": zero["state_index"],
                    "nspectrum": zero["nspectrum"],
                    "rcavity": zero["rcavity"],
                }
            )
    return rows


def write_outputs(rows: list[dict[str, Any]], import_route: str) -> dict[str, Any]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "diagnostic": "mass_derivative_sign_check",
        "purpose": (
            "Compare direct magnetic response to +/- (4/3) mass-derivative "
            "predictions for a controlled non-weak scalar test potential."
        ),
        "potential": {
            "formula": "V_test(r) = lambda_code(mu) * exp[-r / R0_code(mu)]",
            "fixed_physical_convention": (
                "R0_code(mu)=mu*R0_code(1), lambda_code(mu)=lambda_code(1)/mu; "
                "no weak mu^2 scaling is applied."
            ),
            "test_r0_fm_mu1": TEST_R0_FM,
            "lambda_rel_abs_mu1": TEST_LAMBDA_REL_ABS,
        },
        "ion": ION,
        "npoints": NPOINTS,
        "h_values": list(H_VALUES),
        "f2py_import_source": import_route,
        "rows": rows,
    }
    JSON_OUT.write_text(json.dumps(payload, indent=2))

    fieldnames = [
        "ion",
        "Z",
        "rms_fm",
        "npoints",
        "h",
        "lambda_sign",
        "lambda_rel_mu1",
        "test_r0_fm_mu1",
        "test_r0_rel_mu1",
        "Delta_E_test_mu1",
        "D_test",
        "delta_g_direct",
        "delta_g_MD_plus",
        "delta_g_MD_minus",
        "relative_difference_direct_vs_plus",
        "relative_difference_direct_vs_minus",
        "energy_total_mu1",
        "binding_energy_mu1",
        "state_index",
        "nspectrum",
        "rcavity",
    ]
    with CSV_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        ""
                        if row.get(key) is None
                        else f"{row[key]:.16e}"
                        if isinstance(row.get(key), float)
                        else row.get(key)
                    )
                    for key in fieldnames
                }
            )
    return payload


def main() -> None:
    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise MassDerivativeSignCheckError("Imported module does not expose `rad_kernel_f2py`.")
    if not hasattr(module.rad_kernel_f2py, "direct_dkb_exponential_test_response_1s"):
        raise MassDerivativeSignCheckError(
            "Imported f2py module does not expose `direct_dkb_exponential_test_response_1s`; rebuild rad_f2py."
        )

    rows = build_rows(module)
    write_outputs(rows, import_route)

    print(f"Saved {JSON_OUT}")
    print(f"Saved {CSV_OUT}")
    for row in rows:
        print(
            f"lambda_sign={row['lambda_sign']:+d}, h={row['h']:.1e}, "
            f"direct={row['delta_g_direct']:.16e}, "
            f"+MD={row['delta_g_MD_plus']:.16e}, "
            f"-MD={row['delta_g_MD_minus']:.16e}, "
            f"rel_plus={row['relative_difference_direct_vs_plus']:.3e}, "
            f"rel_minus={row['relative_difference_direct_vs_minus']:.3e}"
        )


if __name__ == "__main__":
    main()
