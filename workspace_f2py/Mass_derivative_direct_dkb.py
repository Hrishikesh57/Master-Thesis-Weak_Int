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


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = WORKSPACE_F2PY_ROOT / "outputs" / "weak_g_benchmark"
JSON_OUTPUT = OUTPUT_DIR / "direct_dkb_mass_derivative.json"
CSV_OUTPUT = OUTPUT_DIR / "direct_dkb_mass_derivative.csv"
SAMPLED_BENCHMARK_JSON = OUTPUT_DIR / "two_ion_benchmark.json"
FORTRAN_ROUTINE_NAME = "rad_kernel_suite.rad_kernel_f2py.direct_dkb_weak_shift_1s"
FORTRAN_ROUTINE_NPOINTS_NAME = "rad_kernel_suite.rad_kernel_f2py.direct_dkb_weak_shift_1s_npoints"
FORTRAN_ROUTINE_PN_SPHERES_NAME = "rad_kernel_suite.rad_kernel_f2py.direct_dkb_weak_shift_1s_pn_spheres"
FORTRAN_ROUTINE_COMMON_FERMI_NAME = "rad_kernel_suite.rad_kernel_f2py.direct_dkb_weak_shift_1s_common_fermi"
FORTRAN_ROUTINE_PN_FERMI_NAME = "rad_kernel_suite.rad_kernel_f2py.direct_dkb_weak_shift_1s_pn_fermi"
FORTRAN_ROUTINE_PN_FERMI_FERMI_COULOMB_NAME = (
    "rad_kernel_suite.rad_kernel_f2py.direct_dkb_weak_shift_1s_pn_fermi_fermi_coulomb"
)
WEAK_CONVENTION_LABEL = "internal_common_radius_Cve_Cvp_1_minus_4sin2_Qimpl_CvpZ_minus_N"
PN_SPHERES_CONVENTION_LABEL = "internal_pn_hard_spheres_Cve_CvpZrho_p_minus_Nrho_n"
COMMON_FERMI_CONVENTION_LABEL = "internal_common_fermi_Cve_Qimpl_CvpZ_minus_N_normalized_density"
PN_FERMI_CONVENTION_LABEL = "internal_pn_fermi_Cve_CvpZrho_p_minus_Nrho_n_normalized_densities"
DEFAULT_H_VALUES = [1.0e-3, 5.0e-4]
DEFAULT_NEUTRON_SKIN_FM_VALUES = [0.0, 0.05, 0.10, 0.20]
DEFAULT_FERMI_A_FM = 2.3 / (4.0 * math.log(3.0))
COMMON_SPHERE_DELTA_G_REFERENCE = {
    "Pb-208": -1.3292211714048398e-11,
    "U-238": -3.1133055493998691e-11,
}
COMMON_FERMI_DELTA_G_REFERENCE = {
    "Pb-208": -1.3308757876147243e-11,
    "U-238": -3.1177873673877604e-11,
}
ION_SET = [
    {"label": "Pb-208", "znucl": 82, "nneutrons": 126, "rms_fm": 5.5012},
    {"label": "U-238", "znucl": 92, "nneutrons": 146, "rms_fm": 5.8571},
]


class DirectDkbMassDerivativeError(RuntimeError):
    """User-facing direct-DKB mass-derivative error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Mass-derivative prototype using direct DKB/intlw weak energies, "
            "not exported sampled radial tables."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for JSON and CSV benchmark outputs.",
    )
    parser.add_argument(
        "--npoints",
        type=int,
        default=1000,
        help="B-spline grid point count for the direct DKB solve; default preserves the validated baseline.",
    )
    parser.add_argument(
        "--weak-model",
        choices=("common-sphere", "pn-spheres", "common-fermi", "pn-fermi"),
        default="common-sphere",
        help="Weak source model. The default preserves the frozen common-sphere benchmark route.",
    )
    parser.add_argument(
        "--coulomb-model",
        choices=("sphere", "fermi"),
        default="sphere",
        help="Finite-size Coulomb model for the direct-DKB solve. Fermi Coulomb is currently supported only with --weak-model pn-fermi.",
    )
    parser.add_argument(
        "--rms-fm",
        type=float,
        default=None,
        help="Override physical rms radius for --weak-model common-fermi. Defaults to each ion's benchmark rms radius.",
    )
    parser.add_argument(
        "--fermi-a-fm",
        type=float,
        default=None,
        help="Physical 2pF diffuseness for --weak-model common-fermi. Defaults to 2.3/(4 log 3) fm.",
    )
    parser.add_argument(
        "--proton-rms-fm",
        type=float,
        default=None,
        help="Physical proton rms radius for --weak-model pn-spheres or pn-fermi. Defaults to each ion's benchmark rms radius.",
    )
    parser.add_argument(
        "--neutron-rms-fm",
        type=float,
        default=None,
        help="Physical neutron rms radius for --weak-model pn-spheres or pn-fermi.",
    )
    parser.add_argument(
        "--neutron-skin-fm",
        type=float,
        default=None,
        help="Physical neutron skin, neutron_rms_fm - proton_rms_fm, for --weak-model pn-spheres or pn-fermi.",
    )
    parser.add_argument(
        "--proton-a-fm",
        type=float,
        default=None,
        help="Physical proton 2pF diffuseness for --weak-model pn-fermi. Defaults to 2.3/(4 log 3) fm.",
    )
    parser.add_argument(
        "--neutron-a-fm",
        type=float,
        default=None,
        help="Physical neutron 2pF diffuseness for --weak-model pn-fermi. Defaults to the proton value.",
    )
    return parser.parse_args()


def stability_verdict(relative_difference_value: float) -> str:
    if relative_difference_value < 1.0e-6:
        return "stable_to_better_than_1e-6"
    if relative_difference_value < 1.0e-4:
        return "stable_to_better_than_1e-4"
    if relative_difference_value < 1.0e-2:
        return "stable_to_percent_level"
    return "not_stable_at_percent_level"


def weak_convention_label(weak_model: str) -> str:
    if weak_model == "common-sphere":
        return WEAK_CONVENTION_LABEL
    if weak_model == "pn-spheres":
        return PN_SPHERES_CONVENTION_LABEL
    if weak_model == "common-fermi":
        return COMMON_FERMI_CONVENTION_LABEL
    if weak_model == "pn-fermi":
        return PN_FERMI_CONVENTION_LABEL
    raise DirectDkbMassDerivativeError(f"Unsupported weak model: {weak_model}")


def fortran_routine_name(weak_model: str, coulomb_model: str = "sphere") -> str:
    if weak_model == "common-sphere":
        return FORTRAN_ROUTINE_NAME
    if weak_model == "pn-spheres":
        return FORTRAN_ROUTINE_PN_SPHERES_NAME
    if weak_model == "common-fermi":
        return FORTRAN_ROUTINE_COMMON_FERMI_NAME
    if weak_model == "pn-fermi":
        if coulomb_model == "fermi":
            return FORTRAN_ROUTINE_PN_FERMI_FERMI_COULOMB_NAME
        return FORTRAN_ROUTINE_PN_FERMI_NAME
    raise DirectDkbMassDerivativeError(f"Unsupported weak model: {weak_model}")


def fermi_c_from_rms(rms_fm: float, a_fm: float) -> float:
    c_sq = (5.0 * rms_fm * rms_fm - 7.0 * (math.pi * a_fm) ** 2) / 3.0
    if c_sq <= 0.0:
        raise DirectDkbMassDerivativeError("Fermi c_fm is not real for the requested rms_fm and a_fm.")
    return math.sqrt(c_sq)


def resolve_pn_radii(
    ion: dict[str, Any],
    proton_rms_fm: float | None,
    neutron_rms_fm: float | None,
    neutron_skin_fm: float | None,
) -> tuple[float, float]:
    if neutron_rms_fm is not None and neutron_skin_fm is not None:
        raise DirectDkbMassDerivativeError("Use either --neutron-rms-fm or --neutron-skin-fm, not both.")

    proton_value = ion["rms_fm"] if proton_rms_fm is None else proton_rms_fm
    if proton_value <= 0.0:
        raise DirectDkbMassDerivativeError("proton rms radius must be positive.")

    if neutron_rms_fm is not None:
        neutron_value = neutron_rms_fm
    elif neutron_skin_fm is not None:
        neutron_value = proton_value + neutron_skin_fm
    else:
        neutron_value = ion["rms_fm"]

    if neutron_value <= 0.0:
        raise DirectDkbMassDerivativeError("neutron rms radius must be positive.")
    return proton_value, neutron_value


def direct_raw_energy(
    module: Any,
    ion: dict[str, Any],
    mu_value: float,
    npoints_requested: int = 1000,
    weak_model: str = "common-sphere",
    coulomb_model: str = "sphere",
) -> dict[str, Any]:
    if weak_model == "common-sphere":
        effective_rms_fm = mu_value * ion["rms_fm"]
        if npoints_requested == 1000:
            routine = module.rad_kernel_f2py.direct_dkb_weak_shift_1s
            routine_result = routine(ion["znucl"], ion["nneutrons"], effective_rms_fm)
        else:
            if not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_npoints"):
                raise DirectDkbMassDerivativeError(
                    "Imported f2py module does not expose `direct_dkb_weak_shift_1s_npoints`; rebuild rad_f2py."
                )
            routine = module.rad_kernel_f2py.direct_dkb_weak_shift_1s_npoints
            routine_result = routine(ion["znucl"], ion["nneutrons"], effective_rms_fm, npoints_requested)
        (
            raw_delta_epsilon,
            effective_radius_fm,
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
        ) = routine_result
        if ierr != 0:
            raise DirectDkbMassDerivativeError(
                f"direct_dkb_weak_shift_1s returned ierr={int(ierr)} for {ion['label']} "
                f"at mu={mu_value}, npoints={npoints_requested}"
            )
        radius_payload = {
            "effective_rms_fm": effective_rms_fm,
            "effective_radius_fm": float(effective_radius_fm),
            "r_rel": float(r_rel),
            "q_impl": float(qweak),
            "rho0": float(rho0),
            "v_w_pc_half_radius": float(pc_half_radius),
        }
    elif weak_model == "pn-spheres":
        effective_proton_rms_fm = mu_value * ion["proton_rms_fm"]
        effective_neutron_rms_fm = mu_value * ion["neutron_rms_fm"]
        routine = module.rad_kernel_f2py.direct_dkb_weak_shift_1s_pn_spheres
        routine_result = routine(
            ion["znucl"],
            ion["nneutrons"],
            effective_proton_rms_fm,
            effective_neutron_rms_fm,
            npoints_requested,
        )
        (
            raw_delta_epsilon,
            proton_radius_fm,
            neutron_radius_fm,
            proton_r_rel,
            neutron_r_rel,
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
            proton_rho0,
            neutron_rho0,
            pc_half_radius,
            ierr,
        ) = routine_result
        if ierr != 0:
            raise DirectDkbMassDerivativeError(
                f"direct_dkb_weak_shift_1s_pn_spheres returned ierr={int(ierr)} for {ion['label']} "
                f"at mu={mu_value}, npoints={npoints_requested}"
            )
        radius_payload = {
            "effective_proton_rms_fm": effective_proton_rms_fm,
            "effective_neutron_rms_fm": effective_neutron_rms_fm,
            "effective_proton_radius_fm": float(proton_radius_fm),
            "effective_neutron_radius_fm": float(neutron_radius_fm),
            "proton_r_rel": float(proton_r_rel),
            "neutron_r_rel": float(neutron_r_rel),
            "q_impl_common_limit": float(qweak),
            "proton_rho0": float(proton_rho0),
            "neutron_rho0": float(neutron_rho0),
            "v_w_pc_half_proton_radius": float(pc_half_radius),
        }
    elif weak_model == "common-fermi":
        effective_rms_fm = mu_value * ion["rms_fm"]
        effective_fermi_a_fm = mu_value * ion["fermi_a_fm"]
        routine = module.rad_kernel_f2py.direct_dkb_weak_shift_1s_common_fermi
        routine_result = routine(
            ion["znucl"],
            ion["nneutrons"],
            effective_rms_fm,
            effective_fermi_a_fm,
            npoints_requested,
        )
        (
            raw_delta_epsilon,
            sphere_radius_fm,
            sphere_r_rel,
            c_fm,
            a_fm_used,
            c_rel,
            a_rel,
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
            normalization_check,
            pc_origin,
            pc_c,
            pc_c_plus_5a,
            ierr,
        ) = routine_result
        if ierr != 0:
            raise DirectDkbMassDerivativeError(
                f"direct_dkb_weak_shift_1s_common_fermi returned ierr={int(ierr)} for {ion['label']} "
                f"at mu={mu_value}, npoints={npoints_requested}"
            )
        radius_payload = {
            "effective_rms_fm": effective_rms_fm,
            "effective_fermi_a_fm": effective_fermi_a_fm,
            "effective_fermi_c_fm": float(c_fm),
            "effective_sphere_radius_fm_for_coulomb": float(sphere_radius_fm),
            "sphere_r_rel_for_coulomb": float(sphere_r_rel),
            "fermi_c_rel": float(c_rel),
            "fermi_a_rel": float(a_rel),
            "q_impl": float(qweak),
            "rho0": float(rho0),
            "normalization_check": float(normalization_check),
            "v_w_pc_origin": float(pc_origin),
            "v_w_pc_c": float(pc_c),
            "v_w_pc_c_plus_5a": float(pc_c_plus_5a),
        }
    elif weak_model == "pn-fermi":
        effective_proton_rms_fm = mu_value * ion["proton_rms_fm"]
        effective_neutron_rms_fm = mu_value * ion["neutron_rms_fm"]
        effective_proton_a_fm = mu_value * ion["proton_a_fm"]
        effective_neutron_a_fm = mu_value * ion["neutron_a_fm"]
        if coulomb_model == "fermi":
            routine = module.rad_kernel_f2py.direct_dkb_weak_shift_1s_pn_fermi_fermi_coulomb
        else:
            routine = module.rad_kernel_f2py.direct_dkb_weak_shift_1s_pn_fermi
        routine_result = routine(
            ion["znucl"],
            ion["nneutrons"],
            effective_proton_rms_fm,
            effective_neutron_rms_fm,
            effective_proton_a_fm,
            effective_neutron_a_fm,
            npoints_requested,
        )
        if coulomb_model == "fermi":
            (
                raw_delta_epsilon,
                proton_sphere_radius_fm,
                proton_sphere_r_rel,
                cp_fm,
                ap_fm_used,
                cn_fm,
                an_fm_used,
                cp_rel,
                ap_rel,
                cn_rel,
                an_rel,
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
                proton_rho0,
                neutron_rho0,
                proton_normalization_check,
                neutron_normalization_check,
                pc_origin,
                pc_cp,
                pc_cn,
                pc_cp_plus_5ap,
                pc_cn_plus_5an,
                fermi_coulomb_v_cp,
                sphere_coulomb_v_cp,
                fermi_coulomb_v_cp_plus_5ap,
                sphere_coulomb_v_cp_plus_5ap,
                ierr,
            ) = routine_result
        else:
            (
                raw_delta_epsilon,
                proton_sphere_radius_fm,
                proton_sphere_r_rel,
                cp_fm,
                ap_fm_used,
                cn_fm,
                an_fm_used,
                cp_rel,
                ap_rel,
                cn_rel,
                an_rel,
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
                proton_rho0,
                neutron_rho0,
                proton_normalization_check,
                neutron_normalization_check,
                pc_origin,
                pc_cp,
                pc_cn,
                pc_cp_plus_5ap,
                pc_cn_plus_5an,
                ierr,
            ) = routine_result
            fermi_coulomb_v_cp = None
            sphere_coulomb_v_cp = None
            fermi_coulomb_v_cp_plus_5ap = None
            sphere_coulomb_v_cp_plus_5ap = None
        if ierr != 0:
            raise DirectDkbMassDerivativeError(
                f"{fortran_routine_name(weak_model, coulomb_model)} returned ierr={int(ierr)} for {ion['label']} "
                f"at mu={mu_value}, npoints={npoints_requested}"
            )
        radius_payload = {
            "effective_proton_rms_fm": effective_proton_rms_fm,
            "effective_neutron_rms_fm": effective_neutron_rms_fm,
            "effective_proton_a_fm": effective_proton_a_fm,
            "effective_neutron_a_fm": effective_neutron_a_fm,
            "effective_proton_c_fm": float(cp_fm),
            "effective_neutron_c_fm": float(cn_fm),
            "effective_proton_sphere_radius_fm_for_coulomb": float(proton_sphere_radius_fm),
            "proton_sphere_r_rel_for_coulomb": float(proton_sphere_r_rel),
            "proton_fermi_c_rel": float(cp_rel),
            "proton_fermi_a_rel": float(ap_rel),
            "neutron_fermi_c_rel": float(cn_rel),
            "neutron_fermi_a_rel": float(an_rel),
            "q_impl_common_limit": float(qweak),
            "proton_rho0": float(proton_rho0),
            "neutron_rho0": float(neutron_rho0),
            "proton_normalization_check": float(proton_normalization_check),
            "neutron_normalization_check": float(neutron_normalization_check),
            "v_w_pc_origin": float(pc_origin),
            "v_w_pc_cp": float(pc_cp),
            "v_w_pc_cn": float(pc_cn),
            "v_w_pc_cp_plus_5ap": float(pc_cp_plus_5ap),
            "v_w_pc_cn_plus_5an": float(pc_cn_plus_5an),
        }
        if coulomb_model == "fermi":
            radius_payload.update(
                {
                    "fermi_coulomb_v_cp": float(fermi_coulomb_v_cp),
                    "sphere_coulomb_v_cp": float(sphere_coulomb_v_cp),
                    "fermi_coulomb_v_cp_plus_5ap": float(fermi_coulomb_v_cp_plus_5ap),
                    "sphere_coulomb_v_cp_plus_5ap": float(sphere_coulomb_v_cp_plus_5ap),
                }
            )
    else:
        raise DirectDkbMassDerivativeError(f"Unsupported weak model: {weak_model}")

    weak_scale_factor = mu_value * mu_value
    delta_epsilon = weak_scale_factor * float(raw_delta_epsilon)
    mu_times_delta_epsilon = mu_value * delta_epsilon

    result = {
        "mu": mu_value,
        "weak_model": weak_model,
        "coulomb_model": coulomb_model,
        "raw_direct_dkb_delta_epsilon": float(raw_delta_epsilon),
        "weak_scale_factor": weak_scale_factor,
        "delta_epsilon": delta_epsilon,
        "mu_times_delta_epsilon": mu_times_delta_epsilon,
        "order": int(order),
        "npoints": int(npoints),
        "nknots": int(nknots),
        "nbspl": int(nbspl),
        "nspectrum": int(nspectrum),
        "rcavity": float(rcavity),
        "state_index": int(state_index),
        "energy_total": float(energy_total),
        "binding_energy": float(binding_energy),
    }
    result.update(radius_payload)
    return result


def evaluate_mu_point(
    module: Any,
    ion: dict[str, Any],
    mu_value: float,
    npoints_requested: int,
    cache: dict[str, dict[str, Any]],
    weak_model: str,
    coulomb_model: str,
) -> dict[str, Any]:
    key = f"{weak_model}|coulomb={coulomb_model}|{mu_value:.10f}|npoints={npoints_requested}"
    if key not in cache:
        cache[key] = direct_raw_energy(module, ion, mu_value, npoints_requested, weak_model, coulomb_model)
    return cache[key]


def evaluate_step(
    module: Any,
    ion: dict[str, Any],
    h_value: float,
    npoints_requested: int,
    cache: dict[str, dict[str, Any]],
    weak_model: str,
    coulomb_model: str,
) -> dict[str, Any]:
    mu_minus = evaluate_mu_point(module, ion, 1.0 - h_value, npoints_requested, cache, weak_model, coulomb_model)
    mu_zero = evaluate_mu_point(module, ion, 1.0, npoints_requested, cache, weak_model, coulomb_model)
    mu_plus = evaluate_mu_point(module, ion, 1.0 + h_value, npoints_requested, cache, weak_model, coulomb_model)
    derivative_dw = (mu_plus["mu_times_delta_epsilon"] - mu_minus["mu_times_delta_epsilon"]) / (2.0 * h_value)
    delta_g_w = (4.0 / 3.0) * derivative_dw
    return {
        "h": h_value,
        "mu_minus": mu_minus,
        "mu_zero": mu_zero,
        "mu_plus": mu_plus,
        "derivative_dw": derivative_dw,
        "delta_g_w": delta_g_w,
    }


def load_sampled_reference() -> dict[str, dict[str, Any]]:
    if not SAMPLED_BENCHMARK_JSON.is_file():
        return {}
    payload = json.loads(SAMPLED_BENCHMARK_JSON.read_text())
    return {ion["isotope_label"]: ion for ion in payload.get("ions", [])}


def summarize_ion(
    module: Any,
    ion: dict[str, Any],
    sampled_reference: dict[str, dict[str, Any]],
    npoints_requested: int,
    weak_model: str,
    coulomb_model: str,
) -> dict[str, Any]:
    cache: dict[str, dict[str, Any]] = {}
    steps = [
        evaluate_step(module, ion, h_value, npoints_requested, cache, weak_model, coulomb_model)
        for h_value in DEFAULT_H_VALUES
    ]
    by_h = {f"{step['h']:.6e}": step for step in steps}
    smallest_h = min(steps, key=lambda item: item["h"])
    largest_h = max(steps, key=lambda item: item["h"])
    derivative_rel_diff = relative_difference(smallest_h["derivative_dw"], largest_h["derivative_dw"])
    g_rel_diff = relative_difference(smallest_h["delta_g_w"], largest_h["delta_g_w"])

    sampled_record = sampled_reference.get(ion["label"], {})
    sampled_delta_g_reference = sampled_record.get("delta_g_w_reference")
    sampled_derivative_reference = sampled_record.get("derivative_dw_reference")
    sampled_delta_epsilon_reference = sampled_record.get("delta_epsilon_mu1")
    direct_delta_g_reference = smallest_h["delta_g_w"]

    record = {
        "isotope_label": ion["label"],
        "Z": ion["znucl"],
        "N": ion["nneutrons"],
        "state": "1s",
        "physical_rms_fm": ion["rms_fm"],
        "energy_route": "direct-dkb",
        "npoints_requested": npoints_requested,
        "radius_scaling_mode": "fixed-physical",
        "weak_model": weak_model,
        "coulomb_model": coulomb_model,
        "weak_convention_label": weak_convention_label(weak_model),
        "h_values": DEFAULT_H_VALUES,
        "raw_direct_dkb_delta_epsilon_mu1": by_h["5.000000e-04"]["mu_zero"]["raw_direct_dkb_delta_epsilon"],
        "delta_epsilon_mu1": by_h["5.000000e-04"]["mu_zero"]["delta_epsilon"],
        "section_28_reference_match": by_h["5.000000e-04"]["mu_zero"]["raw_direct_dkb_delta_epsilon"],
        "derivative_dw_by_h": {h_key: step["derivative_dw"] for h_key, step in by_h.items()},
        "delta_g_w_by_h": {h_key: step["delta_g_w"] for h_key, step in by_h.items()},
        "reference_h": smallest_h["h"],
        "derivative_dw_reference": smallest_h["derivative_dw"],
        "delta_g_w_reference": direct_delta_g_reference,
        "step_size_relative_difference_dw": derivative_rel_diff,
        "step_size_relative_difference_delta_g": g_rel_diff,
        "stability_verdict": stability_verdict(derivative_rel_diff),
        "sampled_route_reference": {
            "delta_epsilon_mu1": sampled_delta_epsilon_reference,
            "derivative_dw_reference": sampled_derivative_reference,
            "delta_g_w_reference": sampled_delta_g_reference,
            "delta_g_ratio_direct_over_sampled": (
                direct_delta_g_reference / sampled_delta_g_reference
                if sampled_delta_g_reference not in (None, 0.0)
                else None
            ),
        },
        "steps": sorted(steps, key=lambda item: item["h"]),
    }
    if weak_model == "common-sphere":
        record["physical_common_sphere_radius_fm"] = math.sqrt(5.0 / 3.0) * ion["rms_fm"]
    elif weak_model == "pn-spheres":
        record["physical_proton_rms_fm"] = ion["proton_rms_fm"]
        record["physical_neutron_rms_fm"] = ion["neutron_rms_fm"]
        record["neutron_skin_fm"] = ion["neutron_rms_fm"] - ion["proton_rms_fm"]
        record["physical_proton_sphere_radius_fm"] = math.sqrt(5.0 / 3.0) * ion["proton_rms_fm"]
        record["physical_neutron_sphere_radius_fm"] = math.sqrt(5.0 / 3.0) * ion["neutron_rms_fm"]
    elif weak_model == "common-fermi":
        mu_zero = by_h["5.000000e-04"]["mu_zero"]
        common_sphere_reference = COMMON_SPHERE_DELTA_G_REFERENCE.get(ion["label"])
        record["physical_fermi_a_fm"] = ion["fermi_a_fm"]
        record["physical_fermi_c_fm"] = mu_zero["effective_fermi_c_fm"]
        record["physical_coulomb_sphere_radius_fm"] = mu_zero["effective_sphere_radius_fm_for_coulomb"]
        record["normalization_check_mu1"] = mu_zero["normalization_check"]
        record["rho0_mu1"] = mu_zero["rho0"]
        record["v_w_pc_origin_mu1"] = mu_zero["v_w_pc_origin"]
        record["v_w_pc_c_mu1"] = mu_zero["v_w_pc_c"]
        record["v_w_pc_c_plus_5a_mu1"] = mu_zero["v_w_pc_c_plus_5a"]
        record["common_sphere_delta_g_w_reference"] = common_sphere_reference
        record["delta_g_w_relative_change_vs_common_sphere"] = (
            (direct_delta_g_reference - common_sphere_reference) / common_sphere_reference
            if common_sphere_reference not in (None, 0.0)
            else None
        )
    elif weak_model == "pn-fermi":
        mu_zero = by_h["5.000000e-04"]["mu_zero"]
        common_fermi_reference = COMMON_FERMI_DELTA_G_REFERENCE.get(ion["label"])
        record["physical_proton_rms_fm"] = ion["proton_rms_fm"]
        record["physical_neutron_rms_fm"] = ion["neutron_rms_fm"]
        record["neutron_skin_fm"] = ion["neutron_rms_fm"] - ion["proton_rms_fm"]
        record["physical_proton_a_fm"] = ion["proton_a_fm"]
        record["physical_neutron_a_fm"] = ion["neutron_a_fm"]
        record["physical_proton_c_fm"] = mu_zero["effective_proton_c_fm"]
        record["physical_neutron_c_fm"] = mu_zero["effective_neutron_c_fm"]
        record["physical_coulomb_sphere_radius_fm"] = mu_zero["effective_proton_sphere_radius_fm_for_coulomb"]
        record["proton_rho0_mu1"] = mu_zero["proton_rho0"]
        record["neutron_rho0_mu1"] = mu_zero["neutron_rho0"]
        record["proton_normalization_check_mu1"] = mu_zero["proton_normalization_check"]
        record["neutron_normalization_check_mu1"] = mu_zero["neutron_normalization_check"]
        record["v_w_pc_origin_mu1"] = mu_zero["v_w_pc_origin"]
        record["v_w_pc_cp_mu1"] = mu_zero["v_w_pc_cp"]
        record["v_w_pc_cn_mu1"] = mu_zero["v_w_pc_cn"]
        record["v_w_pc_cp_plus_5ap_mu1"] = mu_zero["v_w_pc_cp_plus_5ap"]
        record["v_w_pc_cn_plus_5an_mu1"] = mu_zero["v_w_pc_cn_plus_5an"]
        if coulomb_model == "fermi":
            record["fermi_coulomb_v_cp_mu1"] = mu_zero["fermi_coulomb_v_cp"]
            record["sphere_coulomb_v_cp_mu1"] = mu_zero["sphere_coulomb_v_cp"]
            record["fermi_coulomb_v_cp_plus_5ap_mu1"] = mu_zero["fermi_coulomb_v_cp_plus_5ap"]
            record["sphere_coulomb_v_cp_plus_5ap_mu1"] = mu_zero["sphere_coulomb_v_cp_plus_5ap"]
        record["common_fermi_delta_g_w_reference"] = common_fermi_reference
        record["delta_g_w_relative_change_vs_common_fermi"] = (
            (direct_delta_g_reference - common_fermi_reference) / common_fermi_reference
            if common_fermi_reference not in (None, 0.0)
            else None
        )
    return record


def write_json(
    records: list[dict[str, Any]],
    output_path: Path,
    import_route: str,
    module: Any,
    weak_model: str,
    coulomb_model: str,
) -> None:
    payload = {
        "benchmark_name": (
            "direct_dkb_mass_derivative"
            if weak_model == "common-sphere"
            else "pn_sphere_validation"
            if weak_model == "pn-spheres"
            else "fermi_coulomb_comparison"
            if weak_model == "pn-fermi" and coulomb_model == "fermi"
            else "pn_fermi_validation"
            if weak_model == "pn-fermi"
            else "common_fermi_validation"
        ),
        "physics_status": (
            "internal_common_radius_prototype_only"
            if weak_model == "common-sphere"
            else "internal_pn_hard_sphere_validation_only"
            if weak_model == "pn-spheres"
            else "internal_pn_fermi_fermi_coulomb_prototype_only"
            if weak_model == "pn-fermi" and coulomb_model == "fermi"
            else "internal_pn_fermi_density_validation_only"
            if weak_model == "pn-fermi"
            else "internal_common_fermi_density_validation_only"
        ),
        "weak_model": weak_model,
        "coulomb_model": coulomb_model,
        "included_coulomb_models": sorted({record["coulomb_model"] for record in records}),
        "weak_convention_label": weak_convention_label(weak_model),
        "workflow": {
            "energy_route": "direct-dkb",
            "fortran_routine": fortran_routine_name(weak_model, coulomb_model),
            "sphere_coulomb_fortran_routine": (
                FORTRAN_ROUTINE_PN_FERMI_NAME
                if weak_model == "pn-fermi" and coulomb_model == "fermi"
                else None
            ),
            "fortran_routine_variable_npoints": FORTRAN_ROUTINE_NPOINTS_NAME,
            "import_route": import_route,
            "module_file": getattr(module, "__file__", "n/a"),
            "mass_derivative_rule": "D_W = [mu_plus Delta(mu_plus) - mu_minus Delta(mu_minus)] / (2h)",
            "g_factor_rule": "delta_g_W(1s) = (4/3) D_W",
            "weak_scaling": "Delta(mu) = mu^2 * raw_direct_DKB(mu)",
            "pn_mass_derivative_radius_rule": (
                "proton_rms_eff = mu * proton_rms_fm; neutron_rms_eff = mu * neutron_rms_fm"
                if weak_model == "pn-spheres"
                else None
            ),
            "pn_fermi_mass_derivative_radius_rule": (
                "proton_rms_eff = mu * proton_rms_fm; neutron_rms_eff = mu * neutron_rms_fm; "
                "proton_a_eff = mu * proton_a_fm; neutron_a_eff = mu * neutron_a_fm"
                if weak_model == "pn-fermi"
                else None
            ),
            "pn_fermi_coulomb_mass_derivative_rule": (
                "For --coulomb-model fermi, the Coulomb proton rms and diffuseness are the same "
                "mu-scaled proton_rms_eff and proton_a_eff used for the proton Fermi weak density."
                if weak_model == "pn-fermi" and coulomb_model == "fermi"
                else None
            ),
            "common_fermi_mass_derivative_radius_rule": (
                "rms_eff = mu * rms_fm; a_eff = mu * a_fm; c_eff recomputed from rms_eff and a_eff"
                if weak_model == "common-fermi"
                else None
            ),
            "pn_fermi_coulomb_model": (
                (
                    "finite-size Coulomb potential uses the relativistic-unit analytic proton Fermi charge potential"
                    if coulomb_model == "fermi"
                    else "finite-size Coulomb potential remains the hard sphere using the proton rms radius"
                )
                if weak_model == "pn-fermi"
                else None
            ),
            "common_fermi_coulomb_model": (
                "finite-size Coulomb potential remains the frozen hard sphere for this weak-density test"
                if weak_model == "common-fermi"
                else None
            ),
        },
        "ions": records,
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(records: list[dict[str, Any]], output_path: Path) -> None:
    fieldnames = [
        "isotope_label",
        "Z",
        "N",
        "weak_model",
        "coulomb_model",
        "physical_rms_fm",
        "physical_common_sphere_radius_fm",
        "physical_proton_rms_fm",
        "physical_neutron_rms_fm",
        "neutron_skin_fm",
        "physical_proton_sphere_radius_fm",
        "physical_neutron_sphere_radius_fm",
        "physical_fermi_a_fm",
        "physical_fermi_c_fm",
        "physical_proton_a_fm",
        "physical_neutron_a_fm",
        "physical_proton_c_fm",
        "physical_neutron_c_fm",
        "physical_coulomb_sphere_radius_fm",
        "rho0_mu1",
        "proton_rho0_mu1",
        "neutron_rho0_mu1",
        "normalization_check_mu1",
        "proton_normalization_check_mu1",
        "neutron_normalization_check_mu1",
        "v_w_pc_origin_mu1",
        "v_w_pc_c_mu1",
        "v_w_pc_c_plus_5a_mu1",
        "v_w_pc_cp_mu1",
        "v_w_pc_cn_mu1",
        "v_w_pc_cp_plus_5ap_mu1",
        "v_w_pc_cn_plus_5an_mu1",
        "fermi_coulomb_v_cp_mu1",
        "sphere_coulomb_v_cp_mu1",
        "fermi_coulomb_v_cp_plus_5ap_mu1",
        "sphere_coulomb_v_cp_plus_5ap_mu1",
        "raw_direct_dkb_delta_epsilon_mu1",
        "npoints_requested",
        "derivative_dw_h_1e-3",
        "delta_g_w_h_1e-3",
        "derivative_dw_h_5e-4",
        "delta_g_w_h_5e-4",
        "delta_g_w_relative_change_vs_skin0",
        "relative_magnitude_change_vs_skin0",
        "delta_g_w_sphere_coulomb_reference",
        "delta_g_w_absolute_difference_vs_sphere_coulomb",
        "delta_g_w_relative_difference_vs_sphere_coulomb",
        "common_sphere_delta_g_w_reference",
        "delta_g_w_relative_change_vs_common_sphere",
        "common_fermi_delta_g_w_reference",
        "delta_g_w_relative_change_vs_common_fermi",
        "step_size_relative_difference_dw",
        "stability_verdict",
        "sampled_delta_g_w_reference",
        "delta_g_ratio_direct_over_sampled",
        "energy_route",
        "radius_scaling_mode",
        "weak_convention_label",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            sampled = record["sampled_route_reference"]
            writer.writerow(
                {
                    "isotope_label": record["isotope_label"],
                    "Z": record["Z"],
                    "N": record["N"],
                    "weak_model": record["weak_model"],
                    "coulomb_model": record["coulomb_model"],
                    "physical_rms_fm": f"{record['physical_rms_fm']:.10f}",
                    "physical_common_sphere_radius_fm": (
                        f"{record['physical_common_sphere_radius_fm']:.10f}"
                        if "physical_common_sphere_radius_fm" in record
                        else ""
                    ),
                    "physical_proton_rms_fm": (
                        f"{record['physical_proton_rms_fm']:.10f}" if "physical_proton_rms_fm" in record else ""
                    ),
                    "physical_neutron_rms_fm": (
                        f"{record['physical_neutron_rms_fm']:.10f}" if "physical_neutron_rms_fm" in record else ""
                    ),
                    "neutron_skin_fm": (
                        f"{record['neutron_skin_fm']:.10f}" if "neutron_skin_fm" in record else ""
                    ),
                    "physical_proton_sphere_radius_fm": (
                        f"{record['physical_proton_sphere_radius_fm']:.10f}"
                        if "physical_proton_sphere_radius_fm" in record
                        else ""
                    ),
                    "physical_neutron_sphere_radius_fm": (
                        f"{record['physical_neutron_sphere_radius_fm']:.10f}"
                        if "physical_neutron_sphere_radius_fm" in record
                        else ""
                    ),
                    "physical_fermi_a_fm": (
                        f"{record['physical_fermi_a_fm']:.10f}" if "physical_fermi_a_fm" in record else ""
                    ),
                    "physical_fermi_c_fm": (
                        f"{record['physical_fermi_c_fm']:.10f}" if "physical_fermi_c_fm" in record else ""
                    ),
                    "physical_proton_a_fm": (
                        f"{record['physical_proton_a_fm']:.10f}" if "physical_proton_a_fm" in record else ""
                    ),
                    "physical_neutron_a_fm": (
                        f"{record['physical_neutron_a_fm']:.10f}" if "physical_neutron_a_fm" in record else ""
                    ),
                    "physical_proton_c_fm": (
                        f"{record['physical_proton_c_fm']:.10f}" if "physical_proton_c_fm" in record else ""
                    ),
                    "physical_neutron_c_fm": (
                        f"{record['physical_neutron_c_fm']:.10f}" if "physical_neutron_c_fm" in record else ""
                    ),
                    "physical_coulomb_sphere_radius_fm": (
                        f"{record['physical_coulomb_sphere_radius_fm']:.10f}"
                        if "physical_coulomb_sphere_radius_fm" in record
                        else ""
                    ),
                    "rho0_mu1": f"{record['rho0_mu1']:.16e}" if "rho0_mu1" in record else "",
                    "proton_rho0_mu1": (
                        f"{record['proton_rho0_mu1']:.16e}" if "proton_rho0_mu1" in record else ""
                    ),
                    "neutron_rho0_mu1": (
                        f"{record['neutron_rho0_mu1']:.16e}" if "neutron_rho0_mu1" in record else ""
                    ),
                    "normalization_check_mu1": (
                        f"{record['normalization_check_mu1']:.16e}" if "normalization_check_mu1" in record else ""
                    ),
                    "proton_normalization_check_mu1": (
                        f"{record['proton_normalization_check_mu1']:.16e}"
                        if "proton_normalization_check_mu1" in record
                        else ""
                    ),
                    "neutron_normalization_check_mu1": (
                        f"{record['neutron_normalization_check_mu1']:.16e}"
                        if "neutron_normalization_check_mu1" in record
                        else ""
                    ),
                    "v_w_pc_origin_mu1": (
                        f"{record['v_w_pc_origin_mu1']:.16e}" if "v_w_pc_origin_mu1" in record else ""
                    ),
                    "v_w_pc_c_mu1": f"{record['v_w_pc_c_mu1']:.16e}" if "v_w_pc_c_mu1" in record else "",
                    "v_w_pc_c_plus_5a_mu1": (
                        f"{record['v_w_pc_c_plus_5a_mu1']:.16e}"
                        if "v_w_pc_c_plus_5a_mu1" in record
                        else ""
                    ),
                    "v_w_pc_cp_mu1": f"{record['v_w_pc_cp_mu1']:.16e}" if "v_w_pc_cp_mu1" in record else "",
                    "v_w_pc_cn_mu1": f"{record['v_w_pc_cn_mu1']:.16e}" if "v_w_pc_cn_mu1" in record else "",
                    "v_w_pc_cp_plus_5ap_mu1": (
                        f"{record['v_w_pc_cp_plus_5ap_mu1']:.16e}"
                        if "v_w_pc_cp_plus_5ap_mu1" in record
                        else ""
                    ),
                    "v_w_pc_cn_plus_5an_mu1": (
                        f"{record['v_w_pc_cn_plus_5an_mu1']:.16e}"
                        if "v_w_pc_cn_plus_5an_mu1" in record
                        else ""
                    ),
                    "fermi_coulomb_v_cp_mu1": (
                        f"{record['fermi_coulomb_v_cp_mu1']:.16e}"
                        if "fermi_coulomb_v_cp_mu1" in record
                        else ""
                    ),
                    "sphere_coulomb_v_cp_mu1": (
                        f"{record['sphere_coulomb_v_cp_mu1']:.16e}"
                        if "sphere_coulomb_v_cp_mu1" in record
                        else ""
                    ),
                    "fermi_coulomb_v_cp_plus_5ap_mu1": (
                        f"{record['fermi_coulomb_v_cp_plus_5ap_mu1']:.16e}"
                        if "fermi_coulomb_v_cp_plus_5ap_mu1" in record
                        else ""
                    ),
                    "sphere_coulomb_v_cp_plus_5ap_mu1": (
                        f"{record['sphere_coulomb_v_cp_plus_5ap_mu1']:.16e}"
                        if "sphere_coulomb_v_cp_plus_5ap_mu1" in record
                        else ""
                    ),
                    "raw_direct_dkb_delta_epsilon_mu1": f"{record['raw_direct_dkb_delta_epsilon_mu1']:.16e}",
                    "npoints_requested": record["npoints_requested"],
                    "derivative_dw_h_1e-3": f"{record['derivative_dw_by_h']['1.000000e-03']:.16e}",
                    "delta_g_w_h_1e-3": f"{record['delta_g_w_by_h']['1.000000e-03']:.16e}",
                    "derivative_dw_h_5e-4": f"{record['derivative_dw_by_h']['5.000000e-04']:.16e}",
                    "delta_g_w_h_5e-4": f"{record['delta_g_w_by_h']['5.000000e-04']:.16e}",
                    "delta_g_w_relative_change_vs_skin0": (
                        f"{record['delta_g_w_relative_change_vs_skin0']:.16e}"
                        if "delta_g_w_relative_change_vs_skin0" in record
                        else ""
                    ),
                    "relative_magnitude_change_vs_skin0": (
                        f"{record['relative_magnitude_change_vs_skin0']:.16e}"
                        if "relative_magnitude_change_vs_skin0" in record
                        else ""
                    ),
                    "delta_g_w_sphere_coulomb_reference": (
                        f"{record['delta_g_w_sphere_coulomb_reference']:.16e}"
                        if "delta_g_w_sphere_coulomb_reference" in record
                        else ""
                    ),
                    "delta_g_w_absolute_difference_vs_sphere_coulomb": (
                        f"{record['delta_g_w_absolute_difference_vs_sphere_coulomb']:.16e}"
                        if "delta_g_w_absolute_difference_vs_sphere_coulomb" in record
                        else ""
                    ),
                    "delta_g_w_relative_difference_vs_sphere_coulomb": (
                        f"{record['delta_g_w_relative_difference_vs_sphere_coulomb']:.16e}"
                        if "delta_g_w_relative_difference_vs_sphere_coulomb" in record
                        else ""
                    ),
                    "common_sphere_delta_g_w_reference": (
                        f"{record['common_sphere_delta_g_w_reference']:.16e}"
                        if record.get("common_sphere_delta_g_w_reference") is not None
                        else ""
                    ),
                    "delta_g_w_relative_change_vs_common_sphere": (
                        f"{record['delta_g_w_relative_change_vs_common_sphere']:.16e}"
                        if record.get("delta_g_w_relative_change_vs_common_sphere") is not None
                        else ""
                    ),
                    "common_fermi_delta_g_w_reference": (
                        f"{record['common_fermi_delta_g_w_reference']:.16e}"
                        if record.get("common_fermi_delta_g_w_reference") is not None
                        else ""
                    ),
                    "delta_g_w_relative_change_vs_common_fermi": (
                        f"{record['delta_g_w_relative_change_vs_common_fermi']:.16e}"
                        if record.get("delta_g_w_relative_change_vs_common_fermi") is not None
                        else ""
                    ),
                    "step_size_relative_difference_dw": f"{record['step_size_relative_difference_dw']:.16e}",
                    "stability_verdict": record["stability_verdict"],
                    "sampled_delta_g_w_reference": (
                        f"{sampled['delta_g_w_reference']:.16e}"
                        if sampled["delta_g_w_reference"] is not None
                        else ""
                    ),
                    "delta_g_ratio_direct_over_sampled": (
                        f"{sampled['delta_g_ratio_direct_over_sampled']:.16e}"
                        if sampled["delta_g_ratio_direct_over_sampled"] is not None
                        else ""
                    ),
                    "energy_route": record["energy_route"],
                    "radius_scaling_mode": record["radius_scaling_mode"],
                    "weak_convention_label": record["weak_convention_label"],
                }
            )


def print_mu_point(point: dict[str, Any]) -> None:
    print(f"mu = {point['mu']:.10f}")
    if point["weak_model"] == "common-sphere":
        print(f"Effective RMS radius passed to direct DKB [fm]: {point['effective_rms_fm']:.10f}")
        print(f"Effective common weak radius returned [fm]: {point['effective_radius_fm']:.10f}")
        print(f"Returned R_rel: {point['r_rel']:.16e}")
    elif point["weak_model"] == "pn-spheres":
        print(f"Effective proton RMS radius passed to direct DKB [fm]: {point['effective_proton_rms_fm']:.10f}")
        print(f"Effective neutron RMS radius passed to direct DKB [fm]: {point['effective_neutron_rms_fm']:.10f}")
        print(f"Effective proton sphere radius returned [fm]: {point['effective_proton_radius_fm']:.10f}")
        print(f"Effective neutron sphere radius returned [fm]: {point['effective_neutron_radius_fm']:.10f}")
        print(f"Returned proton R_rel: {point['proton_r_rel']:.16e}")
        print(f"Returned neutron R_rel: {point['neutron_r_rel']:.16e}")
    elif point["weak_model"] == "common-fermi":
        print(f"Effective RMS radius passed to direct DKB [fm]: {point['effective_rms_fm']:.10f}")
        print(f"Effective Fermi diffuseness a [fm]: {point['effective_fermi_a_fm']:.10f}")
        print(f"Effective Fermi half-density radius c [fm]: {point['effective_fermi_c_fm']:.10f}")
        print(f"Effective Coulomb sphere radius returned [fm]: {point['effective_sphere_radius_fm_for_coulomb']:.10f}")
        print(f"Returned Fermi c_rel: {point['fermi_c_rel']:.16e}")
        print(f"Returned Fermi a_rel: {point['fermi_a_rel']:.16e}")
        print(f"Density normalization check: {point['normalization_check']:.16e}")
    elif point["weak_model"] == "pn-fermi":
        print(f"Effective proton RMS radius passed to direct DKB [fm]: {point['effective_proton_rms_fm']:.10f}")
        print(f"Effective neutron RMS radius passed to direct DKB [fm]: {point['effective_neutron_rms_fm']:.10f}")
        print(f"Effective proton Fermi diffuseness a [fm]: {point['effective_proton_a_fm']:.10f}")
        print(f"Effective neutron Fermi diffuseness a [fm]: {point['effective_neutron_a_fm']:.10f}")
        print(f"Effective proton Fermi half-density c [fm]: {point['effective_proton_c_fm']:.10f}")
        print(f"Effective neutron Fermi half-density c [fm]: {point['effective_neutron_c_fm']:.10f}")
        print(f"Coulomb sphere radius from proton rms [fm]: {point['effective_proton_sphere_radius_fm_for_coulomb']:.10f}")
        print(f"Proton density normalization check: {point['proton_normalization_check']:.16e}")
        print(f"Neutron density normalization check: {point['neutron_normalization_check']:.16e}")
        if point["coulomb_model"] == "fermi":
            print(f"Fermi Coulomb V(cp): {point['fermi_coulomb_v_cp']:.16e}")
            print(f"Sphere Coulomb V(cp): {point['sphere_coulomb_v_cp']:.16e}")
    print(f"Raw direct-DKB Delta E_W before explicit mu^2 scaling: {point['raw_direct_dkb_delta_epsilon']:.16e}")
    print(f"Explicit weak normalization factor applied for fixed physical G_F [mu^2]: {point['weak_scale_factor']:.16e}")
    print(f"Delta_epsilon_W(mu) after mu^2 scaling: {point['delta_epsilon']:.16e}")
    print(f"mu * Delta_epsilon_W(mu): {point['mu_times_delta_epsilon']:.16e}")


def print_record(record: dict[str, Any]) -> None:
    print("")
    print(f"{record['isotope_label']} direct-DKB mass-derivative prototype")
    print(f"Energy route: {record['energy_route']}")
    print(f"Requested direct-DKB npoints: {record['npoints_requested']}")
    print(f"Z = {record['Z']}, N = {record['N']}, physical rms_fm = {record['physical_rms_fm']:.10f}")
    print(f"Weak model: {record['weak_model']}")
    print(f"Coulomb model: {record['coulomb_model']}")
    if record["weak_model"] == "common-sphere":
        print(f"Physical common weak radius [fm]: {record['physical_common_sphere_radius_fm']:.10f}")
    elif record["weak_model"] == "pn-spheres":
        print(f"Physical proton rms radius [fm]: {record['physical_proton_rms_fm']:.10f}")
        print(f"Physical neutron rms radius [fm]: {record['physical_neutron_rms_fm']:.10f}")
        print(f"Physical neutron skin [fm]: {record['neutron_skin_fm']:.10f}")
    elif record["weak_model"] == "common-fermi":
        print(f"Physical Fermi diffuseness a [fm]: {record['physical_fermi_a_fm']:.10f}")
        print(f"Physical Fermi half-density radius c [fm]: {record['physical_fermi_c_fm']:.10f}")
        print(f"Coulomb model: frozen hard sphere with radius [fm] {record['physical_coulomb_sphere_radius_fm']:.10f}")
    elif record["weak_model"] == "pn-fermi":
        print(f"Physical proton rms radius [fm]: {record['physical_proton_rms_fm']:.10f}")
        print(f"Physical neutron rms radius [fm]: {record['physical_neutron_rms_fm']:.10f}")
        print(f"Physical neutron skin [fm]: {record['neutron_skin_fm']:.10f}")
        print(f"Physical proton Fermi diffuseness a [fm]: {record['physical_proton_a_fm']:.10f}")
        print(f"Physical neutron Fermi diffuseness a [fm]: {record['physical_neutron_a_fm']:.10f}")
        if record["coulomb_model"] == "fermi":
            print("Coulomb model: proton Fermi charge distribution")
        else:
            print(f"Coulomb model: frozen hard sphere with proton-radius [fm] {record['physical_coulomb_sphere_radius_fm']:.10f}")
    for step in sorted(record["steps"], key=lambda item: item["h"]):
        print("")
        print(f"=== Step size h = {step['h']:.6e} ===")
        print_mu_point(step["mu_minus"])
        print("")
        print_mu_point(step["mu_zero"])
        print("")
        print_mu_point(step["mu_plus"])
        print("")
        print(f"Finite-difference derivative D_W = {step['derivative_dw']:.16e}")
        print(f"delta_g_W(1s) = (4/3) * D_W = {step['delta_g_w']:.16e}")
    print("")
    print(
        f"Step-size relative difference in D_W: "
        f"{record['step_size_relative_difference_dw']:.16e}"
    )
    print(f"Stability verdict: {record['stability_verdict']}")
    sampled = record["sampled_route_reference"]
    if sampled["delta_g_w_reference"] is not None:
        print(
            "Old sampled-route delta_g_W reference: "
            f"{sampled['delta_g_w_reference']:.16e}"
        )
        print(
            "Direct/sample delta_g ratio: "
            f"{sampled['delta_g_ratio_direct_over_sampled']:.16e}"
        )


def add_skin_relative_changes(records: list[dict[str, Any]]) -> None:
    skin_zero_by_ion_model: dict[tuple[str, str], float] = {}
    for record in records:
        if abs(record.get("neutron_skin_fm", 1.0)) < 1.0e-14:
            skin_zero_by_ion_model[(record["isotope_label"], record["coulomb_model"])] = record["delta_g_w_reference"]

    for record in records:
        baseline = skin_zero_by_ion_model.get((record["isotope_label"], record["coulomb_model"]))
        if baseline not in (None, 0.0):
            record["delta_g_w_relative_change_vs_skin0"] = (
                record["delta_g_w_reference"] - baseline
            ) / baseline
            record["relative_magnitude_change_vs_skin0"] = (
                abs(record["delta_g_w_reference"]) - abs(baseline)
            ) / abs(baseline)


def add_coulomb_comparison_changes(records: list[dict[str, Any]]) -> None:
    sphere_by_ion_skin: dict[tuple[str, float], float] = {}
    for record in records:
        if record["coulomb_model"] == "sphere":
            sphere_by_ion_skin[(record["isotope_label"], record.get("neutron_skin_fm", 0.0))] = record[
                "delta_g_w_reference"
            ]

    for record in records:
        sphere_value = sphere_by_ion_skin.get((record["isotope_label"], record.get("neutron_skin_fm", 0.0)))
        if sphere_value in (None, 0.0):
            continue
        record["delta_g_w_sphere_coulomb_reference"] = sphere_value
        record["delta_g_w_absolute_difference_vs_sphere_coulomb"] = abs(record["delta_g_w_reference"] - sphere_value)
        record["delta_g_w_relative_difference_vs_sphere_coulomb"] = (
            record["delta_g_w_reference"] - sphere_value
        ) / sphere_value


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.expanduser()
    if args.coulomb_model == "fermi" and args.weak_model != "pn-fermi":
        raise DirectDkbMassDerivativeError("--coulomb-model fermi is currently supported only with --weak-model pn-fermi.")
    if args.coulomb_model == "sphere" and args.weak_model != "pn-fermi":
        # Existing non-p/n-Fermi routes are hard-sphere Coulomb routes by construction.
        pass
    if args.weak_model == "common-sphere":
        json_output = output_dir / "direct_dkb_mass_derivative.json"
        csv_output = output_dir / "direct_dkb_mass_derivative.csv"
    elif args.weak_model == "pn-spheres":
        json_output = output_dir / "pn_sphere_validation.json"
        csv_output = output_dir / "pn_sphere_validation.csv"
    elif args.weak_model == "pn-fermi":
        if args.coulomb_model == "fermi":
            json_output = output_dir / "fermi_coulomb_comparison.json"
            csv_output = output_dir / "fermi_coulomb_comparison.csv"
        else:
            json_output = output_dir / "pn_fermi_validation.json"
            csv_output = output_dir / "pn_fermi_validation.csv"
    else:
        json_output = output_dir / "common_fermi_validation.json"
        csv_output = output_dir / "common_fermi_validation.csv"

    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise DirectDkbMassDerivativeError("Imported module does not expose `rad_kernel_f2py`.")
    if not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s"):
        raise DirectDkbMassDerivativeError("Imported f2py module does not expose `direct_dkb_weak_shift_1s`; rebuild rad_f2py.")
    if args.npoints != 1000 and not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_npoints"):
        raise DirectDkbMassDerivativeError("Imported f2py module does not expose `direct_dkb_weak_shift_1s_npoints`; rebuild rad_f2py.")
    if args.weak_model == "pn-spheres" and not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_pn_spheres"):
        raise DirectDkbMassDerivativeError(
            "Imported f2py module does not expose `direct_dkb_weak_shift_1s_pn_spheres`; rebuild rad_f2py."
        )
    if args.weak_model == "common-fermi" and not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_common_fermi"):
        raise DirectDkbMassDerivativeError(
            "Imported f2py module does not expose `direct_dkb_weak_shift_1s_common_fermi`; rebuild rad_f2py."
        )
    if args.weak_model == "pn-fermi" and not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_pn_fermi"):
        raise DirectDkbMassDerivativeError(
            "Imported f2py module does not expose `direct_dkb_weak_shift_1s_pn_fermi`; rebuild rad_f2py."
        )
    if (
        args.weak_model == "pn-fermi"
        and args.coulomb_model == "fermi"
        and not hasattr(module.rad_kernel_f2py, "direct_dkb_weak_shift_1s_pn_fermi_fermi_coulomb")
    ):
        raise DirectDkbMassDerivativeError(
            "Imported f2py module does not expose `direct_dkb_weak_shift_1s_pn_fermi_fermi_coulomb`; rebuild rad_f2py."
        )
    if args.weak_model == "common-sphere":
        if (
            args.proton_rms_fm is not None
            or args.neutron_rms_fm is not None
            or args.neutron_skin_fm is not None
            or args.rms_fm is not None
            or args.fermi_a_fm is not None
            or args.proton_a_fm is not None
            or args.neutron_a_fm is not None
        ):
            raise DirectDkbMassDerivativeError("Radius/Fermi override options require a non-common-sphere weak model.")
        ion_inputs = [dict(ion) for ion in ION_SET]
    elif args.weak_model == "pn-spheres":
        if args.rms_fm is not None or args.fermi_a_fm is not None or args.proton_a_fm is not None or args.neutron_a_fm is not None:
            raise DirectDkbMassDerivativeError("Fermi options require --weak-model common-fermi or pn-fermi.")
        ion_inputs = []
        if args.neutron_rms_fm is None and args.neutron_skin_fm is None:
            for ion in ION_SET:
                for skin_fm in DEFAULT_NEUTRON_SKIN_FM_VALUES:
                    proton_value, neutron_value = resolve_pn_radii(
                        ion,
                        args.proton_rms_fm,
                        None,
                        skin_fm,
                    )
                    ion_inputs.append(
                        {
                            **ion,
                            "proton_rms_fm": proton_value,
                            "neutron_rms_fm": neutron_value,
                        }
                    )
        else:
            for ion in ION_SET:
                proton_value, neutron_value = resolve_pn_radii(
                    ion,
                    args.proton_rms_fm,
                    args.neutron_rms_fm,
                    args.neutron_skin_fm,
                )
                ion_inputs.append(
                    {
                        **ion,
                        "proton_rms_fm": proton_value,
                        "neutron_rms_fm": neutron_value,
                    }
                )
    elif args.weak_model == "common-fermi":
        if (
            args.proton_rms_fm is not None
            or args.neutron_rms_fm is not None
            or args.neutron_skin_fm is not None
            or args.proton_a_fm is not None
            or args.neutron_a_fm is not None
        ):
            raise DirectDkbMassDerivativeError("Proton/neutron radius or diffuseness options require --weak-model pn-spheres or pn-fermi.")
        fermi_a_fm = DEFAULT_FERMI_A_FM if args.fermi_a_fm is None else args.fermi_a_fm
        if fermi_a_fm <= 0.0:
            raise DirectDkbMassDerivativeError("Fermi diffuseness must be positive.")
        ion_inputs = []
        for ion in ION_SET:
            rms_value = ion["rms_fm"] if args.rms_fm is None else args.rms_fm
            if rms_value <= 0.0:
                raise DirectDkbMassDerivativeError("Fermi rms radius must be positive.")
            fermi_c_from_rms(rms_value, fermi_a_fm)
            ion_inputs.append({**ion, "rms_fm": rms_value, "fermi_a_fm": fermi_a_fm})
    else:
        if args.rms_fm is not None or args.fermi_a_fm is not None:
            raise DirectDkbMassDerivativeError("Common Fermi options require --weak-model common-fermi.")
        proton_a_fm = DEFAULT_FERMI_A_FM if args.proton_a_fm is None else args.proton_a_fm
        neutron_a_fm = proton_a_fm if args.neutron_a_fm is None else args.neutron_a_fm
        if proton_a_fm <= 0.0 or neutron_a_fm <= 0.0:
            raise DirectDkbMassDerivativeError("Proton and neutron Fermi diffuseness values must be positive.")
        ion_inputs = []
        if args.neutron_rms_fm is None and args.neutron_skin_fm is None:
            for ion in ION_SET:
                for skin_fm in DEFAULT_NEUTRON_SKIN_FM_VALUES:
                    proton_value, neutron_value = resolve_pn_radii(ion, args.proton_rms_fm, None, skin_fm)
                    fermi_c_from_rms(proton_value, proton_a_fm)
                    fermi_c_from_rms(neutron_value, neutron_a_fm)
                    ion_inputs.append(
                        {
                            **ion,
                            "proton_rms_fm": proton_value,
                            "neutron_rms_fm": neutron_value,
                            "proton_a_fm": proton_a_fm,
                            "neutron_a_fm": neutron_a_fm,
                        }
                    )
        else:
            for ion in ION_SET:
                proton_value, neutron_value = resolve_pn_radii(
                    ion,
                    args.proton_rms_fm,
                    args.neutron_rms_fm,
                    args.neutron_skin_fm,
                )
                fermi_c_from_rms(proton_value, proton_a_fm)
                fermi_c_from_rms(neutron_value, neutron_a_fm)
                ion_inputs.append(
                    {
                        **ion,
                        "proton_rms_fm": proton_value,
                        "neutron_rms_fm": neutron_value,
                        "proton_a_fm": proton_a_fm,
                        "neutron_a_fm": neutron_a_fm,
                    }
                )

    print("Direct-DKB weak mass-derivative prototype")
    print("Physics route: virial / Karshenboim-Lee-Milstein mass derivative")
    print("Energy route: direct-dkb")
    print("Direct magnetic-response route status: not used here")
    print(f"Import route: {import_route}")
    print(f"Imported module: {module.__name__}")
    print(f"Module file: {getattr(module, '__file__', 'n/a')}")
    print(f"Weak model: {args.weak_model}")
    print(f"Coulomb model: {args.coulomb_model}")
    print(f"Fortran routine used for weak energies: {fortran_routine_name(args.weak_model, args.coulomb_model)}")
    print(f"Requested direct-DKB npoints: {args.npoints}")
    print(f"Weak convention label: {weak_convention_label(args.weak_model)}")

    sampled_reference = load_sampled_reference()
    if args.weak_model == "pn-fermi" and args.coulomb_model == "fermi":
        records = []
        for ion in ion_inputs:
            records.append(summarize_ion(module, ion, sampled_reference, args.npoints, args.weak_model, "sphere"))
            records.append(summarize_ion(module, ion, sampled_reference, args.npoints, args.weak_model, "fermi"))
    else:
        records = [
            summarize_ion(module, ion, sampled_reference, args.npoints, args.weak_model, args.coulomb_model)
            for ion in ion_inputs
        ]
    if args.weak_model in ("pn-spheres", "pn-fermi"):
        add_skin_relative_changes(records)
    if args.weak_model == "pn-fermi" and args.coulomb_model == "fermi":
        add_coulomb_comparison_changes(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(records, json_output, import_route, module, args.weak_model, args.coulomb_model)
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
    except (DirectDkbMassDerivativeError, FirstOrderError, MassDerivativeError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
