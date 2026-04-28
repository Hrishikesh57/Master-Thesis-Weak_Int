#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import math
import sys
from pathlib import Path


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parent
PYTHON_API_DIR = WORKSPACE_F2PY_ROOT / "python_api"
RAD_F2PY_DIR = WORKSPACE_F2PY_ROOT / "programs" / "rad_f2py"

# Constants mirror workspace_f2py/codes/bspqed/src/bspqed_physconsts.f90.
HBARC_MEV_FM = 197.3269788
ELECTRON_REST_ENERGY_MEV = 0.5109989461
ELECTRON_REST_ENERGY_GEV = 0.0005109989461
GFERMI_GEV_MINUS2 = 1.1663787e-5
SIN2_THETA_W = 0.2312


class WeakSourceDiagnosticError(RuntimeError):
    """User-facing weak-source diagnostic error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Print the current common-radius weak-source convention used by the "
            "sampled-rad f2py prototype. This is a convention diagnostic only."
        )
    )
    parser.add_argument("--Z", type=int, required=True, dest="znucl", help="Atomic number Z.")
    parser.add_argument("--N", type=int, required=True, dest="nneutrons", help="Neutron number N.")
    radius_group = parser.add_mutually_exclusive_group(required=True)
    radius_group.add_argument(
        "--rms-fm",
        type=float,
        help="Nuclear RMS radius in fm; converted to the current common sphere radius sqrt(5/3) * rms.",
    )
    radius_group.add_argument(
        "--radius-fm",
        type=float,
        help="Common sphere radius R in fm, passed directly to the current weak source.",
    )
    return parser.parse_args()


def resolve_radius(args: argparse.Namespace) -> tuple[float, str]:
    if args.radius_fm is not None:
        if args.radius_fm <= 0.0:
            raise WeakSourceDiagnosticError("--radius-fm must be positive.")
        return float(args.radius_fm), "radius_fm"

    if args.rms_fm is None or args.rms_fm <= 0.0:
        raise WeakSourceDiagnosticError("--rms-fm must be positive.")
    return math.sqrt(5.0 / 3.0) * float(args.rms_fm), "rms_fm"


def load_rad_kernel_suite():
    diagnostics: list[str] = []
    if PYTHON_API_DIR.is_dir():
        python_api_dir = str(PYTHON_API_DIR)
        if python_api_dir not in sys.path:
            sys.path.insert(0, python_api_dir)
        try:
            from workspace_f2py_api import load_rad_kernel_suite_module

            return load_rad_kernel_suite_module(), "workspace_f2py_api.load_rad_kernel_suite_module()"
        except Exception as exc:  # pragma: no cover - import environment dependent
            diagnostics.append(f"workspace_f2py_api route failed: {type(exc).__name__}: {exc}")

    if RAD_F2PY_DIR.is_dir():
        rad_f2py_dir = str(RAD_F2PY_DIR)
        if rad_f2py_dir not in sys.path:
            sys.path.insert(0, rad_f2py_dir)
        try:
            return importlib.import_module("rad_kernel_suite"), "direct import from programs/rad_f2py"
        except Exception as exc:  # pragma: no cover - import environment dependent
            diagnostics.append(f"direct import route failed: {type(exc).__name__}: {exc}")

    return None, "; ".join(diagnostics) if diagnostics else "compiled f2py module not found"


def formula_values(znucl: int, nneutrons: int, radius_fm: float) -> dict[str, float]:
    electron_reduced_compton = HBARC_MEV_FM / ELECTRON_REST_ENERGY_MEV
    radius_rel = radius_fm / electron_reduced_compton
    gfermi_relativistic = GFERMI_GEV_MINUS2 * ELECTRON_REST_ENERGY_GEV * ELECTRON_REST_ENERGY_GEV
    cv_e = 1.0 - 4.0 * SIN2_THETA_W
    cv_p = 1.0 - 4.0 * SIN2_THETA_W
    source_factor = cv_p * float(znucl) - float(nneutrons)
    rho0 = (
        (3.0 / (4.0 * math.pi))
        * (1.0 / radius_rel) ** 3
        * source_factor
        * (gfermi_relativistic / (2.0 * math.sqrt(2.0)))
    )
    return {
        "electron_reduced_compton": electron_reduced_compton,
        "radius_rel": radius_rel,
        "gfermi_relativistic": gfermi_relativistic,
        "cv_e": cv_e,
        "cv_p": cv_p,
        "source_factor": source_factor,
        "rho0": rho0,
        "pc_half_radius": cv_e * rho0,
        "pv_half_radius": rho0,
    }


def main() -> int:
    args = parse_args()
    if args.znucl < 1:
        raise WeakSourceDiagnosticError("--Z must be positive.")
    if args.nneutrons < 0:
        raise WeakSourceDiagnosticError("--N must be nonnegative.")

    radius_fm, radius_mode = resolve_radius(args)
    values = formula_values(args.znucl, args.nneutrons, radius_fm)
    module, import_route = load_rad_kernel_suite()

    print("Weak-source convention diagnostic")
    print("Diagnostic status: convention check only; no wavefunction integral is evaluated")
    print(f"Z: {args.znucl}")
    print(f"N: {args.nneutrons}")
    if radius_mode == "rms_fm":
        print(f"Input RMS radius [fm]: {args.rms_fm:.10f}")
        print(f"Common sphere radius passed to weak kernel [fm]: {radius_fm:.10f}")
    else:
        print(f"Input common sphere radius [fm]: {radius_fm:.10f}")
    print(f"electron_reduced_compton [fm]: {values['electron_reduced_compton']:.16e}")
    print(f"radius_rel = radius_fm / electron_reduced_compton: {values['radius_rel']:.16e}")
    print(f"Gfermi_relativistic = G_F * m_e^2: {values['gfermi_relativistic']:.16e}")
    print(f"Cv_e implemented value: {values['cv_e']:.16e}")
    print(f"Cv_p implemented value: {values['cv_p']:.16e}")
    print(f"source factor Cv_p * Z - N: {values['source_factor']:.16e}")
    print(f"rho0: {values['rho0']:.16e}")
    print(f"V_W^pc(0.5 R) = Cv_e * rho0: {values['pc_half_radius']:.16e}")
    print(f"V_W^pv(0.5 R) = rho0: {values['pv_half_radius']:.16e}")

    if module is None or not hasattr(module, "rad_kernel_f2py"):
        print(f"f2py weak_potential_summary route: unavailable ({import_route})")
        print("Printed values are the Python mirror of the inspected Fortran constants and formula.")
        return 0

    summary = module.rad_kernel_f2py.weak_potential_summary
    r_rel_f, rho0_f, qweak_f, pc_f, pv_f, ierr = summary(args.znucl, args.nneutrons, radius_fm)
    if ierr != 0:
        raise WeakSourceDiagnosticError(f"weak_potential_summary returned ierr={int(ierr)}")
    print(f"f2py weak_potential_summary route: {import_route}")
    print(f"f2py radius_rel: {float(r_rel_f):.16e}")
    print(f"f2py source factor: {float(qweak_f):.16e}")
    print(f"f2py rho0: {float(rho0_f):.16e}")
    print(f"f2py V_W^pc(0.5 R): {float(pc_f):.16e}")
    print(f"f2py V_W^pv(0.5 R): {float(pv_f):.16e}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except WeakSourceDiagnosticError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
