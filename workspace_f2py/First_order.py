#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import math
import sys
from pathlib import Path

import numpy as np


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parent
PYTHON_API_DIR = WORKSPACE_F2PY_ROOT / "python_api"
RAD_F2PY_DIR = WORKSPACE_F2PY_ROOT / "programs" / "rad_f2py"
RAD_F2PY_BUILD_SCRIPT = RAD_F2PY_DIR / "build_f2py.sh"
FORTRAN_ROUTINE_NAME = "rad_kernel_suite.rad_kernel_f2py.first_order_weak_shift"
F2PY_SIGNATURE = "e_pot, r_rel, ierr = first_order_weak_shift(znucl, nneutrons, r_fm, r, g, f)"
NORMALIZATION_TOLERANCE = 5.0e-2

FORTRAN_IERR_MESSAGES = {
    1: "sampled wavefunction must contain at least two radial points",
    2: "sampled r, G, and F arrays must have the same length",
    3: "sampled radial grid must be strictly increasing",
    11: "znucl must be positive",
    12: "nneutrons must be nonnegative",
    13: "r_fm must be positive",
}


class FirstOrderError(RuntimeError):
    """User-facing wrapper error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute the sampled-rad first-order parity-conserving weak energy "
            "shift by calling the existing f2py rad kernel."
        )
    )
    parser.add_argument("--Z", type=int, required=True, dest="znucl", help="Atomic number Z.")
    parser.add_argument("--N", type=int, required=True, dest="nneutrons", help="Neutron number N.")
    radius_group = parser.add_mutually_exclusive_group(required=True)
    radius_group.add_argument(
        "--rms-fm",
        type=float,
        help=(
            "Nuclear RMS radius in fm. For the current common-radius weak kernel this is "
            "converted to the sphere radius R = sqrt(5/3) * rms before calling Fortran."
        ),
    )
    radius_group.add_argument(
        "--radius-fm",
        type=float,
        help="Common sphere radius R in fm passed directly to the current weak kernel.",
    )
    parser.add_argument(
        "--wavefunction",
        type=Path,
        required=True,
        help=(
            "Path to sampled radial data with columns r, G(r), F(r) in relativistic units. "
            "These are the radial Dirac components used by the Fortran rad path, i.e. "
            "G = r g and F = r f, not the reduced little-g/little-f functions."
        ),
    )
    parser.add_argument(
        "--allow-unnormalized",
        action="store_true",
        help=(
            "Allow a badly normalized sampled wavefunction. Use this only for legacy "
            "path reproduction, not for production weak shifts."
        ),
    )
    return parser.parse_args()


def resolve_radius_input(args: argparse.Namespace) -> tuple[float, str]:
    if args.radius_fm is not None:
        if args.radius_fm <= 0.0:
            raise FirstOrderError("--radius-fm must be positive.")
        return float(args.radius_fm), "radius_fm"

    if args.rms_fm is None or args.rms_fm <= 0.0:
        raise FirstOrderError("--rms-fm must be positive.")

    # The current sampled weak kernel takes a single common sphere radius R.
    return math.sqrt(5.0 / 3.0) * float(args.rms_fm), "rms_fm"


def read_sampled_wavefunction(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    wavefunction_path = path.expanduser()
    if not wavefunction_path.is_file():
        raise FirstOrderError(f"Wavefunction file was not found: {wavefunction_path}")

    r_values: list[float] = []
    g_values: list[float] = []
    f_values: list[float] = []

    for raw_line in wavefunction_path.read_text().splitlines():
        stripped = raw_line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 3:
            continue
        try:
            rval = float(parts[0])
            gval = float(parts[1])
            fval = float(parts[2])
        except ValueError:
            continue
        r_values.append(rval)
        g_values.append(gval)
        f_values.append(fval)

    if len(r_values) < 2:
        raise FirstOrderError(
            "No sampled radial data were found in the wavefunction file. "
            "This wrapper needs sampled columns r, G(r), F(r) with G = r g and F = r f; "
            "basis/spectrum files are not enough."
        )

    r_array = np.asarray(r_values, dtype=np.float64)
    g_array = np.asarray(g_values, dtype=np.float64)
    f_array = np.asarray(f_values, dtype=np.float64)

    if not np.all(np.diff(r_array) > 0.0):
        raise FirstOrderError("Wavefunction radial grid must be strictly increasing.")

    return r_array, g_array, f_array


def import_rad_kernel_suite():
    diagnostics: list[str] = []

    if PYTHON_API_DIR.is_dir():
        python_api_dir = str(PYTHON_API_DIR)
        if python_api_dir not in sys.path:
            sys.path.insert(0, python_api_dir)
        try:
            from workspace_f2py_api import load_rad_kernel_suite_module
        except Exception as exc:  # pragma: no cover - import environment dependent
            diagnostics.append(f"Import via workspace_f2py_api failed: {type(exc).__name__}: {exc}")
        else:
            try:
                module = load_rad_kernel_suite_module()
                return module, "workspace_f2py_api.load_rad_kernel_suite_module()"
            except Exception as exc:  # pragma: no cover - import environment dependent
                diagnostics.append(f"Loading rad_kernel_suite via workspace_f2py_api failed: {type(exc).__name__}: {exc}")

    if RAD_F2PY_DIR.is_dir():
        rad_f2py_dir = str(RAD_F2PY_DIR)
        if rad_f2py_dir not in sys.path:
            sys.path.insert(0, rad_f2py_dir)
        try:
            module = importlib.import_module("rad_kernel_suite")
            return module, "direct import from workspace_f2py/programs/rad_f2py"
        except Exception as exc:  # pragma: no cover - import environment dependent
            diagnostics.append(f"Direct import of rad_kernel_suite failed: {type(exc).__name__}: {exc}")

    module_candidates = sorted(RAD_F2PY_DIR.glob("rad_kernel_suite*.so")) + sorted(RAD_F2PY_DIR.glob("rad_kernel_suite*.pyd"))
    message_lines = ["Could not import the f2py module `rad_kernel_suite` needed for the sampled-rad weak shift."]
    if module_candidates:
        message_lines.append("Compiled module candidates found:")
        message_lines.extend(f"  - {candidate}" for candidate in module_candidates)
    else:
        message_lines.append(f"No compiled `rad_kernel_suite` extension was found under {RAD_F2PY_DIR}.")
    message_lines.append(f"Build hint: bash {RAD_F2PY_BUILD_SCRIPT}")
    if diagnostics:
        message_lines.append("Import diagnostics:")
        message_lines.extend(f"  - {item}" for item in diagnostics)
    raise FirstOrderError("\n".join(message_lines))


def normalization_trapz(r_array: np.ndarray, g_array: np.ndarray, f_array: np.ndarray) -> float:
    density = g_array * g_array + f_array * f_array
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(density, x=r_array))
    return float(np.trapz(density, x=r_array))


def normalization_warning_message(normalization: float, wavefunction_path: Path) -> str:
    lines = [
        (
            f"Wavefunction normalization check failed: int(G^2 + F^2) dr = {normalization:.16e}, "
            f"expected about 1 within {NORMALIZATION_TOLERANCE:.1e}."
        ),
        (
            "Physical reason: Delta E_W scales with the sampled density, so a badly normalized "
            "input wavefunction would silently scale the weak-energy shift."
        ),
        (
            "Code reason: the current wrapper accepts pre-sampled r, G, F data and does not "
            "renormalize them automatically, to avoid silently changing the physics."
        ),
    ]
    if wavefunction_path.name == "radial_function_relunits.txt":
        lines.append(
            "Repository note: `data/reference/radial_functions/radial_function_relunits.txt` is a legacy "
            "reference file used by the old `something.f90` energy path and is not unit-normalized."
        )
        lines.append(
            "Use a normalized sampled export from `programs/rad/outputs/wavefunctions/...`, for example "
            "`92_polynomial_table_1s.txt`, or rerun with `--allow-unnormalized` only if you explicitly "
            "want legacy-path reproduction."
        )
    return "\n".join(lines)


def format_fortran_error(ierr: int) -> str:
    detail = FORTRAN_IERR_MESSAGES.get(ierr, "unknown Fortran-side validation error")
    return f"Fortran routine returned ierr={ierr}: {detail}"


def main() -> int:
    args = parse_args()
    radius_fm, radius_mode = resolve_radius_input(args)
    wavefunction_path = args.wavefunction.expanduser()
    r_array, g_array, f_array = read_sampled_wavefunction(wavefunction_path)
    normalization = normalization_trapz(r_array, g_array, f_array)
    if abs(normalization - 1.0) > NORMALIZATION_TOLERANCE and not args.allow_unnormalized:
        raise FirstOrderError(normalization_warning_message(normalization, wavefunction_path))

    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise FirstOrderError("Imported module does not expose `rad_kernel_f2py`.")

    routine = module.rad_kernel_f2py.first_order_weak_shift
    e_pot, r_rel, ierr = routine(args.znucl, args.nneutrons, radius_fm, r_array, g_array, f_array)
    if ierr != 0:
        raise FirstOrderError(format_fortran_error(int(ierr)))

    print("Sampled-rad first-order weak-energy wrapper")
    print(f"Import route: {import_route}")
    print(f"Imported module: {module.__name__}")
    print(f"Module file: {getattr(module, '__file__', 'n/a')}")
    print(f"Fortran routine called: {FORTRAN_ROUTINE_NAME}")
    print(f"Discovered f2py signature: {F2PY_SIGNATURE}")
    print(f"Wavefunction file: {wavefunction_path}")
    print("Wavefunction convention: sampled columns are r, G(r), F(r) with G = r g and F = r f")
    print(f"Z: {args.znucl}")
    print(f"N: {args.nneutrons}")
    if radius_mode == "radius_fm":
        print(f"Input radius in fm (passed directly): {radius_fm:.10f}")
    else:
        print(f"Input RMS radius in fm: {args.rms_fm:.10f}")
        print(f"Derived common sphere radius in fm passed to Fortran [R = sqrt(5/3) * rms]: {radius_fm:.10f}")
    print(f"Number of radial grid points: {r_array.size}")
    print(f"min(r): {r_array[0]:.10e}")
    print(f"max(r): {r_array[-1]:.10e}")
    print(f"Sampled normalization check int(G^2 + F^2) dr [trapz]: {normalization:.16e}")
    if abs(normalization - 1.0) > NORMALIZATION_TOLERANCE:
        print("WARNING: proceeding with an unnormalized sampled wavefunction because --allow-unnormalized was set.")
    print(f"Weak-kernel radius in relativistic units returned by Fortran: {r_rel:.16e}")
    print(f"Delta E_W from Fortran sampled-rad first-order weak shift: {e_pot:.16e}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FirstOrderError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
