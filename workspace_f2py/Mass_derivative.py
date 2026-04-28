#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from First_order import (
    FORTRAN_ROUTINE_NAME,
    NORMALIZATION_TOLERANCE,
    FirstOrderError,
    format_fortran_error,
    import_rad_kernel_suite,
    normalization_trapz,
    read_sampled_wavefunction,
)


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parent
RAD_DIR = WORKSPACE_F2PY_ROOT / "programs" / "rad"
PRACTISE_EXECUTABLE = RAD_DIR / "bin" / "practise"
STANDARD_WAVEFUNCTION_DIR = RAD_DIR / "outputs" / "wavefunctions" / "single_electron"
ARCHIVE_WAVEFUNCTION_DIR = RAD_DIR / "outputs" / "wavefunctions" / "mass_derivative"
STATE_LABEL = "1s"


class MassDerivativeError(RuntimeError):
    """User-facing mass-derivative wrapper error."""


@dataclass
class MuPointResult:
    mu: float
    effective_rms_fm: float
    effective_radius_fm: float
    wavefunction_mode: str
    wavefunction_path: Path
    normalization: float
    used_unnormalized_diagnostic: bool
    npoints: int
    r_min: float
    r_max: float
    raw_delta_epsilon: float
    weak_scale_factor: float
    delta_epsilon: float
    mu_times_delta_epsilon: float
    r_rel: float


@dataclass
class StepResult:
    h: float
    mu_minus: MuPointResult
    mu_zero: MuPointResult
    mu_plus: MuPointResult
    derivative_dw: float
    delta_g_w: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prototype mass-derivative wrapper for the sampled-rad weak 1s shift. "
            "This uses the Karshenboim-Lee-Milstein / virial route, not g_uehling."
        )
    )
    parser.add_argument("--Z", type=int, required=True, dest="znucl", help="Atomic number Z.")
    parser.add_argument("--N", type=int, required=True, dest="nneutrons", help="Neutron number N.")
    radius_group = parser.add_mutually_exclusive_group(required=True)
    radius_group.add_argument(
        "--rms-fm",
        type=float,
        help=(
            "Physical nuclear RMS radius in fm. The single-electron solver uses this directly, "
            "while the current weak kernel uses the derived common sphere radius R = sqrt(5/3) * rms."
        ),
    )
    radius_group.add_argument(
        "--radius-fm",
        type=float,
        help=(
            "Physical common sphere radius in fm for the current weak kernel. "
            "The single-electron solver then uses the derived RMS radius sqrt(3/5) * R."
        ),
    )
    parser.add_argument(
        "--h",
        type=float,
        action="append",
        dest="h_values",
        help="Finite-difference step size h. Repeat the option to test multiple step sizes.",
    )
    parser.add_argument(
        "--wavefunction-mode",
        choices=("regenerate", "reuse"),
        default="regenerate",
        help=(
            "Use `regenerate` for the main physics path. "
            "Use `reuse` only as an incomplete negative control."
        ),
    )
    parser.add_argument(
        "--baseline-wavefunction",
        type=Path,
        help=(
            "Optional sampled 1s wavefunction to reuse for all mu values in `reuse` mode. "
            "If omitted, the script generates a mu=1 baseline and reuses it."
        ),
    )
    parser.add_argument(
        "--radius-scaling-mode",
        choices=("fixed-physical", "fixed-code"),
        default="fixed-physical",
        help=(
            "Use `fixed-physical` for the intended derivative at fixed physical radius, "
            "implemented by passing mu-scaled fm inputs to the existing code. "
            "`fixed-code` is an incorrect control with code radius held fixed."
        ),
    )
    parser.add_argument(
        "--allow-unnormalized-diagnostic",
        action="store_true",
        help=(
            "Allow an unnormalized sampled wavefunction only for explicit provenance diagnostics, "
            "for example when reusing the legacy radial_function_relunits.txt table. "
            "Do not use this for the physical regenerated-wavefunction derivative."
        ),
    )
    return parser.parse_args()


def resolve_physical_radii(args: argparse.Namespace) -> tuple[float, float, str]:
    if args.rms_fm is not None:
        if args.rms_fm <= 0.0:
            raise MassDerivativeError("--rms-fm must be positive.")
        physical_rms_fm = float(args.rms_fm)
        physical_radius_fm = math.sqrt(5.0 / 3.0) * physical_rms_fm
        return physical_rms_fm, physical_radius_fm, "rms_fm"

    if args.radius_fm is None or args.radius_fm <= 0.0:
        raise MassDerivativeError("--radius-fm must be positive.")

    physical_radius_fm = float(args.radius_fm)
    physical_rms_fm = math.sqrt(3.0 / 5.0) * physical_radius_fm
    return physical_rms_fm, physical_radius_fm, "radius_fm"


def validate_h_values(h_values: list[float] | None) -> list[float]:
    if not h_values:
        return [1.0e-3]

    cleaned: list[float] = []
    for h_value in h_values:
        if h_value <= 0.0:
            raise MassDerivativeError("Each --h value must be positive.")
        if h_value >= 1.0:
            raise MassDerivativeError("Each --h value must be smaller than 1.")
        cleaned.append(float(h_value))
    return cleaned


def standard_wavefunction_path(znucl: int) -> Path:
    return STANDARD_WAVEFUNCTION_DIR / f"{znucl}_polynomial_table_{STATE_LABEL}.txt"


def ensure_runtime_dependencies() -> None:
    if not PRACTISE_EXECUTABLE.is_file():
        raise MassDerivativeError(
            f"Missing wavefunction generator: {PRACTISE_EXECUTABLE}\n"
            f"Build hint: make -C {RAD_DIR} practise"
        )


@contextmanager
def preserve_standard_wavefunction(path: Path):
    existed = path.exists()
    backup = path.read_bytes() if existed else None
    try:
        yield
    finally:
        if existed and backup is not None:
            path.write_bytes(backup)
        elif path.exists():
            path.unlink()


def run_practise_for_1s(znucl: int, rms_fm: float) -> Path:
    input_payload = f"1\n1\n{znucl}\n{STATE_LABEL}\n{rms_fm:.16f}\n"
    completed = subprocess.run(
        [str(PRACTISE_EXECUTABLE)],
        cwd=RAD_DIR,
        input=input_payload,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        stdout_tail = "\n".join(completed.stdout.strip().splitlines()[-20:])
        stderr_tail = "\n".join(completed.stderr.strip().splitlines()[-20:])
        message_lines = [
            "Wavefunction regeneration with `practise` failed.",
            f"Command: {PRACTISE_EXECUTABLE}",
            f"Return code: {completed.returncode}",
        ]
        if stdout_tail:
            message_lines.append("Last stdout lines:")
            message_lines.append(stdout_tail)
        if stderr_tail:
            message_lines.append("Last stderr lines:")
            message_lines.append(stderr_tail)
        raise MassDerivativeError("\n".join(message_lines))

    wavefunction_path = standard_wavefunction_path(znucl)
    if not wavefunction_path.is_file():
        raise MassDerivativeError(
            "Wavefunction regeneration finished without producing the expected sampled file:\n"
            f"  {wavefunction_path}"
        )
    return wavefunction_path


def mu_tag(mu_value: float) -> str:
    return f"{mu_value:.8f}".replace("-", "m").replace(".", "p")


def archive_wavefunction(source_path: Path, znucl: int, mu_value: float, scaling_mode: str) -> Path:
    ARCHIVE_WAVEFUNCTION_DIR.mkdir(parents=True, exist_ok=True)
    archive_name = f"{znucl}_polynomial_table_{STATE_LABEL}_mu_{mu_tag(mu_value)}_{scaling_mode}.txt"
    archive_path = ARCHIVE_WAVEFUNCTION_DIR / archive_name
    shutil.copy2(source_path, archive_path)
    return archive_path


def effective_radius_inputs(
    mu_value: float,
    physical_rms_fm: float,
    physical_radius_fm: float,
    scaling_mode: str,
) -> tuple[float, float]:
    if scaling_mode == "fixed-physical":
        return mu_value * physical_rms_fm, mu_value * physical_radius_fm
    if scaling_mode == "fixed-code":
        return physical_rms_fm, physical_radius_fm
    raise MassDerivativeError(f"Unsupported radius scaling mode: {scaling_mode}")


def compute_delta_epsilon(
    module,
    znucl: int,
    nneutrons: int,
    radius_fm: float,
    wavefunction_path: Path,
    allow_unnormalized_diagnostic: bool,
) -> tuple[float, float, float, bool, int, float, float]:
    r_array, g_array, f_array = read_sampled_wavefunction(wavefunction_path)
    normalization = normalization_trapz(r_array, g_array, f_array)
    used_unnormalized_diagnostic = abs(normalization - 1.0) > NORMALIZATION_TOLERANCE
    if used_unnormalized_diagnostic and not allow_unnormalized_diagnostic:
        raise MassDerivativeError(
            "Wavefunction normalization check failed for the mass-derivative run: "
            f"int(G^2 + F^2) dr = {normalization:.16e}. "
            "This prototype requires normalized sampled 1s wavefunctions."
        )

    routine = module.rad_kernel_f2py.first_order_weak_shift
    raw_delta_epsilon, r_rel, ierr = routine(znucl, nneutrons, radius_fm, r_array, g_array, f_array)
    if ierr != 0:
        raise MassDerivativeError(format_fortran_error(int(ierr)))
    return (
        float(raw_delta_epsilon),
        float(r_rel),
        normalization,
        used_unnormalized_diagnostic,
        int(r_array.size),
        float(r_array[0]),
        float(r_array[-1]),
    )


def prepare_reuse_wavefunction(
    args: argparse.Namespace,
    physical_rms_fm: float,
    physical_radius_fm: float,
    scaling_mode: str,
    standard_path: Path,
) -> Path:
    if args.baseline_wavefunction is not None:
        baseline_path = args.baseline_wavefunction.expanduser()
        if not baseline_path.is_file():
            raise MassDerivativeError(f"Baseline wavefunction was not found: {baseline_path}")
        return baseline_path

    baseline_rms_fm, _ = effective_radius_inputs(1.0, physical_rms_fm, physical_radius_fm, scaling_mode)
    generated_path = run_practise_for_1s(args.znucl, baseline_rms_fm)
    archive_path = archive_wavefunction(generated_path, args.znucl, 1.0, scaling_mode)
    if standard_path.exists():
        shutil.copy2(archive_path, standard_path)
    return archive_path


def evaluate_mu_point(
    module,
    args: argparse.Namespace,
    mu_value: float,
    physical_rms_fm: float,
    physical_radius_fm: float,
    scaling_mode: str,
    reuse_wavefunction_path: Path | None,
) -> MuPointResult:
    effective_rms_fm, effective_radius_fm = effective_radius_inputs(
        mu_value, physical_rms_fm, physical_radius_fm, scaling_mode
    )

    if reuse_wavefunction_path is None:
        generated_path = run_practise_for_1s(args.znucl, effective_rms_fm)
        wavefunction_path = archive_wavefunction(generated_path, args.znucl, mu_value, scaling_mode)
        wavefunction_mode = "regenerated"
    else:
        wavefunction_path = reuse_wavefunction_path
        wavefunction_mode = "reused"

    raw_delta_epsilon, r_rel, normalization, used_unnormalized_diagnostic, npoints, r_min, r_max = compute_delta_epsilon(
        module,
        args.znucl,
        args.nneutrons,
        effective_radius_fm,
        wavefunction_path,
        args.allow_unnormalized_diagnostic,
    )
    weak_scale_factor = mu_value * mu_value
    delta_epsilon = weak_scale_factor * raw_delta_epsilon
    mu_times_delta_epsilon = mu_value * delta_epsilon

    return MuPointResult(
        mu=mu_value,
        effective_rms_fm=effective_rms_fm,
        effective_radius_fm=effective_radius_fm,
        wavefunction_mode=wavefunction_mode,
        wavefunction_path=wavefunction_path,
        normalization=normalization,
        used_unnormalized_diagnostic=used_unnormalized_diagnostic,
        npoints=npoints,
        r_min=r_min,
        r_max=r_max,
        raw_delta_epsilon=raw_delta_epsilon,
        weak_scale_factor=weak_scale_factor,
        delta_epsilon=delta_epsilon,
        mu_times_delta_epsilon=mu_times_delta_epsilon,
        r_rel=r_rel,
    )


def evaluate_step_size(
    module,
    args: argparse.Namespace,
    h_value: float,
    physical_rms_fm: float,
    physical_radius_fm: float,
    scaling_mode: str,
    reuse_wavefunction_path: Path | None,
) -> StepResult:
    mu_minus = evaluate_mu_point(
        module, args, 1.0 - h_value, physical_rms_fm, physical_radius_fm, scaling_mode, reuse_wavefunction_path
    )
    mu_zero = evaluate_mu_point(
        module, args, 1.0, physical_rms_fm, physical_radius_fm, scaling_mode, reuse_wavefunction_path
    )
    mu_plus = evaluate_mu_point(
        module, args, 1.0 + h_value, physical_rms_fm, physical_radius_fm, scaling_mode, reuse_wavefunction_path
    )

    derivative_dw = (mu_plus.mu_times_delta_epsilon - mu_minus.mu_times_delta_epsilon) / (2.0 * h_value)
    delta_g_w = (4.0 / 3.0) * derivative_dw

    return StepResult(
        h=h_value,
        mu_minus=mu_minus,
        mu_zero=mu_zero,
        mu_plus=mu_plus,
        derivative_dw=derivative_dw,
        delta_g_w=delta_g_w,
    )


def print_mu_result(mu_result: MuPointResult) -> None:
    print(f"mu = {mu_result.mu:.10f}")
    print(f"Effective RMS radius passed to practise [fm]: {mu_result.effective_rms_fm:.10f}")
    print(f"Effective common weak radius passed to kernel [fm]: {mu_result.effective_radius_fm:.10f}")
    print(f"Wavefunction handling: {mu_result.wavefunction_mode}")
    print(f"Wavefunction file used: {mu_result.wavefunction_path}")
    print(f"Number of radial grid points: {mu_result.npoints}")
    print(f"min(r): {mu_result.r_min:.10e}")
    print(f"max(r): {mu_result.r_max:.10e}")
    print(f"Sampled normalization check int(G^2 + F^2) dr [trapz]: {mu_result.normalization:.16e}")
    if mu_result.used_unnormalized_diagnostic:
        print(
            "WARNING: proceeding with an unnormalized sampled wavefunction because "
            "--allow-unnormalized-diagnostic was set. This is a provenance diagnostic only."
        )
    print(
        "Raw sampled Delta E_W from f2py kernel before explicit mu^2 weak scaling: "
        f"{mu_result.raw_delta_epsilon:.16e}"
    )
    print(f"Explicit weak normalization factor applied for fixed physical G_F [mu^2]: {mu_result.weak_scale_factor:.16e}")
    print(f"Delta_epsilon_W(mu) in relativistic code units: {mu_result.delta_epsilon:.16e}")
    print(f"mu * Delta_epsilon_W(mu): {mu_result.mu_times_delta_epsilon:.16e}")
    print(f"Weak-kernel radius returned by Fortran [relativistic units]: {mu_result.r_rel:.16e}")


def relative_difference(reference_value: float, trial_value: float) -> float:
    if reference_value == 0.0:
        return math.inf if trial_value != 0.0 else 0.0
    return abs(trial_value - reference_value) / abs(reference_value)


def print_step_result(step_result: StepResult) -> None:
    print("")
    print(f"=== Step size h = {step_result.h:.6e} ===")
    print_mu_result(step_result.mu_minus)
    print("")
    print_mu_result(step_result.mu_zero)
    print("")
    print_mu_result(step_result.mu_plus)
    print("")
    print(f"Finite-difference derivative D_W = {step_result.derivative_dw:.16e}")
    print(f"delta_g_W(1s) = (4/3) * D_W = {step_result.delta_g_w:.16e}")


def print_step_size_stability(step_results: list[StepResult]) -> None:
    if len(step_results) < 2:
        return

    sorted_results = sorted(step_results, key=lambda item: item.h, reverse=True)
    reference = sorted_results[-1]
    comparison = sorted_results[0]
    derivative_rel_diff = relative_difference(reference.derivative_dw, comparison.derivative_dw)
    g_rel_diff = relative_difference(reference.delta_g_w, comparison.delta_g_w)

    print("")
    print("=== Step-size stability summary ===")
    for step_result in sorted(step_results, key=lambda item: item.h):
        print(
            f"h = {step_result.h:.6e}: D_W = {step_result.derivative_dw:.16e}, "
            f"delta_g_W = {step_result.delta_g_w:.16e}"
        )
    print(f"Relative difference in D_W between largest and smallest h: {derivative_rel_diff:.16e}")
    print(f"Relative difference in delta_g_W between largest and smallest h: {g_rel_diff:.16e}")
    if derivative_rel_diff < 1.0e-6:
        print("Stability verdict: stable to better than 1e-6 relative across tested h values.")
    elif derivative_rel_diff < 1.0e-4:
        print("Stability verdict: stable to better than 1e-4 relative across tested h values.")
    elif derivative_rel_diff < 1.0e-2:
        print("Stability verdict: stable to percent level.")
    else:
        print("Stability verdict: not stable at percent level.")


def main() -> int:
    args = parse_args()
    ensure_runtime_dependencies()
    h_values = validate_h_values(args.h_values)
    physical_rms_fm, physical_radius_fm, radius_input_mode = resolve_physical_radii(args)
    module, import_route = import_rad_kernel_suite()
    if not hasattr(module, "rad_kernel_f2py"):
        raise MassDerivativeError("Imported module does not expose `rad_kernel_f2py`.")

    standard_path = standard_wavefunction_path(args.znucl)
    step_results: list[StepResult] = []

    with preserve_standard_wavefunction(standard_path):
        reuse_wavefunction_path = None
        if args.wavefunction_mode == "reuse":
            reuse_wavefunction_path = prepare_reuse_wavefunction(
                args, physical_rms_fm, physical_radius_fm, args.radius_scaling_mode, standard_path
            )

        for h_value in h_values:
            step_results.append(
                evaluate_step_size(
                    module,
                    args,
                    h_value,
                    physical_rms_fm,
                    physical_radius_fm,
                    args.radius_scaling_mode,
                    reuse_wavefunction_path,
                )
            )

    print("Sampled-rad weak mass-derivative prototype")
    print("Physics route: virial / Karshenboim-Lee-Milstein mass derivative")
    print("Direct magnetic-response route status: not used here")
    print(f"Import route: {import_route}")
    print(f"Imported module: {module.__name__}")
    print(f"Module file: {getattr(module, '__file__', 'n/a')}")
    print(f"Fortran routine used for sampled weak shifts: {FORTRAN_ROUTINE_NAME}")
    print(f"Z: {args.znucl}")
    print(f"N: {args.nneutrons}")
    if radius_input_mode == "rms_fm":
        print(f"Physical RMS radius held fixed [fm]: {physical_rms_fm:.10f}")
        print(f"Derived physical common weak radius [fm]: {physical_radius_fm:.10f}")
    else:
        print(f"Physical common weak radius held fixed [fm]: {physical_radius_fm:.10f}")
        print(f"Derived physical RMS radius for wavefunction generation [fm]: {physical_rms_fm:.10f}")
    print(f"Radius-scaling mode: {args.radius_scaling_mode}")
    if args.radius_scaling_mode == "fixed-physical":
        print("Code-level radius rule: pass mu * radius_fm to emulate R_code(mu) = mu * R_code(1) at fixed physical radius.")
    else:
        print("WARNING: fixed-code-radius control selected. This is an incorrect derivative convention kept only as a diagnostic.")
    print(f"Wavefunction mode: {args.wavefunction_mode}")
    if args.wavefunction_mode == "reuse":
        print("WARNING: the same sampled 1s wavefunction is being reused for all mu values. This is an incomplete negative control.")
    else:
        print("Wavefunction provenance: regenerate mode uses practise-generated sampled 1s files from the copied rad workflow.")
    if args.allow_unnormalized_diagnostic:
        print("WARNING: unnormalized sampled-wavefunction diagnostics are enabled. Use this only for legacy provenance checks.")
    print("Weak normalization rule: the f2py kernel uses the built-in G_F code units, so this prototype multiplies each sampled shift by mu^2 externally.")
    print("Derivative formula used: D_W = [(1+h) Delta_epsilon_W(1+h) - (1-h) Delta_epsilon_W(1-h)] / (2h)")
    print("1s g-factor formula used: delta_g_W(1s) = (4/3) * D_W")

    for step_result in step_results:
        print_step_result(step_result)
    print_step_size_stability(step_results)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FirstOrderError, MassDerivativeError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
    except (FirstOrderError, MassDerivativeError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
