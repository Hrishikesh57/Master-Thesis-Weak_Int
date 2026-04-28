#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from First_order import NORMALIZATION_TOLERANCE, import_rad_kernel_suite, normalization_trapz, read_sampled_wavefunction


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parent
RAD_DIR = WORKSPACE_F2PY_ROOT / "programs" / "rad"
RAD_RESULTS_DIR = RAD_DIR / "results"
RAD_F2PY_BUILD_SCRIPT = WORKSPACE_F2PY_ROOT / "programs" / "rad_f2py" / "build_f2py.sh"
PRACTISE_EXECUTABLE = RAD_DIR / "bin" / "practise"
WAVEFUNCTION_DIR = RAD_DIR / "outputs" / "wavefunctions" / "single_electron"

CSV_PATH = RAD_RESULTS_DIR / "weak_energy_corrections_basis_sums.csv"
LATEX_PATH = RAD_RESULTS_DIR / "weak_energy_corrections_basis_table.tex"
RUNLOG_PATH = RAD_RESULTS_DIR / "weak_energy_corrections_runlog.md"

ENERGY_REGEX = re.compile(r"Eigenvalue \(energy\) for chosen state =\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)")
FLOAT_REGEX = re.compile(r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)")

REQUESTED_CASES = (
    (55, 78, 4.8041, "1s"),
    (55, 78, 4.8041, "2s"),
    (87, 125, 5.5915, "1s"),
    (87, 125, 5.5915, "2s"),
    (92, 146, 5.8571, "1s"),
    (92, 146, 5.8571, "2s"),
)


class WeakEnergyError(RuntimeError):
    """User-facing failure for the basis-sum regeneration driver."""


@dataclass(frozen=True)
class RequestedCase:
    znucl: int
    nneutrons: int
    rms_fm: float
    ref_state: str

    @property
    def ref_nmain(self) -> int:
        return int(self.ref_state[:-1])

    @property
    def ref_kappa(self) -> int:
        if self.ref_state.endswith("s"):
            return -1
        raise WeakEnergyError(f"Unsupported reference state for this driver: {self.ref_state}")


@dataclass
class GeneratedReference:
    path: Path
    raw_eigenvalue_rel: float
    binding_energy_rel: float
    stdout: str


@dataclass
class BasisSpectrum:
    path: Path
    channel_label: str
    kappa: int
    order: int
    rcavity_rel: float
    units: str
    knots: np.ndarray
    spectrum: np.ndarray
    coefmat: np.ndarray

    @property
    def nspectrum(self) -> int:
        return int(self.spectrum.size)

    @property
    def nknots(self) -> int:
        return int(self.knots.size)

    @property
    def ncoeffs(self) -> int:
        return int(self.coefmat.shape[0])

    @property
    def nbsplines(self) -> int:
        return self.nknots - self.order


@dataclass
class ChannelResult:
    channel_label: str
    basis: BasisSpectrum
    delta_e_rel: float
    num_positive_available: int
    num_states_used: int
    excluded_basis_index: int
    energies_rel: np.ndarray
    matrix_elements: np.ndarray
    contributions: np.ndarray


@dataclass
class CaseResult:
    case: RequestedCase
    reference: GeneratedReference
    normalization: float
    delta_e1_pc_rel: float
    r_rel: float
    np12: ChannelResult
    np32: ChannelResult
    delta_e_total_rel: float
    notes: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate first-order parity-conserving weak shifts and second-order "
            "basis-sum parity-violating weak shifts for the requested H-like ions."
        )
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Reuse existing executables and the existing rad_f2py extension instead of rebuilding first.",
    )
    return parser.parse_args()


def ensure_results_dir() -> None:
    RAD_RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def format_command(cmd: list[str], cwd: Path | None = None, input_text: str | None = None) -> str:
    cmd_text = " ".join(_shell_quote(piece) for piece in cmd)
    lines = []
    if cwd is not None:
        lines.append(f"(cd {cwd} && {cmd_text})")
    else:
        lines.append(cmd_text)
    if input_text is not None:
        lines.append("stdin:")
        lines.append("```text")
        lines.append(input_text.rstrip())
        lines.append("```")
    return "\n".join(lines)


def _shell_quote(text: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_./:-]+", text):
        return text
    return "'" + text.replace("'", "'\"'\"'") + "'"


def run_logged_command(
    cmd: list[str],
    *,
    cwd: Path | None = None,
    input_text: str | None = None,
    log: list[dict[str, Any]],
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd is not None else None,
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )
    log.append(
        {
            "cmd": list(cmd),
            "cwd": cwd,
            "input": input_text,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        }
    )
    if completed.returncode != 0:
        lines = [
            "Command failed.",
            format_command(cmd, cwd=cwd, input_text=input_text),
            f"returncode = {completed.returncode}",
        ]
        stdout_tail = "\n".join(completed.stdout.strip().splitlines()[-20:])
        stderr_tail = "\n".join(completed.stderr.strip().splitlines()[-20:])
        if stdout_tail:
            lines.append("Last stdout lines:")
            lines.append(stdout_tail)
        if stderr_tail:
            lines.append("Last stderr lines:")
            lines.append(stderr_tail)
        raise WeakEnergyError("\n".join(lines))
    return completed


def build_required_targets(log: list[dict[str, Any]]) -> None:
    run_logged_command(["make", "-B", "-C", str(RAD_DIR), "practise", "something"], log=log)
    run_logged_command(["bash", str(RAD_F2PY_BUILD_SCRIPT)], cwd=WORKSPACE_F2PY_ROOT, log=log)


def get_git_hash(log: list[dict[str, Any]]) -> str:
    completed = run_logged_command(["git", "rev-parse", "--short", "HEAD"], cwd=WORKSPACE_F2PY_ROOT.parent, log=log)
    return completed.stdout.strip()


def get_gfortran_version(log: list[dict[str, Any]]) -> str:
    completed = run_logged_command(["gfortran", "--version"], log=log)
    return completed.stdout.strip()


def generate_reference_wavefunction(case: RequestedCase, log: list[dict[str, Any]]) -> GeneratedReference:
    input_text = f"1\n1\n{case.znucl}\n{case.ref_state}\n{case.rms_fm:.16f}\n"
    completed = run_logged_command([str(PRACTISE_EXECUTABLE)], cwd=RAD_DIR, input_text=input_text, log=log)
    matches = ENERGY_REGEX.findall(completed.stdout)
    if not matches:
        raise WeakEnergyError(
            "Could not parse the reference-state eigenvalue from `practise` output.\n"
            f"Case: Z={case.znucl}, state={case.ref_state}"
        )

    raw_eigenvalue_rel = float(matches[-1])
    binding_energy_rel = raw_eigenvalue_rel - 1.0
    path = WAVEFUNCTION_DIR / f"{case.znucl}_polynomial_table_{case.ref_state}.txt"
    if not path.is_file():
        raise WeakEnergyError(f"Expected reference wavefunction was not produced: {path}")
    return GeneratedReference(
        path=path,
        raw_eigenvalue_rel=raw_eigenvalue_rel,
        binding_energy_rel=binding_energy_rel,
        stdout=completed.stdout,
    )


def generate_basis_spectrum(znucl: int, rms_fm: float, channel_label: str, log: list[dict[str, Any]]) -> Path:
    input_text = f"1\n2\n{znucl}\n{channel_label}\n{rms_fm:.16f}\n"
    run_logged_command([str(PRACTISE_EXECUTABLE)], cwd=RAD_DIR, input_text=input_text, log=log)
    path = WAVEFUNCTION_DIR / f"{znucl}_spectrum_{channel_label}.txt"
    if not path.is_file():
        raise WeakEnergyError(f"Expected basis spectrum was not produced: {path}")
    return path


def read_basis_spectrum(path: Path) -> BasisSpectrum:
    lines = path.read_text().splitlines()
    metadata: dict[str, str] = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "knots_begin":
            break
        if line.startswith("#") and "=" in line:
            key, value = line[1:].split("=", 1)
            metadata[key.strip().lower()] = value.strip()
        i += 1

    if i >= len(lines) or lines[i].strip() != "knots_begin":
        raise WeakEnergyError(f"Malformed basis spectrum file, missing `knots_begin`: {path}")
    i += 1

    knots: list[float] = []
    while i < len(lines):
        line = lines[i].strip()
        if line == "knots_end":
            i += 1
            break
        parts = line.split()
        if len(parts) != 2:
            raise WeakEnergyError(f"Malformed knot entry in basis file: {path}\n{line}")
        knots.append(float(parts[1]))
        i += 1

    if i >= len(lines) or lines[i].strip() != "states_begin":
        raise WeakEnergyError(f"Malformed basis spectrum file, missing `states_begin`: {path}")
    i += 1

    spectrum: list[float] = []
    coefficient_columns: list[list[float]] = []
    while i < len(lines):
        line = lines[i].strip()
        if line == "states_end":
            break
        if line != "state_begin":
            raise WeakEnergyError(f"Malformed state block in basis file: {path}\n{line}")

        i += 1
        if i >= len(lines) or not lines[i].strip().startswith("index "):
            raise WeakEnergyError(f"Missing state index in basis file: {path}")
        i += 1
        if i >= len(lines) or not lines[i].strip().startswith("eigenvalue "):
            raise WeakEnergyError(f"Missing eigenvalue in basis file: {path}")
        eigenvalue_match = FLOAT_REGEX.search(lines[i])
        if eigenvalue_match is None:
            raise WeakEnergyError(f"Could not parse eigenvalue in basis file: {path}\n{lines[i]}")
        spectrum.append(float(eigenvalue_match.group(1)))
        i += 1

        if i >= len(lines) or lines[i].strip() != "coefficients_begin":
            raise WeakEnergyError(f"Missing `coefficients_begin` in basis file: {path}")
        i += 1

        coefficients: list[float] = []
        while i < len(lines):
            coeff_line = lines[i].strip()
            if coeff_line == "coefficients_end":
                i += 1
                break
            parts = coeff_line.split()
            if len(parts) != 2:
                raise WeakEnergyError(f"Malformed coefficient entry in basis file: {path}\n{coeff_line}")
            coefficients.append(float(parts[1]))
            i += 1

        if i >= len(lines) or lines[i].strip() != "state_end":
            raise WeakEnergyError(f"Missing `state_end` in basis file: {path}")
        i += 1
        coefficient_columns.append(coefficients)

    if not spectrum:
        raise WeakEnergyError(f"No basis states were parsed from: {path}")

    ncoeffs = len(coefficient_columns[0])
    if any(len(column) != ncoeffs for column in coefficient_columns):
        raise WeakEnergyError(f"Inconsistent coefficient counts across basis states: {path}")

    coefmat = np.asfortranarray(np.column_stack(coefficient_columns), dtype=np.float64)
    return BasisSpectrum(
        path=path,
        channel_label=metadata.get("channel", ""),
        kappa=int(metadata["kappa"]),
        order=int(metadata["order"]),
        rcavity_rel=float(metadata["rcavity"]),
        units=metadata.get("units", "relativistic_units"),
        knots=np.asarray(knots, dtype=np.float64),
        spectrum=np.asarray(spectrum, dtype=np.float64),
        coefmat=coefmat,
    )


def ensure_relativistic_single_electron_basis(basis: BasisSpectrum) -> None:
    normalized_units = basis.units.strip().lower()
    if normalized_units and "atomic" in normalized_units:
        raise WeakEnergyError(
            "The requested table uses H-like single-electron basis spectra in relativistic units. "
            f"Encountered atomic-unit basis metadata instead: {basis.path}"
        )


def compute_first_order_shift(module: Any, case: RequestedCase, reference_path: Path) -> tuple[float, float, float]:
    r, g, f = read_sampled_wavefunction(reference_path)
    normalization = normalization_trapz(r, g, f)
    if abs(normalization - 1.0) > NORMALIZATION_TOLERANCE:
        raise WeakEnergyError(
            "Reference wavefunction normalization check failed.\n"
            f"file = {reference_path}\n"
            f"integral = {normalization:.16e}\n"
            f"tolerance = {NORMALIZATION_TOLERANCE:.1e}"
        )
    delta_e1, r_rel, ierr = module.rad_kernel_f2py.first_order_weak_shift(
        case.znucl,
        case.nneutrons,
        case.rms_fm,
        np.asarray(r, dtype=np.float64),
        np.asarray(g, dtype=np.float64),
        np.asarray(f, dtype=np.float64),
    )
    if int(ierr) != 0:
        raise WeakEnergyError(f"first_order_weak_shift returned ierr={int(ierr)} for {reference_path}")
    return float(delta_e1), float(r_rel), float(normalization)


def compute_channel_sum(
    module: Any,
    case: RequestedCase,
    reference: GeneratedReference,
    basis: BasisSpectrum,
    max_positive_states: int = 0,
) -> ChannelResult:
    ensure_relativistic_single_electron_basis(basis)
    r_initial, g_initial, f_initial = read_sampled_wavefunction(reference.path)
    outputs = module.rad_kernel_f2py.second_order_basis_sum_rel(
        case.znucl,
        case.nneutrons,
        case.rms_fm,
        case.ref_nmain,
        case.ref_kappa,
        reference.binding_energy_rel,
        basis.kappa,
        basis.order,
        basis.rcavity_rel,
        max_positive_states,
        np.asarray(basis.knots, dtype=np.float64),
        np.asarray(basis.spectrum, dtype=np.float64),
        np.asfortranarray(basis.coefmat, dtype=np.float64),
        np.asarray(r_initial, dtype=np.float64),
        np.asarray(g_initial, dtype=np.float64),
        np.asarray(f_initial, dtype=np.float64),
    )
    (
        delta_e_total,
        num_positive_available,
        num_states_used,
        excluded_basis_index,
        energies_rel,
        matrix_elements,
        contributions,
        _r_rel,
        ierr,
    ) = outputs
    if int(ierr) != 0:
        raise WeakEnergyError(
            "second_order_basis_sum_rel failed.\n"
            f"case = Z={case.znucl}, state={case.ref_state}, channel={basis.channel_label}\n"
            f"ierr = {int(ierr)}"
        )
    return ChannelResult(
        channel_label=basis.channel_label,
        basis=basis,
        delta_e_rel=float(delta_e_total),
        num_positive_available=int(num_positive_available),
        num_states_used=int(num_states_used),
        excluded_basis_index=int(excluded_basis_index),
        energies_rel=np.asarray(energies_rel, dtype=np.float64),
        matrix_elements=np.asarray(matrix_elements, dtype=np.float64),
        contributions=np.asarray(contributions, dtype=np.float64),
    )


def discrete_2p_index(basis: BasisSpectrum) -> int:
    l_value = 1
    state_index_fortran = basis.nspectrum // 2 + (2 - l_value)
    return state_index_fortran - 1


def latex_scientific(value: float) -> str:
    if value == 0.0:
        return "$0.000000\\times 10^{0}$"
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10.0 ** exponent)
    return f"${mantissa:.6f}\\times 10^{{{exponent}}}$"


def build_notes(np12: ChannelResult, np32: ChannelResult) -> str:
    notes: list[str] = []
    if np32.num_positive_available > 0 and np32.num_states_used == 0 and abs(np32.delta_e_rel) == 0.0:
        notes.append(
            f"np3/2 basis file contained {np32.num_positive_available} positive-energy pseudostates, "
            "but the existing code-level weak angular rule suppresses kappa=-2 contributions exactly"
        )
    if np12.excluded_basis_index > 0:
        notes.append(f"reference-state index {np12.excluded_basis_index} excluded from kappa={np12.basis.kappa} basis channel")
    if not notes:
        notes.append("full positive-energy basis branch used")
    return "; ".join(notes)


def compute_requested_results(module: Any, log: list[dict[str, Any]]) -> tuple[list[CaseResult], dict[int, dict[str, Any]]]:
    generated_by_z: dict[int, dict[str, Any]] = {}
    results: list[CaseResult] = []

    for raw_case in REQUESTED_CASES:
        case = RequestedCase(*raw_case)
        cache = generated_by_z.setdefault(case.znucl, {"references": {}, "basis_paths": {}, "basis_data": {}})
        if case.ref_state not in cache["references"]:
            cache["references"][case.ref_state] = generate_reference_wavefunction(case, log)
        reference: GeneratedReference = cache["references"][case.ref_state]

        for channel_label in ("p", "p3"):
            if channel_label not in cache["basis_paths"]:
                cache["basis_paths"][channel_label] = generate_basis_spectrum(case.znucl, case.rms_fm, channel_label, log)
            if channel_label not in cache["basis_data"]:
                cache["basis_data"][channel_label] = read_basis_spectrum(cache["basis_paths"][channel_label])

        delta_e1_pc_rel, r_rel, normalization = compute_first_order_shift(module, case, reference.path)
        np12 = compute_channel_sum(module, case, reference, cache["basis_data"]["p"])
        np32 = compute_channel_sum(module, case, reference, cache["basis_data"]["p3"])
        delta_e_total_rel = delta_e1_pc_rel + np12.delta_e_rel + np32.delta_e_rel

        results.append(
            CaseResult(
                case=case,
                reference=reference,
                normalization=normalization,
                delta_e1_pc_rel=delta_e1_pc_rel,
                r_rel=r_rel,
                np12=np12,
                np32=np32,
                delta_e_total_rel=delta_e_total_rel,
                notes=build_notes(np12, np32),
            )
        )

    return results, generated_by_z


def write_csv(results: list[CaseResult]) -> None:
    basis = results[0].np12.basis
    basis_label = f"single_electron_DKB_order{basis.order}_nspectrum{basis.nspectrum}"
    header = (
        "Z,N,rms_fm,ref_state,DeltaE1_PC_rel,DeltaE2_np12_rel,DeltaE2_np32_rel,"
        "DeltaE_total_rel,num_np12_states,num_np32_states,basis_label,notes"
    )
    lines = [header]
    for result in results:
        row = [
            str(result.case.znucl),
            str(result.case.nneutrons),
            f"{result.case.rms_fm:.4f}",
            result.case.ref_state,
            f"{result.delta_e1_pc_rel:.16e}",
            f"{result.np12.delta_e_rel:.16e}",
            f"{result.np32.delta_e_rel:.16e}",
            f"{result.delta_e_total_rel:.16e}",
            str(result.np12.num_states_used),
            str(result.np32.num_states_used),
            basis_label,
            result.notes.replace(",", ";"),
        ]
        lines.append(",".join(row))
    CSV_PATH.write_text("\n".join(lines) + "\n")


def write_latex_table(results: list[CaseResult]) -> None:
    lines = [
        "\\begin{table}[ht]",
        "\\centering",
        "\\resizebox{1.0\\textwidth}{!}{",
        "\\renewcommand{\\arraystretch}{1.2}",
        "\\begin{tabular}{c c c c c c}",
        "\\toprule",
        "\\textbf{Z} & \\textbf{Ref. state} & $\\Delta E^{(1)}_{\\rm PC}$ [rel]",
        "& $\\Delta E^{(2)}_{np_{1/2}}$ [rel]",
        "& $\\Delta E^{(2)}_{np_{3/2}}$ [rel]",
        "& $\\Delta E_{\\rm tot}$ [rel] \\\\",
        "\\midrule",
    ]
    for result in results:
        lines.append(
            f"{result.case.znucl} & {result.case.ref_state} & "
            f"{latex_scientific(result.delta_e1_pc_rel)} & "
            f"{latex_scientific(result.np12.delta_e_rel)} & "
            f"{latex_scientific(result.np32.delta_e_rel)} & "
            f"{latex_scientific(result.delta_e_total_rel)} \\\\"
        )
    lines.extend(
        [
            "\\bottomrule",
            "\\end{tabular}",
            "}",
            "\\caption{First-order parity-conserving weak energy shifts and second-order parity-violating basis-sum corrections for hydrogen-like ions. The second-order terms are evaluated over the finite-cavity B-spline pseudostate spectra in the indicated $np_j$ channels.}",
            "\\label{tab:weak_energy_corrections_basis}",
            "\\end{table}",
        ]
    )
    LATEX_PATH.write_text("\n".join(lines) + "\n")


def top_contribution_lines(channel: ChannelResult, limit: int = 5) -> list[str]:
    nonzero_indices = np.flatnonzero(np.abs(channel.contributions) > 0.0)
    if nonzero_indices.size == 0:
        return ["No nonzero intermediate-state contributions were produced for this channel."]

    ranked = sorted(
        nonzero_indices.tolist(),
        key=lambda idx: abs(float(channel.contributions[idx])),
        reverse=True,
    )[:limit]
    lines = []
    for rank, idx in enumerate(ranked, start=1):
        lines.append(
            f"{rank}. basis_index={idx + 1}, E_n={channel.energies_rel[idx]:+.16e}, "
            f"M_nref={channel.matrix_elements[idx]:+.16e}, contribution={channel.contributions[idx]:+.16e}"
        )
    return lines


def validate_results(
    module: Any,
    results: list[CaseResult],
    generated_by_z: dict[int, dict[str, Any]],
) -> dict[str, Any]:
    lookup = {(result.case.znucl, result.case.ref_state): result for result in results}
    target = lookup[(92, "1s")]
    rerun_np12 = compute_channel_sum(module, target.case, target.reference, generated_by_z[92]["basis_data"]["p"])
    rerun_np32 = compute_channel_sum(module, target.case, target.reference, generated_by_z[92]["basis_data"]["p3"])
    rerun_total = target.delta_e1_pc_rel + rerun_np12.delta_e_rel + rerun_np32.delta_e_rel

    deterministic_ok = (
        abs(rerun_np12.delta_e_rel - target.np12.delta_e_rel) <= 1.0e-24
        and abs(rerun_np32.delta_e_rel - target.np32.delta_e_rel) <= 1.0e-30
        and abs(rerun_total - target.delta_e_total_rel) <= 1.0e-24
    )

    discrete_idx = discrete_2p_index(target.np12.basis)
    discrete_2p_contribution = float(target.np12.contributions[discrete_idx])
    full_basis_contribution = target.np12.delta_e_rel

    truncated_count = max(1, target.np12.num_states_used - 100)
    truncated_np12 = compute_channel_sum(
        module,
        target.case,
        target.reference,
        generated_by_z[92]["basis_data"]["p"],
        max_positive_states=truncated_count,
    )

    return {
        "deterministic_ok": deterministic_ok,
        "rerun_np12_rel": rerun_np12.delta_e_rel,
        "rerun_np32_rel": rerun_np32.delta_e_rel,
        "rerun_total_rel": rerun_total,
        "discrete_2p_only_rel": discrete_2p_contribution,
        "full_basis_np12_rel": full_basis_contribution,
        "truncated_positive_states": truncated_count,
        "truncated_np12_rel": truncated_np12.delta_e_rel,
        "top_contributions": top_contribution_lines(target.np12),
    }


def write_runlog(
    results: list[CaseResult],
    log: list[dict[str, Any]],
    git_hash: str,
    gfortran_version: str,
    validation: dict[str, Any],
) -> None:
    basis = results[0].np12.basis
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S %Z")
    lines = [
        "# Weak Energy Corrections Basis-Sum Runlog",
        "",
        f"- Generated at: {now}",
        f"- Repository commit: `{git_hash}`",
        "- Compiler:",
        "```text",
        gfortran_version,
        "```",
        "",
        "## Build And Run Commands",
        "",
    ]

    for entry in log:
        lines.append("```text")
        lines.append(format_command(entry["cmd"], cwd=entry["cwd"], input_text=entry["input"]))
        lines.append("```")
        lines.append(f"- returncode = {entry['returncode']}")
        if entry["stdout"].strip():
            lines.append("- stdout tail:")
            lines.append("```text")
            lines.extend(entry["stdout"].strip().splitlines()[-12:])
            lines.append("```")
        if entry["stderr"].strip():
            lines.append("- stderr tail:")
            lines.append("```text")
            lines.extend(entry["stderr"].strip().splitlines()[-12:])
            lines.append("```")
        lines.append("")

    lines.extend(
        [
            "## Basis Parameters",
            "",
            "- Single-electron generation path: `workspace_f2py/programs/rad/src/practise.f90`",
            "- Reference wavefunction mode: `practise` single-electron sampled export (`wavefunction_type=1`, `output_data_type=1`)",
            "- Basis-spectrum mode: `practise` single-electron basis export (`wavefunction_type=1`, `output_data_type=2`)",
            f"- Spline order: `{basis.order}`",
            f"- Number of knots: `{basis.nknots}`",
            f"- Number of B-splines: `{basis.nbsplines}`",
            f"- Total basis eigenstates per channel: `{basis.nspectrum}`",
            f"- Basis cavity radius from spectrum metadata [rel]: `{basis.rcavity_rel:.16e}`",
            "- Knot-generation convention from source: `kntseq1(rnucl, rcavity, npoints/10, order, knots)` with `npoints=1000`, so `ninside=100`",
            "- Weak-integration radius convention: `R_rel = rms_fm / electron_reduced_compton` and denominators use binding energies in `m_e = c = \\hbar = 1` units",
            "",
            "## Per-Case State Counts",
            "",
        ]
    )

    for result in results:
        lines.append(
            f"- `Z={result.case.znucl}`, `{result.case.ref_state}`:"
            f" normalization = `{result.normalization:.16e}`;"
            f" E_ref(raw) = `{result.reference.raw_eigenvalue_rel:.16e}`;"
            f" E_ref(binding) = `{result.reference.binding_energy_rel:.16e}`;"
            f" np1/2 used = `{result.np12.num_states_used}` of `{result.np12.num_positive_available}` positive states;"
            f" np3/2 used = `{result.np32.num_states_used}` of `{result.np32.num_positive_available}` positive states"
        )

    lines.extend(
        [
            "",
            "## Validation",
            "",
            f"- Deterministic rerun for `Z=92`, `1s`: `{validation['deterministic_ok']}`",
            f"- Rerun `DeltaE2_np12_rel`: `{validation['rerun_np12_rel']:.16e}`",
            f"- Rerun `DeltaE2_np32_rel`: `{validation['rerun_np32_rel']:.16e}`",
            f"- Rerun total correction: `{validation['rerun_total_rel']:.16e}`",
            f"- Discrete `2p`-only diagnostic contribution [rel]: `{validation['discrete_2p_only_rel']:.16e}`",
            f"- Full `np1/2` basis-sum contribution [rel]: `{validation['full_basis_np12_rel']:.16e}`",
            f"- Truncated `np1/2` cutoff diagnostic with first `{validation['truncated_positive_states']}` positive states [rel]: `{validation['truncated_np12_rel']:.16e}`",
            "- Largest `Z=92`, `1s`, `np1/2` intermediate-state contributions:",
        ]
    )
    lines.extend([f"  {line}" for line in validation["top_contributions"]])

    warnings = [
        result.notes for result in results if "np3/2 basis file contained" in result.notes
    ]
    if warnings:
        lines.extend(
            [
                "",
                "## Warnings",
                "",
            ]
        )
        lines.extend([f"- {warning}" for warning in warnings])

    RUNLOG_PATH.write_text("\n".join(lines) + "\n")


def print_final_table(results: list[CaseResult]) -> str:
    lines = []
    for result in results:
        lines.append(
            f"Z={result.case.znucl:2d}  ref={result.case.ref_state:>2s}  "
            f"DeltaE1_PC={result.delta_e1_pc_rel:+.16e}  "
            f"DeltaE2_np12={result.np12.delta_e_rel:+.16e}  "
            f"DeltaE2_np32={result.np32.delta_e_rel:+.16e}  "
            f"DeltaE_tot={result.delta_e_total_rel:+.16e}"
        )
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    ensure_results_dir()

    command_log: list[dict[str, Any]] = []
    git_hash = get_git_hash(command_log)
    gfortran_version = get_gfortran_version(command_log)

    if not args.skip_build:
        build_required_targets(command_log)

    module, import_route = import_rad_kernel_suite()
    if not hasattr(module.rad_kernel_f2py, "second_order_basis_sum_rel"):
        raise WeakEnergyError(
            "The imported rad_kernel_suite module does not expose `second_order_basis_sum_rel`. "
            "Rebuild with the updated `workspace_f2py/programs/rad_f2py/build_f2py.sh`."
        )

    results, generated_by_z = compute_requested_results(module, command_log)
    validation = validate_results(module, results, generated_by_z)

    write_csv(results)
    write_latex_table(results)
    write_runlog(results, command_log, git_hash, gfortran_version, validation)

    print("Imported rad kernel module via:", import_route)
    print("CSV written to:", CSV_PATH)
    print("LaTeX table written to:", LATEX_PATH)
    print("Runlog written to:", RUNLOG_PATH)
    print("")
    print("Final numerical table")
    print(print_final_table(results))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except WeakEnergyError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
