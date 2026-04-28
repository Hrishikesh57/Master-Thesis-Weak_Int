from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAD_DIR = ROOT.parent / "rad"
sys.path.insert(0, str(ROOT))

import rad_kernel_suite


def read_sampled_wavefunction(path: Path) -> tuple[float | None, np.ndarray, np.ndarray, np.ndarray]:
    energy_rel = None
    r_vals: list[float] = []
    g_vals: list[float] = []
    f_vals: list[float] = []

    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("# state_energy_rel ="):
            energy_rel = float(stripped.split("=", 1)[1])
        parts = stripped.split()
        if len(parts) != 3:
            continue
        try:
            rval, gval, fval = (float(item) for item in parts)
        except ValueError:
            continue
        r_vals.append(rval)
        g_vals.append(gval)
        f_vals.append(fval)

    return energy_rel, np.array(r_vals), np.array(g_vals), np.array(f_vals)


def parse_last_float(pattern: str, text: str) -> float:
    match = re.search(pattern, text)
    if match is None:
        raise RuntimeError(f"Could not parse pattern: {pattern!r}")
    return float(match.group(1))


def run_first_order_reference() -> float:
    proc = subprocess.run(
        ["bin/something"],
        cwd=RAD_DIR,
        input="energy\n1\n92\n146\n5.8571\n",
        text=True,
        capture_output=True,
        check=True,
    )
    return parse_last_float(r"E_pot .*?=\s*([+-]?\d+\.\d+E[+-]\d+)", proc.stdout)


def run_amplitude_reference(e_initial: float, e_final: float, e_intermediate: float) -> tuple[float, float, float]:
    stdin = (
        "amplitudes\n"
        "55\n"
        "78\n"
        "4.8041\n"
        "2\n"
        "6s\n"
        "7s\n"
        "1\n"
        "1\n"
        "6p\n"
        f"{e_initial}\n"
        f"{e_final}\n"
        "-0.1\n"
        f"{e_intermediate}\n"
    )
    proc = subprocess.run(
        ["bin/something"],
        cwd=RAD_DIR,
        input=stdin,
        text=True,
        capture_output=True,
        check=True,
    )
    a_s = parse_last_float(r"Stark Amplitude A_S =\s*([+-]?\d+\.\d+(?:E[+-]\d+)?)", proc.stdout)
    a_w = parse_last_float(r"Weak Amplitude A_W =\s*([+-]?\d+\.\d+(?:E[+-]\d+)?)", proc.stdout)
    eta = parse_last_float(r"Asymmetry Factor .*?=\s*([+-]?\d+\.\d+(?:E[+-]\d+)?)", proc.stdout)
    return a_s, a_w, eta


def abs_rel_diff(value: float, reference: float) -> tuple[float, float]:
    abs_diff = abs(value - reference)
    rel_diff = abs_diff / abs(reference) if reference != 0.0 else float("nan")
    return abs_diff, rel_diff


def main() -> None:
    ref_energy_file = RAD_DIR / "data/reference/radial_functions/radial_function_relunits.txt"
    _, r_ref, g_ref, f_ref = read_sampled_wavefunction(ref_energy_file)

    first_order = rad_kernel_suite.rad_kernel_f2py.first_order_weak_shift(
        92, 146, 5.8571, r_ref, g_ref, f_ref
    )
    first_order_ref = run_first_order_reference()
    first_abs, first_rel = abs_rel_diff(first_order[0], first_order_ref)

    scf_dir = RAD_DIR / "outputs/wavefunctions/scf"
    e_initial, r_i, g_i, f_i = read_sampled_wavefunction(scf_dir / "55_scf_polynomial_table_6s.txt")
    e_intermediate, r_n, g_n, f_n = read_sampled_wavefunction(scf_dir / "55_scf_polynomial_table_6p.txt")
    e_final, r_f, g_f, f_f = read_sampled_wavefunction(scf_dir / "55_scf_polynomial_table_7s.txt")
    if e_initial is None or e_intermediate is None or e_final is None:
        raise RuntimeError("Missing `state_energy_rel` metadata in SCF wavefunction headers.")

    amp = rad_kernel_suite.rad_kernel_f2py.discrete_amp_sampled(
        55,
        78,
        4.8041,
        e_initial,
        e_final,
        e_intermediate,
        -0.1,
        1.0,
        1,
        r_i,
        g_i,
        f_i,
        r_n,
        g_n,
        f_n,
        r_f,
        g_f,
        f_f,
    )
    amp_ref = run_amplitude_reference(e_initial, e_final, e_intermediate)
    a_s_abs, a_s_rel = abs_rel_diff(amp[2], amp_ref[0])
    a_w_abs, a_w_rel = abs_rel_diff(amp[3], amp_ref[1])
    eta_abs, eta_rel = abs_rel_diff(amp[4], amp_ref[2])

    print("First-order weak shift")
    print("  wrapper:", first_order[0])
    print("  reference:", first_order_ref)
    print("  abs diff:", first_abs)
    print("  rel diff:", first_rel)
    print("")
    print("Discrete 6s -> 7s amplitude via 6p")
    print("  wrapper A_S:", amp[2])
    print("  reference A_S:", amp_ref[0])
    print("  A_S abs diff:", a_s_abs)
    print("  A_S rel diff:", a_s_rel)
    print("  wrapper A_W:", amp[3])
    print("  reference A_W:", amp_ref[1])
    print("  A_W abs diff:", a_w_abs)
    print("  A_W rel diff:", a_w_rel)
    print("  wrapper eta:", amp[4])
    print("  reference eta:", amp_ref[2])
    print("  eta abs diff:", eta_abs)
    print("  eta rel diff:", eta_rel)


if __name__ == "__main__":
    main()
