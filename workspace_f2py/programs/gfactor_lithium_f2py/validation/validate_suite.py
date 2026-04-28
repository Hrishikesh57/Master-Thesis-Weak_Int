#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import numpy as np

MODULE_DIR = pathlib.Path(__file__).resolve().parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))


def run_reference(executable: pathlib.Path, z: int, rms_fm: float) -> str:
    completed = subprocess.run(
        [str(executable)],
        input=f"{z}\n{rms_fm}\n",
        text=True,
        capture_output=True,
        check=True,
    )
    return completed.stdout


def parse_tables(text: str, expected_columns: int, expect_coulomb: bool) -> dict[str, np.ndarray]:
    current = "feynman"
    rows: dict[str, list[list[float]]] = {"feynman": [], "coulomb": []}

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if "Results in Feynman gauge" in line:
            current = "feynman"
            continue
        if "Results in Coulomb gauge" in line:
            current = "coulomb"
            continue
        if line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) != expected_columns:
            continue
        rows[current].append([float(value) for value in fields])

    if not rows["feynman"]:
        raise RuntimeError("Failed to parse the Feynman table from the reference output.")
    if expect_coulomb and not rows["coulomb"]:
        raise RuntimeError("Failed to parse the Coulomb table from the reference output.")

    parsed = {"feynman": np.asarray(rows["feynman"], dtype=np.float64)}
    if expect_coulomb:
        parsed["coulomb"] = np.asarray(rows["coulomb"], dtype=np.float64)
    return parsed


def compare_columns(prefix: str, wrapped_columns: list[np.ndarray], reference: np.ndarray, labels: list[str]) -> None:
    ref_npoints = reference[:, 0].astype(np.int64)
    if not np.array_equal(wrapped_columns[0].astype(np.int64), ref_npoints):
        raise RuntimeError(f"{prefix}: npoints mismatch between wrapper and reference.")

    for idx, label in enumerate(labels, start=1):
        wrapped = np.asarray(wrapped_columns[idx], dtype=np.float64)
        ref = reference[:, idx]
        abs_diff = np.abs(wrapped - ref)
        rel_diff = abs_diff / np.maximum(np.abs(ref), 1.0e-300)
        print(f"{prefix} {label}: max_abs_diff={np.max(abs_diff):.16e} max_rel_diff={np.max(rel_diff):.16e}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the copied-workspace lithium f2py suite.")
    parser.add_argument("--z", type=int, default=92, help="Nuclear charge number Z.")
    parser.add_argument("--rms-fm", type=float, default=5.8571, help="Nuclear RMS radius in fm.")
    parser.add_argument("--order", type=int, default=9, help="B-spline order.")
    parser.add_argument(
        "--bin-dir",
        type=pathlib.Path,
        default=MODULE_DIR / "validation" / "bin",
        help="Directory containing the compiled copied-workspace reference executables.",
    )
    args = parser.parse_args()

    import gfactor_lithium_suite

    wrapper = gfactor_lithium_suite.gfactor_lithium_suite_f2py

    print(f"Validation case: Z={args.z} rms_fm={args.rms_fm:.6f} order={args.order}")

    result = wrapper.compute_dirac_gfactor_tables(args.z, args.rms_fm, args.order)
    npoints, irr_f, red_f, tot_f, irr_c, red_c, tot_c, zalpha, rnucl, rcavity, ierr = result
    if ierr != 0:
        raise RuntimeError(f"Dirac wrapper returned ierr={ierr}")
    print(f"Derived geometry: zalpha={zalpha:.16e} rnucl={rnucl:.16e} rcavity={rcavity:.16e}")
    parsed = parse_tables(run_reference(args.bin_dir / "g_dirac", args.z, args.rms_fm), 4, True)
    compare_columns("g_dirac Feynman", [npoints, irr_f, red_f, tot_f], parsed["feynman"], ["irr", "red", "tot"])
    compare_columns("g_dirac Coulomb", [npoints, irr_c, red_c, tot_c], parsed["coulomb"], ["irr", "red", "tot"])

    result = wrapper.compute_delbrueck_gfactor_tables(args.z, args.rms_fm, args.order)
    npoints, irr_f, red_f, tot_f, irr_c, red_c, tot_c, _, _, _, ierr = result
    if ierr != 0:
        raise RuntimeError(f"Delbrueck wrapper returned ierr={ierr}")
    parsed = parse_tables(run_reference(args.bin_dir / "g_delbrueck", args.z, args.rms_fm), 4, True)
    compare_columns("g_delbrueck Feynman", [npoints, irr_f, red_f, tot_f], parsed["feynman"], ["irr", "red", "tot"])
    compare_columns("g_delbrueck Coulomb", [npoints, irr_c, red_c, tot_c], parsed["coulomb"], ["irr", "red", "tot"])

    result = wrapper.compute_propagator_loop_feynman_table(args.z, args.rms_fm, args.order)
    npoints, irr, red, tot, _, _, _, ierr = result
    if ierr != 0:
        raise RuntimeError(f"Propagator-loop wrapper returned ierr={ierr}")
    parsed = parse_tables(run_reference(args.bin_dir / "g_propagator_loop", args.z, args.rms_fm), 4, False)
    compare_columns("g_propagator_loop Feynman", [npoints, irr, red, tot], parsed["feynman"], ["irr", "red", "tot"])

    result = wrapper.compute_uehling_gfactor_tables(args.z, args.rms_fm, args.order)
    if result[-1] != 0:
        raise RuntimeError(f"Uehling wrapper returned ierr={result[-1]}")
    npoints = result[0]
    parsed = parse_tables(run_reference(args.bin_dir / "g_uehling", args.z, args.rms_fm), 8, True)
    compare_columns(
        "g_uehling Feynman",
        [npoints, result[1], result[2], result[3], result[4], result[5], result[6], result[7]],
        parsed["feynman"],
        ["a", "b", "e", "g1", "g2", "g3", "tot"],
    )
    compare_columns(
        "g_uehling Coulomb",
        [npoints, result[8], result[9], result[10], result[11], result[12], result[13], result[14]],
        parsed["coulomb"],
        ["a", "b", "e", "g1", "g2", "g3", "tot"],
    )

    result = wrapper.compute_uehling_all_tables(args.z, args.rms_fm, args.order)
    if result[-1] != 0:
        raise RuntimeError(f"Uehling-all wrapper returned ierr={result[-1]}")
    npoints = result[0]
    parsed = parse_tables(run_reference(args.bin_dir / "g_uehling_all", args.z, args.rms_fm), 4, True)
    compare_columns(
        "g_uehling_all Feynman",
        [npoints, result[1], result[2], result[3]],
        parsed["feynman"],
        ["g_base", "g_with_uehling", "delta"],
    )
    compare_columns(
        "g_uehling_all Coulomb",
        [npoints, result[4], result[5], result[6]],
        parsed["coulomb"],
        ["g_base", "g_with_uehling", "delta"],
    )

    result = wrapper.compute_wichmann_kroll_tables(args.z, args.rms_fm, args.order)
    if result[-1] != 0:
        raise RuntimeError(f"Wichmann-Kroll wrapper returned ierr={result[-1]}")
    npoints = result[0]
    parsed = parse_tables(run_reference(args.bin_dir / "g_wichmann_kroll", args.z, args.rms_fm), 8, True)
    compare_columns(
        "g_wichmann_kroll Feynman",
        [npoints, result[1], result[2], result[3], result[4], result[5], result[6], result[7]],
        parsed["feynman"],
        ["a", "b", "e", "g1", "g2", "g3", "tot"],
    )
    compare_columns(
        "g_wichmann_kroll Coulomb",
        [npoints, result[8], result[9], result[10], result[11], result[12], result[13], result[14]],
        parsed["coulomb"],
        ["a", "b", "e", "g1", "g2", "g3", "tot"],
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
