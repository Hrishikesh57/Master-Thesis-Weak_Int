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


def parse_driver_output(text: str) -> dict[str, np.ndarray]:
    gauge = None
    feynman_rows: list[list[float]] = []
    coulomb_rows: list[list[float]] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if "Results in Feynman gauge" in line:
            gauge = "feynman"
            continue
        if "Results in Coulomb gauge" in line:
            gauge = "coulomb"
            continue
        if line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) != 4:
            continue

        row = [float(value) for value in fields]
        if gauge == "feynman":
            feynman_rows.append(row)
        elif gauge == "coulomb":
            coulomb_rows.append(row)

    if not feynman_rows or not coulomb_rows:
        raise RuntimeError("Failed to parse both gauge tables from the reference driver output.")

    return {
        "feynman": np.asarray(feynman_rows, dtype=np.float64),
        "coulomb": np.asarray(coulomb_rows, dtype=np.float64),
    }


def compare_arrays(label: str, wrapped: np.ndarray, reference: np.ndarray) -> tuple[float, float]:
    abs_diff = np.abs(wrapped - reference)
    rel_diff = abs_diff / np.maximum(np.abs(reference), 1.0e-300)
    print(f"{label}: max_abs_diff={np.max(abs_diff):.16e} max_rel_diff={np.max(rel_diff):.16e}")
    return float(np.max(abs_diff)), float(np.max(rel_diff))


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the f2py Dirac wrapper against the original driver.")
    parser.add_argument("--z", type=int, default=92, help="Nuclear charge number Z.")
    parser.add_argument("--rms-fm", type=float, default=5.8571, help="Nuclear RMS radius in fm.")
    parser.add_argument("--order", type=int, default=9, help="B-spline order.")
    parser.add_argument(
        "--reference-exe",
        type=pathlib.Path,
        default=MODULE_DIR / "validation" / "bin" / "g_dirac_reference",
        help="Path to the compiled reference executable.",
    )
    args = parser.parse_args()

    import gfactor_lithium_dirac

    wrapper = gfactor_lithium_dirac.gfactor_lithium_dirac_f2py
    result = wrapper.compute_dirac_gfactor_tables(args.z, args.rms_fm, args.order)
    (
        npoints,
        irr_feynman,
        red_feynman,
        tot_feynman,
        irr_coulomb,
        red_coulomb,
        tot_coulomb,
        zalpha,
        rnucl,
        rcavity,
        ierr,
    ) = result

    if ierr != 0:
        raise RuntimeError(f"Fortran wrapper returned ierr={ierr}")

    completed = subprocess.run(
        [str(args.reference_exe)],
        input=f"{args.z}\n{args.rms_fm}\n",
        text=True,
        capture_output=True,
        check=True,
    )
    parsed = parse_driver_output(completed.stdout)

    ref_feynman = parsed["feynman"]
    ref_coulomb = parsed["coulomb"]

    ref_npoints_f = ref_feynman[:, 0].astype(np.int64)
    ref_npoints_c = ref_coulomb[:, 0].astype(np.int64)

    if not np.array_equal(npoints, ref_npoints_f) or not np.array_equal(npoints, ref_npoints_c):
        raise RuntimeError("npoints mismatch between wrapper and reference driver.")

    print(f"Validation case: Z={args.z} rms_fm={args.rms_fm:.6f} order={args.order}")
    print(f"Derived geometry: zalpha={zalpha:.16e} rnucl={rnucl:.16e} rcavity={rcavity:.16e}")

    compare_arrays("Feynman irr", irr_feynman, ref_feynman[:, 1])
    compare_arrays("Feynman red", red_feynman, ref_feynman[:, 2])
    compare_arrays("Feynman tot", tot_feynman, ref_feynman[:, 3])
    compare_arrays("Coulomb irr", irr_coulomb, ref_coulomb[:, 1])
    compare_arrays("Coulomb red", red_coulomb, ref_coulomb[:, 2])
    compare_arrays("Coulomb tot", tot_coulomb, ref_coulomb[:, 3])

    return 0


if __name__ == "__main__":
    sys.exit(main())
