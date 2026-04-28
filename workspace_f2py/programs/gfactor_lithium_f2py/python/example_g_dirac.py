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


def parse_reference_output(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rows: list[list[float]] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 4:
            continue
        rows.append([float(value) for value in fields])

    if not rows:
        raise RuntimeError("No numeric rows were parsed from the reference executable output.")

    data = np.asarray(rows, dtype=np.float64)
    return data[:, 0].astype(np.int64), data[:, 1], data[:, 2], data[:, 3]


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the lithium-like Dirac g-factor f2py wrapper.")
    parser.add_argument("--z", type=int, default=92, help="Nuclear charge number Z.")
    parser.add_argument("--rms-fm", type=float, default=5.8571, help="Nuclear RMS radius in fm.")
    parser.add_argument("--order", type=int, default=9, help="B-spline order passed to the wrapper.")
    parser.add_argument(
        "--reference-exe",
        type=pathlib.Path,
        default=MODULE_DIR / "validation" / "bin" / "g_dirac_reference",
        help="Optional path to the original Fortran g_dirac executable for comparison.",
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

    print(f"Z = {args.z}")
    print(f"rms_fm = {args.rms_fm:.6f}")
    print(f"order = {args.order}")
    print(f"zalpha = {zalpha:.16e}")
    print(f"rnucl = {rnucl:.16e}")
    print(f"rcavity = {rcavity:.16e}")
    print("Last-row totals:")
    print(f"  npoints = {int(npoints[-1])}")
    print(f"  Feynman total = {tot_feynman[-1]:.16e}")
    print(f"  Coulomb total = {tot_coulomb[-1]:.16e}")

    if args.reference_exe.is_file():
        completed = subprocess.run(
            [str(args.reference_exe)],
            input=f"{args.z}\n{args.rms_fm}\n",
            text=True,
            capture_output=True,
            check=True,
        )
        numeric_lines = []
        gauge = None
        feynman_lines = []
        coulomb_lines = []
        for raw_line in completed.stdout.splitlines():
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
            if gauge == "feynman":
                feynman_lines.append(line)
            elif gauge == "coulomb":
                coulomb_lines.append(line)

        ref_npoints_f, ref_irr_f, ref_red_f, ref_tot_f = parse_reference_output("\n".join(feynman_lines))
        ref_npoints_c, ref_irr_c, ref_red_c, ref_tot_c = parse_reference_output("\n".join(coulomb_lines))

        if not np.array_equal(npoints, ref_npoints_f) or not np.array_equal(npoints, ref_npoints_c):
            raise RuntimeError("Reference npoints grid does not match the wrapper output.")

        print("Reference comparison:")
        print(f"  max |Feynman total diff| = {np.max(np.abs(tot_feynman - ref_tot_f)):.16e}")
        print(f"  max |Coulomb total diff| = {np.max(np.abs(tot_coulomb - ref_tot_c)):.16e}")
    else:
        print(f"Reference executable not found at {args.reference_exe}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
