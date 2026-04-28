#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import sys

MODULE_DIR = pathlib.Path(__file__).resolve().parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one copied-workspace lithium g-factor workflow through f2py.")
    parser.add_argument(
        "--workflow",
        choices=["dirac", "delbrueck", "propagator_loop", "uehling", "uehling_all", "wichmann_kroll"],
        default="dirac",
        help="Which lithium workflow to evaluate.",
    )
    parser.add_argument("--z", type=int, default=92, help="Nuclear charge number Z.")
    parser.add_argument("--rms-fm", type=float, default=5.8571, help="Nuclear RMS radius in fm.")
    parser.add_argument("--order", type=int, default=9, help="B-spline order.")
    args = parser.parse_args()

    import gfactor_lithium_suite

    wrapper = gfactor_lithium_suite.gfactor_lithium_suite_f2py

    print(f"workflow = {args.workflow}")
    print(f"Z = {args.z}")
    print(f"rms_fm = {args.rms_fm:.6f}")
    print(f"order = {args.order}")

    if args.workflow == "dirac":
        result = wrapper.compute_dirac_gfactor_tables(args.z, args.rms_fm, args.order)
        npoints, _, _, tot_f, _, _, tot_c, zalpha, rnucl, rcavity, ierr = result
        print(f"zalpha = {zalpha:.16e}")
        print(f"rnucl = {rnucl:.16e}")
        print(f"rcavity = {rcavity:.16e}")
        print(f"last_npoints = {int(npoints[-1])}")
        print(f"feynman_total = {tot_f[-1]:.16e}")
        print(f"coulomb_total = {tot_c[-1]:.16e}")
        print(f"ierr = {ierr}")
        return 0

    if args.workflow == "delbrueck":
        result = wrapper.compute_delbrueck_gfactor_tables(args.z, args.rms_fm, args.order)
        npoints, _, _, tot_f, _, _, tot_c, zalpha, rnucl, rcavity, ierr = result
        print(f"zalpha = {zalpha:.16e}")
        print(f"rnucl = {rnucl:.16e}")
        print(f"rcavity = {rcavity:.16e}")
        print(f"last_npoints = {int(npoints[-1])}")
        print(f"feynman_total = {tot_f[-1]:.16e}")
        print(f"coulomb_total = {tot_c[-1]:.16e}")
        print(f"ierr = {ierr}")
        return 0

    if args.workflow == "propagator_loop":
        result = wrapper.compute_propagator_loop_feynman_table(args.z, args.rms_fm, args.order)
        npoints, _, _, tot, zalpha, rnucl, rcavity, ierr = result
        print(f"zalpha = {zalpha:.16e}")
        print(f"rnucl = {rnucl:.16e}")
        print(f"rcavity = {rcavity:.16e}")
        print(f"first_npoints = {int(npoints[0])}")
        print(f"last_npoints = {int(npoints[-1])}")
        print(f"feynman_total = {tot[-1]:.16e}")
        print(f"ierr = {ierr}")
        return 0

    if args.workflow == "uehling":
        result = wrapper.compute_uehling_gfactor_tables(args.z, args.rms_fm, args.order)
        npoints = result[0]
        tot_f = result[7]
        tot_c = result[14]
        zalpha, rnucl, rcavity, ierr = result[15:]
        print(f"zalpha = {zalpha:.16e}")
        print(f"rnucl = {rnucl:.16e}")
        print(f"rcavity = {rcavity:.16e}")
        print(f"last_npoints = {int(npoints[-1])}")
        print(f"feynman_total = {tot_f[-1]:.16e}")
        print(f"coulomb_total = {tot_c[-1]:.16e}")
        print(f"ierr = {ierr}")
        return 0

    if args.workflow == "uehling_all":
        result = wrapper.compute_uehling_all_tables(args.z, args.rms_fm, args.order)
        npoints = result[0]
        g_base_f, g_with_f, delta_f, g_base_c, g_with_c, delta_c, zalpha, rnucl, rcavity, ierr = result[1:]
        print(f"zalpha = {zalpha:.16e}")
        print(f"rnucl = {rnucl:.16e}")
        print(f"rcavity = {rcavity:.16e}")
        print(f"first_npoints = {int(npoints[0])}")
        print(f"last_npoints = {int(npoints[-1])}")
        print(f"feynman_base = {g_base_f[-1]:.16e}")
        print(f"feynman_with_uehling = {g_with_f[-1]:.16e}")
        print(f"feynman_delta = {delta_f[-1]:.16e}")
        print(f"coulomb_base = {g_base_c[-1]:.16e}")
        print(f"coulomb_with_uehling = {g_with_c[-1]:.16e}")
        print(f"coulomb_delta = {delta_c[-1]:.16e}")
        print(f"ierr = {ierr}")
        return 0

    result = wrapper.compute_wichmann_kroll_tables(args.z, args.rms_fm, args.order)
    npoints = result[0]
    tot_f = result[7]
    tot_c = result[14]
    zalpha, rnucl, rcavity, ierr = result[15:]
    print(f"zalpha = {zalpha:.16e}")
    print(f"rnucl = {rnucl:.16e}")
    print(f"rcavity = {rcavity:.16e}")
    print(f"last_npoints = {int(npoints[-1])}")
    print(f"feynman_total = {tot_f[-1]:.16e}")
    print(f"coulomb_total = {tot_c[-1]:.16e}")
    print(f"ierr = {ierr}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
