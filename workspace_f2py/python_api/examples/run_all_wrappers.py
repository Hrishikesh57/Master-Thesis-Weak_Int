#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PYTHON_API_DIR = SCRIPT_DIR.parent
if str(PYTHON_API_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_API_DIR))

from workspace_f2py_api import (  # noqa: E402
    build_all_f2py,
    compute_hartree_shf,
    compute_lithium_dirac_tables,
    compute_rad_weak_potential_summary,
    compute_selfen_zero,
    run_selfen_zero,
)


def main() -> None:
    build_all_f2py()

    rad_summary = compute_rad_weak_potential_summary(92, 146, 5.8571)
    print("rad weak-potential summary:", rad_summary)

    zero = compute_selfen_zero(92, 5.8571, -1, 1, 80, 1.0e-4, 1.0e-4)
    print("selfen zero:", zero[:2], "ierr=", zero[-1])

    core_energy, valence_energy, total_energy, ierr = compute_hartree_shf(2, [1], [0], [2], [0])
    print("hartree shf:", core_energy, valence_energy, total_energy, "ierr=", ierr)

    npoints, irr_f, red_f, tot_f = compute_lithium_dirac_tables(92, 5.8571, 9)
    print("lithium dirac last row:", int(npoints[-1]), irr_f[-1], red_f[-1], tot_f[-1])

    zero_reference = run_selfen_zero(92, 5.8571, -1, 1, 80, 1.0e-4, 1.0e-4)
    print("selfen zero reference executable:", zero_reference["result"])


if __name__ == "__main__":
    main()
