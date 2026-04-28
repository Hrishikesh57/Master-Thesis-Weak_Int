#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PYTHON_API_DIR = SCRIPT_DIR.parent
if str(PYTHON_API_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_API_DIR))

from workspace_f2py_api import (  # noqa: E402
    build_lithium_f2py,
    load_lithium_suite_module,
    run_gfactor_hydrogen_program,
    run_selfen_zero,
)


def main() -> None:
    build_lithium_f2py()
    lithium_module = load_lithium_suite_module()
    wrapper = lithium_module.gfactor_lithium_suite_f2py

    npoints, irr_f, red_f, tot_f, irr_c, red_c, tot_c, zalpha, rnucl, rcavity, ierr = (
        wrapper.compute_dirac_gfactor_tables(92, 5.8571, 9)
    )
    if ierr != 0:
        raise RuntimeError(f"Lithium f2py call failed with ierr={ierr}")
    print("Lithium-like Dirac last row:", int(npoints[-1]), tot_f[-1], tot_c[-1])
    print("Lithium geometry:", zalpha, rnucl, rcavity)

    hydrogen = run_gfactor_hydrogen_program("dirac", 92, 5.8571, -1, 1)
    print("Hydrogen-like Dirac last row:", hydrogen["rows"][-1])

    zero = run_selfen_zero(92, 5.8571, -1, 1, 80, 1.0e-4, 1.0e-4)
    print("Zero-potential self-energy:", zero["result"])


if __name__ == "__main__":
    main()
