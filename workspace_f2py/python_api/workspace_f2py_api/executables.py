from __future__ import annotations

from typing import Sequence

from ._common import WORKSPACE_F2PY_ROOT, format_stdin, numeric_rows, run_command, split_gauge_sections
from .build import build_program_target


def _single_feynman_section(stdout: str) -> str:
    marker_f = "# Results in Feynman gauge"
    if marker_f not in stdout:
        raise ValueError("Feynman-gauge marker was not found in the program output.")
    return stdout.split(marker_f, 1)[1]


def _run_executable(relative_executable: str, stdin_values: Sequence[object]) -> str:
    executable = WORKSPACE_F2PY_ROOT / relative_executable
    if not executable.exists():
        raise FileNotFoundError(f"Executable not found: {executable}")
    result = run_command([str(executable)], stdin_text=format_stdin(stdin_values))
    return result.stdout


def run_gfactor_hydrogen_program(workflow: str, znucl: int, rms_fm: float, kappa: int, nmain: int) -> dict[str, object]:
    build_program_target("gfactor", f"g_{workflow}")
    stdout = _run_executable(f"programs/gfactor/bin/g_{workflow}", (znucl, rms_fm, kappa, nmain))
    expected_columns = 4 if workflow == "uehling_all" else 2
    return {"stdout": stdout, "rows": numeric_rows(stdout, expected_columns)}


def run_gfactor_boron_program(workflow: str, znucl: int, rms_fm: float) -> dict[str, object]:
    build_program_target("gfactor_boron", f"g_{workflow}")
    stdout = _run_executable(f"programs/gfactor_boron/bin/g_{workflow}", (znucl, rms_fm))
    if workflow == "propagator_loop":
        return {"stdout": stdout, "feynman": numeric_rows(_single_feynman_section(stdout), 4)}
    expected_columns = 8 if workflow in {"uehling", "wichmann_kroll"} else 4
    feynman_text, coulomb_text = split_gauge_sections(stdout)
    return {
        "stdout": stdout,
        "feynman": numeric_rows(feynman_text, expected_columns),
        "coulomb": numeric_rows(coulomb_text, expected_columns),
    }


def run_gfactor_lithium_program(workflow: str, znucl: int, rms_fm: float) -> dict[str, object]:
    build_program_target("gfactor_lithium", f"g_{workflow}")
    stdout = _run_executable(f"programs/gfactor_lithium/bin/g_{workflow}", (znucl, rms_fm))
    if workflow == "propagator_loop":
        return {"stdout": stdout, "feynman": numeric_rows(_single_feynman_section(stdout), 4)}
    expected_columns = 8 if workflow in {"uehling", "wichmann_kroll"} else 4
    feynman_text, coulomb_text = split_gauge_sections(stdout)
    return {
        "stdout": stdout,
        "feynman": numeric_rows(feynman_text, expected_columns),
        "coulomb": numeric_rows(coulomb_text, expected_columns),
    }


def run_selfen_zero(znucl: int, rms_fm: float, kappa: int, nmain: int, npoints: int, epsabs: float, epsrel: float) -> dict[str, object]:
    build_program_target("selfen", "zero")
    stdout = _run_executable("programs/selfen/bin/zero", (znucl, rms_fm, kappa, nmain, npoints, epsabs, epsrel))
    rows = numeric_rows(stdout, 2)
    return {"stdout": stdout, "result": rows[-1] if rows else None}


def run_selfen_one(znucl: int, rms_fm: float, kappa: int, nmain: int, npoints: int, epsabs: float, epsrel: float) -> dict[str, object]:
    build_program_target("selfen", "one")
    stdout = _run_executable("programs/selfen/bin/one", (znucl, rms_fm, kappa, nmain, npoints, epsabs, epsrel))
    rows = numeric_rows(stdout, 2)
    return {"stdout": stdout, "result": rows[-1] if rows else None}


def run_selfen_many(
    znucl: int,
    rms_fm: float,
    kappa: int,
    nmain: int,
    npoints: int,
    epsabs: float,
    epsrel: float,
    modkapilow: int,
    modkapiupp: int,
) -> dict[str, object]:
    build_program_target("selfen", "many")
    stdout = _run_executable(
        "programs/selfen/bin/many",
        (znucl, rms_fm, kappa, nmain, npoints, epsabs, epsrel, modkapilow, modkapiupp),
    )
    return {"stdout": stdout, "rows": numeric_rows(stdout, 5)}


def run_hartree_program(program: str, stdin_values: Sequence[object]) -> dict[str, object]:
    if program not in {"dhfcb", "dhfc", "shf"}:
        raise ValueError(f"Unknown Hartree program: {program}")
    build_program_target("hartree", program)
    stdout = _run_executable(f"programs/hartree/bin/{program}", stdin_values)
    return {"stdout": stdout}


def run_rad_program(program: str, stdin_values: Sequence[object]) -> dict[str, object]:
    if program not in {"practise", "something", "theoretical_radial_functions_1s_2s", "bspwf_h"}:
        raise ValueError(f"Unknown RAD program: {program}")
    build_program_target("rad", program)
    stdout = _run_executable(f"programs/rad/bin/{program}", stdin_values)
    return {"stdout": stdout}
