from __future__ import annotations

from pathlib import Path

from ._common import WORKSPACE_F2PY_ROOT, run_command


PROGRAM_TARGETS: dict[str, tuple[str, ...]] = {
    "gfactor": ("g_dirac", "g_delbrueck", "g_uehling", "g_uehling_all", "g_wichmann_kroll"),
    "gfactor_boron": ("g_dirac", "g_delbrueck", "g_propagator_loop", "g_uehling", "g_uehling_all", "g_wichmann_kroll"),
    "gfactor_lithium": ("g_dirac", "g_delbrueck", "g_propagator_loop", "g_uehling", "g_uehling_all", "g_wichmann_kroll"),
    "selfen": ("zero", "one", "many"),
    "hartree": ("dhfcb", "dhfc", "shf"),
    "rad": ("practise", "something", "theoretical_radial_functions_1s_2s", "bspwf_h"),
}

F2PY_BUILD_SCRIPTS: dict[str, Path] = {
    "gfactor_lithium": WORKSPACE_F2PY_ROOT / "programs" / "gfactor_lithium_f2py" / "build_f2py.sh",
    "selfen": WORKSPACE_F2PY_ROOT / "programs" / "selfen_f2py" / "build_f2py.sh",
    "hartree": WORKSPACE_F2PY_ROOT / "programs" / "hartree_f2py" / "build_f2py.sh",
    "rad": WORKSPACE_F2PY_ROOT / "programs" / "rad_f2py" / "build_f2py.sh",
}


def build_program_target(program_name: str, *targets: str) -> None:
    makefile_dir = WORKSPACE_F2PY_ROOT / "programs" / program_name
    if not makefile_dir.exists():
        raise FileNotFoundError(f"Unknown copied-workspace program family: {program_name}")
    run_command(["make", "-C", str(makefile_dir), *targets])


def build_lithium_f2py() -> None:
    build_f2py_target("gfactor_lithium")


def build_selfen_f2py() -> None:
    build_f2py_target("selfen")


def build_hartree_f2py() -> None:
    build_f2py_target("hartree")


def build_rad_f2py() -> None:
    build_f2py_target("rad")


def build_f2py_target(target_name: str) -> None:
    script = F2PY_BUILD_SCRIPTS.get(target_name)
    if script is None:
        raise ValueError(f"Unknown copied-workspace f2py target: {target_name}")
    run_command(["bash", str(script)])


def build_all_f2py() -> None:
    for target_name in ("gfactor_lithium", "selfen", "hartree", "rad"):
        build_f2py_target(target_name)


def build_python_backend() -> None:
    script = WORKSPACE_F2PY_ROOT / "build_python_backend.sh"
    run_command(["bash", str(script)])
