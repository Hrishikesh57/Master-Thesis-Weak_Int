from __future__ import annotations

from pathlib import Path

from ._paths import CONFIG_DIR
from .cases import IonCase, load_ion_case


DEFAULT_CASE_FILES = {
    "hlike_pb": CONFIG_DIR / "hlike_pb.toml",
    "hlike_u": CONFIG_DIR / "hlike_u.toml",
    "lilike_pb": CONFIG_DIR / "lilike_pb.toml",
    "lilike_u": CONFIG_DIR / "lilike_u.toml",
}

DEFAULT_SPECIFIC_DIFFERENCE_FILES = {
    "pb": CONFIG_DIR / "specific_difference_pb.toml",
    "u": CONFIG_DIR / "specific_difference_u.toml",
}

DOCUMENTED_WEAK_TARGETS = {
    "m1_hlike_weak": {
        "status": "to_extract",
        "source": "Reading/Thesisvf2.pdf",
        "note": "A local documented weak target must be extracted before M1 is considered complete.",
    },
    "reduced_g_context": {
        "status": "context_only",
        "source": "Reading/Cakir_PhD_Thesis.pdf",
        "note": "Local reduced-g-factor context exists, but it is not yet a specific-difference target value.",
    },
}


def load_default_case(name: str) -> IonCase:
    path = DEFAULT_CASE_FILES.get(name)
    if path is None:
        raise KeyError(f"Unknown default case {name!r}.")
    return load_ion_case(path)


def list_default_case_paths() -> tuple[Path, ...]:
    return tuple(DEFAULT_CASE_FILES.values())
