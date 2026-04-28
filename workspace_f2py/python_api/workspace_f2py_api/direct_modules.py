from __future__ import annotations

import importlib
import sys

from ._common import WORKSPACE_F2PY_ROOT


def _load_module(module_dir_relative: str, module_name: str):
    module_dir = WORKSPACE_F2PY_ROOT / module_dir_relative
    module_dir_str = str(module_dir)
    if module_dir_str not in sys.path:
        sys.path.insert(0, module_dir_str)
    return importlib.import_module(module_name)


def load_lithium_suite_module():
    return _load_module("programs/gfactor_lithium_f2py", "gfactor_lithium_suite")


def load_selfen_suite_module():
    return _load_module("programs/selfen_f2py", "selfen_suite")


def load_hartree_suite_module():
    return _load_module("programs/hartree_f2py", "hartree_suite")


def load_rad_kernel_suite_module():
    return _load_module("programs/rad_f2py", "rad_kernel_suite")
