from __future__ import annotations

import importlib
import sys

from ._common import WORKSPACE_F2PY_ROOT


def load_lithium_suite_module():
    module_dir = WORKSPACE_F2PY_ROOT / "programs" / "gfactor_lithium_f2py"
    module_dir_str = str(module_dir)
    if module_dir_str not in sys.path:
        sys.path.insert(0, module_dir_str)
    return importlib.import_module("gfactor_lithium_suite")
