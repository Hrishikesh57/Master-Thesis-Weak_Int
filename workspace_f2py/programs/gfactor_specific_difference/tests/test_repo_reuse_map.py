from __future__ import annotations

import unittest
from pathlib import Path

from _support import PROJECT_ROOT


class RepoReuseMapTest(unittest.TestCase):
    def test_documented_reuse_paths_exist(self) -> None:
        repo_root = PROJECT_ROOT.parents[2]
        expected_paths = [
            "workspace_f2py/codes/bspqed/src/bspqed_dkb.f90",
            "workspace_f2py/codes/bspqed/src/bspqed_discretization.f90",
            "workspace_f2py/codes/bspqed/src/bspqed_physconsts.f90",
            "workspace_f2py/codes/bspqed/src/bspqed_potentials.f90",
            "workspace_f2py/python_api/workspace_f2py_api/executables.py",
            "workspace_f2py/python_api/workspace_f2py_api/direct_wrappers.py",
            "workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90",
            "workspace_f2py/programs/rad_f2py/src/rad_kernel_f2py.f90",
        ]
        for relative_path in expected_paths:
            with self.subTest(path=relative_path):
                self.assertTrue((repo_root / relative_path).exists())


if __name__ == "__main__":
    unittest.main()
