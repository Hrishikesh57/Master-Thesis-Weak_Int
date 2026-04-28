from __future__ import annotations

import unittest

from _support import PROJECT_ROOT
from gfactor_specific_difference import DEFAULT_CASE_FILES, DEFAULT_SPECIFIC_DIFFERENCE_FILES, load_default_case, load_named_run_config


class ConfigLoadingTest(unittest.TestCase):
    def test_default_case_configs_load(self) -> None:
        for name in DEFAULT_CASE_FILES:
            with self.subTest(case=name):
                case = load_default_case(name)
                self.assertGreater(case.znucl, 0)
                self.assertGreater(case.rms_fm, 0.0)
                self.assertIn(case.state_label, {"1s", "2s"})

    def test_specific_difference_configs_load(self) -> None:
        for path in DEFAULT_SPECIFIC_DIFFERENCE_FILES.values():
            with self.subTest(config=path.name):
                run_cfg = load_named_run_config(path)
                self.assertIn("hydrogen_case", run_cfg)
                self.assertIn("lithium_case", run_cfg)
                self.assertGreaterEqual(len(run_cfg["fns_rms_scan_fm"]), 3)


if __name__ == "__main__":
    unittest.main()
