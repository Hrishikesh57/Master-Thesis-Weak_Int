from __future__ import annotations

import unittest

from _support import PROJECT_ROOT
from gfactor_specific_difference import run_lilike_dirac
from gfactor_specific_difference.reference_data import load_default_case
from gfactor_specific_difference.dirac_bspline_interface import run_lilike_correction


class LiLikeBackendSmokeTest(unittest.TestCase):
    def test_lilike_dirac_matches_selected_gauge_table(self) -> None:
        case = load_default_case("lilike_u")
        components = run_lilike_dirac(case)
        tables = run_lilike_correction(case, "dirac")
        selected = tables[case.backend.lithium_selected_gauge]
        self.assertAlmostEqual(components.dirac_value, selected.selected_value, places=14)
        self.assertGreater(len(selected.rows), 0)


if __name__ == "__main__":
    unittest.main()
