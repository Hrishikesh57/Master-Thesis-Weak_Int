from __future__ import annotations

import unittest

from _support import PROJECT_ROOT
from gfactor_specific_difference import run_hlike_dirac
from gfactor_specific_difference.reference_data import load_default_case
from gfactor_specific_difference.dirac_bspline_interface import run_hlike_correction


class HLikeBackendSmokeTest(unittest.TestCase):
    def test_hlike_dirac_matches_parsed_driver_output(self) -> None:
        case = load_default_case("hlike_u")
        components = run_hlike_dirac(case)
        direct_table = run_hlike_correction(case, "dirac")
        self.assertAlmostEqual(components.dirac_value, direct_table.selected_value, places=14)
        self.assertGreater(len(direct_table.rows), 0)


if __name__ == "__main__":
    unittest.main()
