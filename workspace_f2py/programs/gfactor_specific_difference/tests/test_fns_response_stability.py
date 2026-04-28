from __future__ import annotations

import unittest

from _support import PROJECT_ROOT
from gfactor_specific_difference import assess_derivative_stability, estimate_fns_response
from gfactor_specific_difference.reference_data import load_default_case


class FnsResponseStabilityTest(unittest.TestCase):
    def test_hlike_fns_derivative_is_reasonably_stable_under_step_halving(self) -> None:
        case = load_default_case("hlike_u")
        coarse = estimate_fns_response(case, "hlike", [case.rms_fm - 0.01, case.rms_fm, case.rms_fm + 0.01])
        fine = estimate_fns_response(case, "hlike", [case.rms_fm - 0.005, case.rms_fm, case.rms_fm + 0.005])
        stability = assess_derivative_stability((coarse.derivative_value, fine.derivative_value), tolerance_rel=2.5e-1)
        self.assertTrue(stability.stable)


if __name__ == "__main__":
    unittest.main()
