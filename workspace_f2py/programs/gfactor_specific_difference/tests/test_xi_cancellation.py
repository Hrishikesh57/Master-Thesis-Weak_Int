from __future__ import annotations

import unittest

from _support import PROJECT_ROOT
from gfactor_specific_difference import determine_xi, estimate_fns_response
from gfactor_specific_difference.reference_data import load_default_case


class XiCancellationTest(unittest.TestCase):
    def test_xi_cancels_the_leading_fns_derivative_by_construction(self) -> None:
        h_case = load_default_case("hlike_u")
        li_case = load_default_case("lilike_u")
        grid = [h_case.rms_fm - 0.01, h_case.rms_fm, h_case.rms_fm + 0.01]
        fns_h = estimate_fns_response(h_case, "hlike", grid)
        fns_li = estimate_fns_response(li_case, "lilike", grid)
        xi = determine_xi(li_case, h_case, fns_li, fns_h)
        cancellation = fns_li.derivative_value - xi.value * fns_h.derivative_value
        scale = max(abs(fns_li.derivative_value), 1.0)
        self.assertLess(abs(cancellation) / scale, 1.0e-12)


if __name__ == "__main__":
    unittest.main()
