from __future__ import annotations

import unittest

from _support import PROJECT_ROOT
from gfactor_specific_difference import convert_r_fm_to_relativistic, electron_reduced_compton_fm, rms_fm_to_spherical_radius_relativistic
from gfactor_specific_difference.weak_contribution import summarize_weak_potential
from gfactor_specific_difference.reference_data import load_default_case


class UnitConventionTest(unittest.TestCase):
    def test_rms_to_rnucl_matches_driver_convention(self) -> None:
        rms_fm = 5.8571
        expected = ((5.0 / 3.0) ** 0.5) * rms_fm / electron_reduced_compton_fm()
        self.assertAlmostEqual(rms_fm_to_spherical_radius_relativistic(rms_fm), expected, places=14)

    def test_weak_radius_conversion_matches_rad_kernel(self) -> None:
        case = load_default_case("hlike_u")
        weak_summary = summarize_weak_potential(case)
        self.assertEqual(weak_summary.status, "ready")
        self.assertAlmostEqual(
            weak_summary.weak_potential_summary["r_rel"],
            convert_r_fm_to_relativistic(case.rms_fm),
            places=14,
        )


if __name__ == "__main__":
    unittest.main()
