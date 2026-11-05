"""Physics and regression checks; run after a clean build with unittest."""
from dataclasses import replace
import subprocess
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad

from spherical import Calculation, FermiDensity, Grid, NuclearInputs, HERE


class SphericalChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        profile = FermiDensity(6.624, 0.523)
        cls.calc = Calculation(NuclearInputs(82, 126, profile, profile, profile,
                                            'Explicit common-profile numerical test'), Grid(npoints=140))
        cls.s, cls.p = cls.calc.basis(-1), cls.calc.basis(1)
        cls.i, cls.j = cls.s.state(1), cls.p.state(2)

    def test_normalization_and_orthogonality(self):
        s = self.s
        ids = [0, self.i, s.state(2), len(s.energies)-1]
        overlap = s.g[:, ids].T @ (s.w[:, None]*s.g[:, ids]) + s.f[:, ids].T @ (s.w[:, None]*s.f[:, ids])
        np.testing.assert_allclose(overlap, np.eye(len(ids)), atol=2e-11)

    def test_density_normalization_weak_charge_and_exact_rms(self):
        c, b = self.calc, self.s
        for name in ('charge', 'proton', 'neutron'):
            profile = getattr(c.nucleus, name)
            integral = 4*np.pi*b.w @ (b.r**2*profile.density(b.r, c.lambda_fm))
            self.assertAlmostEqual(integral, 1.0, delta=2e-10)
        q = 4*np.pi*b.w @ (b.r**2*c.source(b))/(c.gf/(2*np.sqrt(2)))
        self.assertAlmostEqual(q, c.cvp*82-126, delta=3e-8)
        target = FermiDensity.from_rms(5.5012, 0.523)
        self.assertAlmostEqual(target.rms_fm, 5.5012, delta=1e-11)
        with self.assertRaises(ValueError):
            FermiDensity.from_rms(0.1, 0.523)

    def test_fermi_coulomb_against_independent_poisson_integral(self):
        c = self.calc; profile = c.nucleus.charge
        for r in (0.0, 1e-5, 3.0, 6.624, 15.0, 50.0):
            inside = 0 if r == 0 else quad(lambda x: x*x*profile.physical_density(x), 0, r,
                                          epsabs=1e-13)[0]/r
            # Finite upper endpoint here leaves an exponentially negligible tail.
            outside = quad(lambda x: x*profile.physical_density(x), r, max(r,profile.c_fm)+80*profile.a_fm,
                           epsabs=1e-13)[0]
            reference = -4*np.pi*c.nucleus.z*c.alpha*c.lambda_fm*(inside+outside)
            actual = c.kernel.charge_potential(82, profile.c_fm, profile.a_fm, r)
            self.assertAlmostEqual(actual/reference, 1.0, delta=2e-10)

    def test_pv_hermiticity_selection_and_electron_coupling(self):
        c, s, p, i, j = self.calc, self.s, self.p, self.i, self.j
        forward = c.pv_column(p, s, i)[j]
        reverse = c.pv_column(s, p, j)[i]
        self.assertAlmostEqual(abs(forward-reverse.conjugate()), 0, delta=1e-25)
        np.testing.assert_array_equal(c.pv_column(s, s, i), 0)
        pc = c.pc_column(s, i)[i]
        cve = c.cve
        try:
            c.cve = 2*cve
            self.assertEqual(c.pv_column(p, s, i)[j], forward)
            self.assertAlmostEqual(c.pc_column(s, i)[i]/pc, 2, delta=1e-13)
        finally:
            c.cve = cve

    def test_operator_against_explicit_four_spinor_angular_integration(self):
        c, s, p, i, j = self.calc, self.s, self.p, self.i, self.j
        mus, angular_weights = np.polynomial.legendre.leggauss(4)
        direct_weak, direct_dipole = 0j, 0j
        for mu, weight in zip(mus, angular_weights):
            # Condon-Shortley Omega_{-1,1/2} and Omega_{+1,1/2}, phi=0.
            os = np.array([1, 0])/np.sqrt(4*np.pi)
            op = np.array([-mu, -np.sqrt(1-mu*mu)])/np.sqrt(4*np.pi)
            psi_s = np.column_stack((s.g[:, i, None]*os, 1j*s.f[:, i, None]*op))
            psi_p = np.column_stack((p.g[:, j, None]*op, 1j*p.f[:, j, None]*os))
            gamma5_s = psi_s[:, [2, 3, 0, 1]]
            direct_weak += 2*np.pi*weight*np.sum(-s.w*c.source(s)*np.sum(psi_p.conj()*gamma5_s, axis=1))
            direct_dipole += 2*np.pi*weight*np.sum(-s.w*s.r*mu*np.sum(psi_p.conj()*psi_s, axis=1))
        computed = c.pv_column(p, s, i)[j]
        self.assertAlmostEqual(abs(direct_weak/computed-1), 0, delta=2e-13)
        radial_dipole = s.w @ (s.r/3*(p.g[:, j]*s.g[:, i]+p.f[:, j]*s.f[:, i]))
        self.assertAlmostEqual(abs(direct_dipole/radial_dipole-1), 0, delta=2e-13)

    def test_weak_profile_does_not_modify_charge_orbitals(self):
        c = self.calc
        changed = Calculation(replace(c.nucleus, neutron=FermiDensity(6.9, 0.55)), c.grid)
        b = changed.basis()
        np.testing.assert_array_equal(b.energies, self.s.energies)
        self.assertGreater(abs(changed.pc_column(b, self.i)[self.i]/c.pc_column(self.s, self.i)[self.i]-1), 1e-3)

    def test_direct_magnetic_response_vs_mass_derivative(self):
        direct = self.calc.magnetic_response(self.s, self.i)['delta_g']
        coarse = self.calc.mass_derivative(h=1e-3)
        fine = self.calc.mass_derivative(h=5e-4)
        self.assertLess(abs(fine/direct-1), 3e-7)
        self.assertLess(abs(fine/coarse-1), 5e-7)

    def test_coulomb_limit_and_nonrelativistic_weak_contact_limit(self):
        tiny = FermiDensity(0.1, 0.01)
        c = Calculation(NuclearInputs(1, 0, tiny, tiny, tiny, 'Small-nucleus analytic limiting test'), Grid(npoints=160, cavity_bohr=40))
        b = c.basis(); i = b.state(1)
        energy = np.sqrt(1-c.alpha**2)
        self.assertAlmostEqual(b.energies[i], energy, delta=2e-10)
        self.assertAlmostEqual(c.g_column(b, i)[i], 2/3*(1+2*energy), delta=2e-8)
        contact = c.cve*c.gf/(2*np.sqrt(2))*c.cvp*c.alpha**3/np.pi
        # Relativistic near-origin enhancement is O(alpha^2 log(lambda/R)).
        self.assertLess(abs(c.pc_column(b, i)[i]/contact-1), 2e-3)

    def test_full_apv_sum_is_phase_invariant_in_magnitude(self):
        from copy import deepcopy
        c, s = self.calc, self.s
        original = c.apv_1s_2s(s, self.p)
        p = deepcopy(self.p)
        p.g[:, ::2] *= -1; p.f[:, ::2] *= -1
        self.assertAlmostEqual(c.apv_1s_2s(s, p)['imag_e_bohr']/original['imag_e_bohr'], 1, delta=2e-13)
        s = deepcopy(s); f = s.state(2)
        s.g[:, f] *= -1; s.f[:, f] *= -1
        self.assertAlmostEqual(c.apv_1s_2s(s, p)['imag_e_bohr']/original['imag_e_bohr'], -1, delta=2e-13)

    def test_fortran_angular_label_regression(self):
        source = HERE.parent/'rad/src/weak_angular_labels.f90'
        program = """program check_labels
use weak_angular_labels
implicit none
if(p_two_j('3p1/2')/=1) stop 1
if(p_two_j('3p')/=1) stop 2
if(p_two_j('13p_{3/2}')/=3) stop 3
if(p_two_j('2P_{1/2}')/=1) stop 4
if(p_two_j('3s1/2')/=0) stop 5
if(p_two_j('2p33')/=0) stop 6
if(j_state_factor('3p1/2')/=1d0) stop 7
if(abs(j_state_factor('2p3/2')-sqrt(2d0))>1d-14) stop 8
end program
"""
        with tempfile.TemporaryDirectory() as work:
            from pathlib import Path
            work = Path(work)
            (work/'check.f90').write_text(program)
            subprocess.run(['gfortran', '-fcheck=all', source, work/'check.f90', '-o', work/'check'], cwd=work, check=True)
            subprocess.run([work/'check'], check=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
