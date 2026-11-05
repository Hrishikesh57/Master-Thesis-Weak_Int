"""Spherical one-electron weak operators, with explicit physical nuclear inputs.

Working (slide) convention: H_PC=Cve U, H_PV=-gamma5 U. This is not a
resolution of the overall electroweak sign discrepancy in Greiner Example 6.10.
Spinor=(G Omega_kappa, i F Omega_-kappa)/r, integral (G^2+F^2)dr=1.
All matrix elements below use dimensionless electron-mass units at mass mu.
"""
from dataclasses import asdict, dataclass
from functools import cached_property
import hashlib
import importlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit

HERE = Path(__file__).resolve().parent


def load_kernel():
    """Refuse a binary after any recorded source or binary has changed."""
    manifest_path = HERE / 'build_manifest.json'
    if not manifest_path.exists():
        raise RuntimeError('Run python3 build.py in weak_spherical first.')
    manifest = json.loads(manifest_path.read_text())
    root = HERE.parents[2]
    for name, expected in manifest['sources_sha256'].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'Source changed: {name}; rebuild with build.py.')
    kernel = importlib.import_module('weak_spherical_kernel')
    if hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest() != manifest['extension_sha256']:
        raise RuntimeError('Kernel fingerprint differs from build manifest; rebuild.')
    return kernel


@dataclass(frozen=True)
class FermiDensity:
    """Point-centre or charge profile; c and a are physical lengths, in fm.

    The profile integrates to ONE, not Z or N. Its physical role is supplied
    by the NuclearInputs field, and must not be inferred from its radius.
    """
    c_fm: float
    a_fm: float

    def __post_init__(self):
        if not np.isfinite([self.c_fm, self.a_fm]).all() or min(self.c_fm, self.a_fm) <= 0:
            raise ValueError('Fermi c and a must be finite positive physical lengths.')

    def moment(self, power):
        # Integrate both regions separately so an arbitrarily narrow surface
        # is not skipped by adaptive quadrature. The infinite tail is included.
        c, a = self.c_fm, self.a_fm
        inner = quad(lambda x: x**power * expit((c-x)/a), 0, c,
                     epsabs=1e-13, epsrel=2e-12)[0]
        outer = quad(lambda t: a*(c+a*t)**power*expit(-t), 0, np.inf,
                     epsabs=1e-13, epsrel=2e-12)[0]
        return inner + outer

    @cached_property
    def normalization(self):
        return 4*np.pi*self.moment(2)

    @cached_property
    def rms_fm(self):
        return np.sqrt(self.moment(4)/self.moment(2))

    @classmethod
    def from_rms(cls, rms_fm, a_fm):
        if not np.isfinite([rms_fm, a_fm]).all() or min(rms_fm, a_fm) <= 0:
            raise ValueError('rms and diffuseness must be finite positive lengths.')
        lower = max(1e-10, a_fm*1e-9)
        if cls(lower, a_fm).rms_fm >= rms_fm:
            raise ValueError('Requested rms is too small for this diffuseness and c>0.')
        # Solve the actual normalized Fermi moments, not the thin-surface approximation.
        return cls(brentq(lambda c: cls(c, a_fm).rms_fm-rms_fm,
                          lower, 2*rms_fm, xtol=1e-12), a_fm)

    def physical_density(self, radius_fm):
        return expit((self.c_fm-np.asarray(radius_fm))/self.a_fm)/self.normalization

    def density(self, r, lambda_fm, mu=1.0):
        # x = mu*r_phys/lambda_ref. rho_x = (lambda_ref/mu)^3 rho_phys.
        length = lambda_fm/mu
        return length**3*self.physical_density(np.asarray(r)*length)


@dataclass(frozen=True)
class NuclearInputs:
    z: int
    n: int
    charge: FermiDensity
    proton: FermiDensity
    neutron: FermiDensity
    provenance: str

    def __post_init__(self):
        if not isinstance(self.z, int) or not isinstance(self.n, int) or self.z <= 0 or self.n < 0:
            raise ValueError('Z must be a positive integer and N a nonnegative integer.')
        if not self.provenance.strip():
            raise ValueError('State the source or model assumption for the nuclear inputs.')


@dataclass(frozen=True)
class Grid:
    npoints: int = 200
    order: int = 9
    nq: int = 20
    cavity_bohr: float = 2.5

    def __post_init__(self):
        if self.npoints < 40 or self.order < 4 or self.order >= self.npoints or self.nq not in (20, 40):
            raise ValueError('Require npoints>=40, 4<=order<npoints, and nq=20 or 40.')
        if not np.isfinite(self.cavity_bohr) or self.cavity_bohr <= 0:
            raise ValueError('Cavity radius must be finite and positive.')


@dataclass
class Basis:
    kappa: int
    mu: float
    energies: np.ndarray
    r: np.ndarray
    w: np.ndarray
    g: np.ndarray
    f: np.ndarray

    def state(self, n):
        if not isinstance(n, int) or n < abs(self.kappa) + (self.kappa > 0):
            raise ValueError('Principal quantum number is incompatible with kappa.')
        bound = np.flatnonzero((self.energies > 0) & (self.energies < 1))
        ordinal = n - abs(self.kappa) - (self.kappa > 0)
        if ordinal >= len(bound):
            raise ValueError('Requested bound state is absent; increase the cavity/basis.')
        return int(bound[ordinal])

    def norm(self, i):
        return float(self.w @ (self.g[:, i]**2+self.f[:, i]**2))


class Calculation:
    def __init__(self, nucleus, grid=Grid()):
        self.nucleus, self.grid = nucleus, grid
        self.kernel = load_kernel()
        self.alpha, self.lambda_fm, self.gf, self.cve, self.cvp = self.kernel.spherical_constants()

    def basis(self, kappa=-1, mu=1.0):
        if not isinstance(kappa, int) or kappa == 0 or not np.isfinite(mu) or mu <= 0:
            raise ValueError('Nonzero integer kappa and finite positive mass ratio required.')
        nu, gr = self.nucleus, self.grid
        e, r, w, g, f, ierr = self.kernel.solve_spherical_basis(
            nu.z, nu.charge.c_fm, nu.charge.a_fm, gr.npoints, gr.order, gr.nq,
            gr.cavity_bohr, mu, kappa)
        if ierr or not all(np.isfinite(a).all() for a in (e, r, w, g, f)):
            raise RuntimeError(f'DKB solve failed (ierr={ierr}).')
        return Basis(kappa, mu, e, r, w, g, f)

    def source(self, basis):
        nu = self.nucleus
        weak_density = (self.cvp*nu.z*nu.proton.density(basis.r, self.lambda_fm, basis.mu)
                        - nu.n*nu.neutron.density(basis.r, self.lambda_fm, basis.mu))
        # G_F*(mu*m_ref)^2: physical G_F and all nuclear lengths stay fixed.
        return self.gf*basis.mu**2/(2*np.sqrt(2))*weak_density

    @staticmethod
    def same_grid(a, b):
        if a.mu != b.mu or not np.array_equal(a.r, b.r) or not np.array_equal(a.w, b.w):
            raise ValueError('Matrix elements require the same radial quadrature and mass.')

    def pc_column(self, basis, i):
        weighted = basis.w*self.cve*self.source(basis)
        return basis.g.T @ (weighted*basis.g[:, i]) + basis.f.T @ (weighted*basis.f[:, i])

    def pv_column(self, bra, ket, i):
        """Full complex <bra states,m|H_PV|ket i,m>, for equal m."""
        self.same_grid(bra, ket)
        if bra.kappa != -ket.kappa:
            return np.zeros(len(bra.energies), dtype=complex)
        weighted = ket.w*self.source(ket)
        return 1j*(bra.f.T @ (weighted*ket.g[:, i]) - bra.g.T @ (weighted*ket.f[:, i]))

    @staticmethod
    def g_column(basis, i):
        j = abs(basis.kappa)-0.5
        weighted = basis.w*basis.r*basis.kappa/(j*(j+1))
        return basis.g.T @ (weighted*basis.f[:, i]) + basis.f.T @ (weighted*basis.g[:, i])

    def magnetic_response(self, basis, i):
        """First-order Zeeman response; includes positive AND negative energy states."""
        v = self.pc_column(basis, i)
        g = self.g_column(basis, i)
        keep = np.arange(len(v)) != i
        den = basis.energies[i]-basis.energies[keep]
        if np.min(abs(den)) < 1e-12:
            raise ValueError('Unresolved same-parity degeneracy in magnetic response.')
        terms = 2*g[keep]*v[keep]/den
        negative = basis.energies[keep] < 0
        return {'delta_g': float(terms.sum()), 'negative_energy_part': float(terms[negative].sum()),
                'positive_energy_part': float(terms[~negative].sum()), 'g_unperturbed': float(g[i])}

    def mass_derivative(self, n=1, kappa=-1, h=5e-4):
        if not 0 < h < 0.05:
            raise ValueError('Require 0<h<0.05; check step convergence separately.')
        energy = []
        for mu in (1-h, 1+h):
            basis = self.basis(kappa, mu)
            i = basis.state(n)
            # E_phys/m_ref = mu*epsilon; do not differentiate epsilon alone.
            energy.append(mu*self.pc_column(basis, i)[i])
        j = abs(kappa)-0.5
        return float(kappa**2/(j*(j+1))*(energy[1]-energy[0])/(2*h))

    def apv_1s_2s(self, s, p):
        """Diagnostic d_z amplitude, m=+1/2, pure finite-charge Dirac denominators.

        d_z=-r cos(theta); Condon-Shortley spinor harmonics give +1/3
        times the s1/2-p1/2 radial integral. No Stark geometry is assumed.
        Includes both time orderings and the complete finite p1/2 spectrum.
        """
        self.same_grid(s, p)
        if s.kappa != -1 or p.kappa != 1:
            raise ValueError('This explicit angular coefficient is only for s1/2-p1/2.')
        i, f = s.state(1), s.state(2)
        def dipole(n):
            weight = s.w*s.r/3
            return p.g.T @ (weight*s.g[:, n]) + p.f.T @ (weight*s.f[:, n])
        wi, wf = self.pv_column(p, s, i), self.pv_column(p, s, f)
        di, df = dipole(i), dipole(f)
        den_i, den_f = s.energies[i]-p.energies, s.energies[f]-p.energies
        if min(np.min(abs(den_i)), np.min(abs(den_f))) < 1e-12:
            raise ValueError('APV denominator unresolved; supply physical level energies in a future route.')
        first, second = df*wi/den_i, wf.conj()*di/den_f
        # Convert length x (Compton units at mu) into reference bohr.
        conversion = self.alpha/s.mu
        amplitude = conversion*(first+second).sum()
        return {'imag_e_bohr': float(amplitude.imag), 'abs_e_bohr': float(abs(amplitude)),
                'first_ordering_imag_e_bohr': float(conversion*first.sum().imag),
                'second_ordering_imag_e_bohr': float(conversion*second.sum().imag),
                'smallest_denominator_mc2': float(min(np.min(abs(den_i)), np.min(abs(den_f)))),
                'status': 'Dirac-only diagnostic; Lamb shifts, linewidths and experimental geometry absent'}

    def inputs(self):
        return {'nucleus': asdict(self.nucleus), 'grid': asdict(self.grid),
                'convention': 'slide: H_PC=+Cve*U, H_PV=-gamma5*U; tree-level couplings',
                'spinor': '(G Omega_kappa, i F Omega_-kappa)/r',
                'constants': {'alpha': self.alpha, 'lambda_fm': self.lambda_fm,
                              'G_F_m_ref_squared': self.gf, 'Cve': self.cve, 'Cvp': self.cvp}}
