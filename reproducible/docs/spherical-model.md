# Source-built spherical weak-interaction calculation

This is the corrected one-electron spherical baseline. It reuses the original
Dirac dual-kinetic-balance (DKB) solver, with a small source-complete wrapper.
It does not use `rad_kernel_suite`, the copied `.o` files, or project archives.
The missing original wrapper has **not** been recovered; this is an explicit
replacement for the spherical calculations, not a drop-in replacement for its
entire undocumented API.

## Run from the repository root

Requires Python with NumPy and SciPy, gfortran, gcc, binutils, and LAPACK/BLAS.

```sh
python3 workspace/programs/weak_spherical/build.py
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover -s workspace/programs/weak_spherical -p 'test_*.py' -v
OPENBLAS_NUM_THREADS=1 python3 workspace/programs/weak_spherical/run.py \
  --inputs workspace/programs/weak_spherical/example_pb208.json \
  --output /tmp/weak-pb208-example.json
```

Use a new output filename on each run: existing records are never overwritten.
The example assigns equal charge/proton/neutron profiles **as an explicit model
assumption**. It is not a measured nuclear-density fit. Supply all three profiles
and their provenance for a physical application; there is no silent fallback
from a missing point-proton density to the electric charge density.

To repeat the convergence study, with a new output filename:

```sh
OPENBLAS_NUM_THREADS=1 python3 workspace/programs/weak_spherical/validate.py \
  --output /tmp/weak-spherical-validation.json
```

The build uses a fresh temporary directory, resolves required SLATEC routines
from repository source, checks runtime symbols before installing the extension,
and records source/binary hashes, compiler flags, compiler version, NumPy and
Python versions and commands in `build_manifest.json`. LAPACK/BLAS are system
dependencies, not rebuilt here. The Python loader rejects changed source or
a binary that does not match this manifest. Numerical result records embed the
manifest and hashes of the Python calculation code.

## Physics and units

The working sign convention is the user's slide:

$$
 U_W=\frac{G_F}{2\sqrt2}
 (C_V^p Z\rho_p-N\rho_n),\qquad
 H_{\rm PC}=C_V^eU_W\,1_4,\qquad H_{\rm PV}=-U_W\gamma^5.
$$

Each nuclear profile obeys `4*pi*integral r^2*rho dr = 1`. Charge, point-proton
and point-neutron distributions are independent inputs. `c_fm` is the Fermi
half-density radius, **not** the rms radius. `FermiDensity.from_rms` solves the
actual Fermi moment equation rather than assuming the thin-surface formula is
exact. Density normalization includes the infinite tail; electronic integrals
extend to the finite cavity, whose convergence must be checked.

This preserves the slide's overall sign. It does not resolve the inconsistency
between Eqs. (2a), (6), (8) and (9) in Greiner and Mueller, *Gauge Theory of Weak
Interactions* (2009), Example 6.10, pp. 242-244. Signed Standard Model predictions
need that matching; a consistent overall reversal preserves first-order
magnitudes. Relative signs and interference are explicitly retained here.

Spinors are `(G Omega_kappa, i F Omega_-kappa)/r`, with real G/F and
`integral (G^2+F^2)dr=1`, and gamma5 has positive off-diagonal unit blocks.
The code evaluates

$$
 \langle a|H_{\rm PV}|b\rangle
 =i\delta_{\kappa_a,-\kappa_b}\delta_{m_a m_b}
 \int U_W(F_aG_b-G_aF_b)\,dr.
$$

There is no electron vector factor in this expression. Whole-state rephasing
always changes G and F together. The PC kernel uses `G_a G_b+F_a F_b`.

At mass ratio `mu=m_e/m_ref`, the solver uses coordinates
`x=mu*r_phys/lambda_ref` and energies in `mu*m_ref*c^2`. Every physical nuclear
length and the physical cavity remain fixed during differentiation; their
dimensionless values scale with mu. The dimensionless weak coupling is
`G_F*m_ref^2*mu^2`. Differentiate `mu*delta_epsilon(mu)`, not delta_epsilon alone:

$$
 \delta g=\frac{\kappa^2}{j(j+1)}
 \frac{d[\mu\delta\epsilon]}{d\mu}.
$$

An independent first-order magnetic response sums over the complete finite
same-kappa DKB spectrum, **including negative-energy states**. Its magnetic
matrix is `kappa/[j(j+1)] * integral r*(G_a F_b+F_a G_b)dr`.

The optional `1s -> 2s` APV output is the length-gauge `d_z` amplitude for
`m=+1/2`, with `d_z=-r*cos(theta)` and the Condon-Shortley angular coefficient
`+1/3` multiplying each s1/2-p1/2 dipole radial integral. Both time orderings
and the full finite p1/2 spectrum are included; matrix elements are summed
before taking the magnitude. It uses pure finite-charge Dirac energies:
**it is a diagnostic, not a precision APV prediction**. Near-degenerate
physical level splittings, Lamb shifts and possibly widths must be supplied
before experimental interpretation. No asymmetry is inferred without geometry.

The active constants and tree-level weak mixing angle are inherited from the
project; this change does not update them to a modern radiative scheme. No
recoil, many-electron correlations, nuclear polarization, finite single-nucleon
weak form factors or nuclear deformation are included.

## Legacy diagnostic and omitted extensions

The build also compiles the inherited `sampled_radial_diagnostic` against the same new objects. Its interactive sampled-table workflow and external tables are not part of this curated release. It is not used by `run.py`, the tests, or the convergence study. Its source remains solely because the inherited build checks it. No historical `.so`, `.o`, executable, or unpublished source PDF is required for the released calculation.

See [the preserved September correction record](operator-correction-record.md) for the changes that distinguish this baseline from the original thesis implementation. Historical filenames in that record refer to the private research workspace; the current reproduction commands are in the repository README.

The later `apv_reference.py`, photon-current, quadrupole, scattering, and molecular extensions are outside this release. This extraction provides the spherical baseline only. A published-value APV comparison, radiative matching, and a complete uncertainty budget have not been established by running its internal tests.
