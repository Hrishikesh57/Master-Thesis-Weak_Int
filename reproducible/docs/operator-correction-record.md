# Spherical operator correction and source-built baseline

22 September 2026. Scope: the one-electron spherical problem before adding
nuclear deformation. An unfinished correction was already present when this
implementation pass began; that work was reviewed, completed and validated.
Existing unrelated edits, thesis/slides, historical output tables and binaries
were preserved.

## Outcome

The new `workspace/programs/weak_spherical` route is runnable from source and
does not link the missing-source weak wrapper. It calculates a parity-even weak
energy shift, a parity-odd matrix element, a direct magnetic response, and a
mass-derivative cross-check. It also provides an explicitly labelled Dirac-only
1s-2s PV transition diagnostic. This is a numerical baseline, not a precision
paper result or a claim that the full legacy application has been validated.

## Physical and numerical changes

1. **Separate PC and PV operators.** The common nuclear source U contains
   `G_F/(2*sqrt(2))` and the weak nuclear density, but no electron vector factor.
   The PC operator multiplies U by Cve; the PV operator is `-gamma5*U` in the
   requested slide convention. The original hard-sphere potential now exposes
   `source_eval` separately from its existing PC `eval` method.
2. **Correct radial structure and phases.** For spinors
   `(G Omega_kappa, i F Omega_-kappa)/r`, the PV matrix element is
   `i*integral U*(F_a*G_b-G_a*F_b)dr`, nonzero only for opposite kappa and equal m.
   The original sampled driver's two weak matrix elements were put in the
   correct bra/ket order. The new route keeps complex values through both time
   orderings and takes a magnitude only after summing. Four-component angular
   integration independently verifies both this weak kernel and the explicit
   s1/2-p1/2 dipole coefficient.
3. **Repair state labels.** A small Fortran label module examines the angular
   suffix after p. `3p1/2` is no longer assigned the p3/2 factor. Existing bare
   `np` shorthand still means np1/2. Unsupported angular suffixes fail explicitly.
4. **Separate nuclear inputs.** Electric charge, point-proton and point-neutron
   Fermi profiles must each be supplied. Every profile integrates to one; Z and
   N appear only in the nuclear weak source. Radius input is unambiguous: c/a
   in fm, or conversion from rms through the actual Fermi moment equation.
   This permits realistic inputs later without silently treating a charge
   radius as a point-proton radius.
5. **Use direct interval quadrature.** The DKB wavefunctions are evaluated at
   the quadrature points used on each spline interval; no exported radial table
   is interpolated in the new route. Repeated origin knots stay zero and ordered;
   Gauss nodes do not touch the singular endpoint. The legacy change that moved
   just the first repeated zero knot was removed in the unfinished patch and
   its replacement is now exercised by the source-built calculation.
6. **Correct the Fermi Coulomb origin limit.** The potential at zero is its
   finite negative electrostatic limit, not zero. Independent integration of
   the spherical Poisson formula checks the origin and five other radii.
7. **Hold physical inputs fixed under a mass derivative.** Nuclear lengths and
   cavity lengths scale with mu in electron-mass coordinates; the dimensionless
   Fermi constant scales with mu squared. Differentiate `mu*delta_epsilon` with
   the positive `kappa^2/[j(j+1)]` coefficient. A separate Zeeman response includes
   the full positive- and negative-energy DKB spectrum.
8. **Make build provenance explicit.** The original `rad_kernel_f2py.f90` is
   still absent. A small new wrapper replaces the required spherical functions.
   Every project and SLATEC object is compiled in a fresh directory; system
   LAPACK/BLAS are linked explicitly. Import is checked before publishing the
   extension, and source/binary fingerprints are checked by the Python loader.
   Legacy loaders require explicit historical-use opt-in; their build script
   refuses to link old objects. This does not attempt to reconstruct the entire
   missing wrapper's historical API.
9. **Keep remaining diagnostics honest.** The sampled amplitude driver retains
   legacy overall angular/optical factors and imported radial data. Its output
   is labelled unvalidated, its asymmetry is NaN, and corrected files have new
   names. New physical results must use the source-built route. The original
   unit-conversion correction `G_F[MeV^-2]=10^-6 G_F[GeV^-2]` was verified; the
   active relativistic coupling was not changed. The already-present legacy
   mass-derivative cache fix includes the full nuclear input in the cache key.

## Validation and model readings

`unit-tests.txt`: **10 tests passed**. These cover normalization/orthogonality,
weak-charge conservation, exact Fermi rms inversion, electrostatic origin and
radial values, weak Hermiticity/selection rules, coupling separation, explicit
four-spinor/angular evaluation, independent nuclear interfaces, magnetic mass
derivatives, analytic point-Coulomb energy/g and nonrelativistic contact limits,
phase invariance of the summed PV amplitude, and compiled angular-label tests.

`validation.json`: **58 checks passed** for Pb-208 and U-238, with 16 grid/model
runs. Breakpoint counts 100, 160, 240 and 320 were examined; the 160-point
baseline was also checked with spline orders 7/9/11, 20/40 Gauss points per
interval and cavities 2.5/5 reference bohr. The 100-point result is a retained
coarse diagnostic; pass/fail comparisons use 160 points and above against 320.

At 160 points, order 9, 20 Gauss points and 2.5 bohr:

| Explicit common-profile model | PC energy shift / (m_ref c^2) | Imaginary 2p1/2–1s PV matrix element / (m_ref c^2) | Direct weak delta-g |
|---|---:|---:|---:|
| Pb-208, rms 5.5012 fm | -3.8884884e-12 | -6.3330329e-12 | -1.3328870e-11 |
| U-238, rms 5.8571 fm | -9.5998305e-12 | -1.9136949e-11 | -3.1233279e-11 |

For these numerical benchmarks, charge=point-proton=neutron is an **explicit
common-profile approximation**, with a=2.3/(4 ln 3) fm. The profiles' c values
are obtained by solving for the specified rms. The example JSON used by run.py
instead directly specifies c=6.624 fm and a=0.523 fm; its result is a different
documented input case and must not be compared as if the inputs were identical.

The smallest mass-difference step, h=0.00025, agrees with the direct magnetic
response to relative 7.35e-9 (Pb) and 5.03e-9 (U). Across the tested refined grids,
the maximum relative change is below 7e-12 for the PC shift, weak matrix element
and direct delta-g, and below 2.8e-10 for the Dirac-only transition diagnostic.
These are internal numerical comparisons, **not total physical uncertainties**.
The wider pass tolerances are recorded explicitly in validation.json.

Additional checks: original sampled driver compiled against fresh modules;
`practise.f90` passed a compiler syntax check using installed project module
interfaces; modified Python sources parsed; standard legacy main/sidecar loaders
refused unapproved reuse as expected; legacy build returned a failure with an
actionable replacement path. Both legacy .so hashes still match the saved
before/sha256.json record. See auxiliary_checks.json for the recorded checks.

## Files changed or completed

Paths are relative to the repository root.

- `workspace/codes/bspqed/src/bspqed_potentials.f90`: split source/PC evaluation,
  corrected Fermi origin (completed existing patch; tested here).
- `workspace/codes/bspqed/src/bspqed_physconsts.f90`: MeV^-2 conversion
  (existing patch, active relativistic coupling unchanged).
- `workspace/programs/rad/src/something.f90`: PV radial kernel and ordering,
  parser use, denominator guards and diagnostic output labels.
- `workspace/programs/rad/src/weak_angular_labels.f90`: new tested angular parser.
- `workspace/programs/rad/src/practise.f90`: valid origin knots (existing patch).
- `workspace/programs/rad/Makefile`: angular-module dependency for sampled driver.
- `workspace/programs/weak_spherical/{kernel.f90,build.py,spherical.py,run.py,
  validate.py,test_spherical.py,example_pb208.json,README.md,.gitignore}`:
  completed/new standalone calculation, documented inputs and validation.
  `build_manifest.json` records the built sources. Generated binaries are ignored.
- `workspace_f2py/programs/rad_f2py/{build_f2py.sh,README.md,
  python/example_rad_kernel.py}`: historical build/usage safeguards and routing.
- `workspace_f2py/python_api/workspace_f2py_api/{_common.py,direct_modules.py}`,
  `workspace_f2py/{First_order.py,Rcavity_sensitivity.py,Weak_source_diagnostic.py}`:
  explicit legacy-binary opt-in, with a pointer to the corrected route.
- `workspace_f2py/Mass_derivative_direct_dkb.py`: complete cache identity
  (existing patch, preserved).
- `research/operator_cleanup_2026-09-22/`: this report, logs, validation and
  example output records. Existing before/ snapshots were preserved.

The copied `workspace_f2py/programs/rad` source/executables and historical
result tables are **not** relabelled as corrected calculations. The thesis,
slides and the earlier report PDF have not been rewritten.

## Remaining physical work before a deformation prediction

The code uses the requested slide sign convention and the repository's fixed
tree-level couplings. The overall electroweak sign discrepancy in Greiner's
Example 6.10 remains documented; it does not justify changing relative signs.
No signed precision Standard Model claim is made here.

The 1s-2s PV amplitude is numerically stable for the finite-charge Dirac model,
but its near-degenerate denominators omit Lamb shifts and other physical level
corrections. An independent published APV benchmark using matched inputs and
level splittings is still needed. Many-electron effects, radiative corrections,
recoil and physical nuclear-density uncertainty are outside this baseline.

The radial infrastructure is ready for a controlled deformation feasibility
extension. A rank-2 weak density requires new angular/hyperfine couplings and
selection-rule checks; varying a spherical Fermi radius is not that extension.
Neither a deformation signal nor experimental feasibility has been established.
