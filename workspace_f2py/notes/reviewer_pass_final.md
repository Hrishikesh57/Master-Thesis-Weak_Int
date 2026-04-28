# Reviewer pass: final internal numerical pipeline

## 1. Overall verdict

**Promising but needs specific numerical/theory fixes.**

The project is internally coherent as a numerical development program. The final internal route is now a reproducible Pb/U benchmark pipeline:

```text
direct-DKB/intlw
Fermi Coulomb
separate proton/neutron Fermi weak density
H-like 1s
order = 9
npoints = 1000
rcavity = 2.5/alpha
nodes20/weights20
h = 5e-4 with h = 1e-3 stability check
external mu^2 weak-normalization scaling retained
```

The development sequence is logical: the sampled-rad route was first reproduced, then rejected after showing nuclear-region under-resolution; the direct-DKB/intlw route was introduced; p/n density generalizations were checked against common-density limits; Fermi Coulomb was independently validated; and the final route was converged in `npoints`.

However, I would not yet describe the result as a paper-ready physical prediction. The remaining risks are not primarily B-spline convergence risks. They are convention and interpretation risks: the mass-derivative sign/prefactor, the mapping of the internal weak source to standard electroweak notation, and the absence of an independent magnetic-response cross-check.

## 2. Strong points

- The project correctly quarantined the sampled-rad route once it was shown to under-resolve the nuclear weak region. That is an important numerical failure mode, and the report documents it rather than hiding it.
- The final weak-energy route uses direct DKB/intlw matrix elements, which is the right numerical direction for a short-range nuclear weak operator.
- The p/n hard-sphere and p/n Fermi implementations were validated against common-density limits. This is a useful guard against normalization, sign, and weighting mistakes.
- The Fermi Coulomb validation is strong for positive radii: the analytic `bspqed` Fermi Coulomb potential agrees with an independent normalized-density Coulomb integration at roundoff away from exactly `r = 0`.
- The final-route `npoints` convergence is convincing for the tested route. The maximum `npoints = 1000` relative difference versus `npoints = 1500` in `delta_g_W` is `7.5919628895751802e-08`, comfortably below the `1e-6` internal threshold.
- The final Pb/U neutron-skin table is smooth and physically interpretable. U is larger in magnitude than Pb, and increasing neutron skin makes the correction less negative.

## 3. Remaining blockers

**Must-fix before final paper claims: mass-derivative sign and prefactor.**

The working convention uses:

```text
D_W = [(1+h) Delta_epsilon_W(1+h) - (1-h) Delta_epsilon_W(1-h)] / (2h)
delta_g_W(1s) = (4/3) D_W
```

This is now internally frozen, but it still needs an independent audit against the literature identity and this codebase's bound-electron `g` convention. Earlier report sections mention a KLM-style formula with an opposite-looking sign. That does not automatically invalidate the numbers, because signs can move between definitions of `Delta_E_W`, `D_W`, electron charge, and `g`; but it is a paper-level blocker until reconciled explicitly.

**Must-fix before final paper claims: weak-source convention notation.**

The frozen internal source is:

```text
Cv_e = 1 - 4 sin^2(theta_W)
Cv_p = 1 - 4 sin^2(theta_W)
V_W^pc = G_F/(2 sqrt(2)) Cv_e [Cv_p Z rho_p - N rho_n]
```

This is executable and traceable, but it is not notation-safe for a paper. Standard electroweak conventions often use different signs and factors for `g_V^e`, `g_V^p`, neutron weak charge, or `Q_W`. The manuscript should not present `Cv_e` and `Cv_p` as if they were universally standard symbols. It should define internal coefficients or translate them into standard weak-charge notation.

**Must-fix before describing the Fermi Coulomb routine as generally production-safe: exact-origin branch.**

Section 37 shows that `fermi_potential_eval(0)` returns zero, while the physical central Fermi Coulomb potential is finite and negative. The DKB route appears safe if quadrature never evaluates exactly at `r = 0`, but this branch remains a sharp caveat. It can be tolerated for the current positive-radius DKB use; it should be fixed or explicitly quarantined before broader reuse.

## 4. Recommended additional checks before paper drafting

**Must-fix or must-document before drafting final results:**

- Reconcile the sign and prefactor in `delta_g_W = (4/3) D_W` against the original bound-electron `g`-factor identity, including the sign convention for the weak energy shift.
- Add a short convention note mapping the internal weak potential to standard neutral-current or weak-charge notation.
- Verify that the selected DKB eigenstate is robustly the H-like `1s` state across `mu = 1 +/- h`, `npoints`, and the Fermi Coulomb route. The current fixed positive-spectrum index is plausible, but a paper should state how the state was identified.

**Should-fix before final submission:**

- Run a small `rcavity` variation around `2.5/alpha` for Pb and U at the frozen route.
- Run a small order or quadrature variation, for example order around 9 or nodes20 versus a higher-order quadrature, for one representative Pb and U case.
- Perform one independent magnetic-response or controlled-potential check of the mass-derivative route. This is the cleanest way to test the sign and prefactor in the actual code conventions.
- Check whether moving the external `mu^2` scaling into the Fortran route gives the same finite-difference result, or document why the external placement is algebraically equivalent.

- Fermi density normalization and Fermi Coulomb validation are documented;
- finite-difference step checks are at roughly the `1e-7` level;
- final `npoints = 1000` values agree with `npoints = 1500` at better than `1e-6` relative by a wide margin.

The main missing numerical checks are not expected to overturn the table, but they are important for paper defensibility: `rcavity`, spline order, quadrature order, and eigenstate identification.

## 7. Nuclear-density/model-dependence assessment

The model progression is sensible:

The largest convention risks are:

- A factor-of-two ambiguity in the weak Hamiltonian prefactor.
- A sign ambiguity between electron weak coupling, proton/neutron weak charge notation, and the sign of the scalar potential used in the Dirac equation.
- A naming ambiguity from using `Cv_e` and `Cv_p` for quantities that may not match common `g_V^e` and `g_V^p` definitions.
- A global sign ambiguity in `delta_g_W` if the mass-derivative identity is imported with a different convention for `Delta_E_W`.

The current results are valid as internally defined numbers. They should not yet be presented as convention-independent physical predictions.
The final p/n Fermi model with common diffuseness is a reasonable internal benchmark model. The neutron-skin values `0.00`, `0.05`, `0.10`, and `0.20 fm` are reasonable sensitivity points, especially because they are presented as a controlled scan rather than a comprehensive nuclear-structure uncertainty model.

The limitation is that this is not yet a realistic nuclear-density error budget. A later paper version should either:

- keep the claim narrow and call this a neutron-skin sensitivity band; or
- add literature-motivated neutron radii and diffuseness uncertainties for Pb-208 and U-238.

## 8. Final-result interpretation

The final internal zero-skin values are:

```text
Pb-208: delta_g_W = -1.3328870635281932e-11
U-238:  delta_g_W = -3.1233277218115705e-11
```

The U/Pb hierarchy is plausible. The higher-Z ion has stronger near-nucleus wavefunction overlap and a larger weak source, so a larger magnitude for U is expected.

The neutron-skin trend is also plausible. Increasing the neutron rms radius moves part of the dominant neutron weak density outward, reducing overlap with the compact 1s electron density. In the internal convention, this makes `delta_g_W` less negative:

```text
Pb-208, skin 0.20 fm: relative magnitude change about -6.65e-3
U-238,  skin 0.20 fm: relative magnitude change about -7.77e-3
```

Nothing in the final table looks numerically erratic. The skin dependence is monotonic, small, and near-linear over the tested range. The caveat is that the absolute sign remains tied to the mass-derivative and weak-source convention audit.

## 9. Suggested paper narrative

Frame the paper as a compact H-like benchmark and numerical-methods result, not as an experimental observability claim.

Suggested main narrative:

1. Define the spin-independent parity-conserving weak potential and the exact internal convention used.
2. Derive or cite the bound-electron `g`-factor mass-derivative relation, with all signs and units reconciled.
3. Explain why short-range nuclear weak operators require direct matrix-element integration rather than sampled radial interpolation.
4. Present the direct-DKB/intlw route.
5. Present the p/n Fermi density model and Fermi Coulomb consistency.
6. Give compact validation: sampled route rejection, Coulomb validation, p/n recovery tests, finite-difference stability, and `npoints` convergence.
7. Present the Pb/U neutron-skin sensitivity table as the final internal benchmark.

Material that should move to appendices:

- detailed sampled-route history;
- long convergence tables;
- Fermi Coulomb analytic formula and numerical validation details;
- p/n common-density recovery tests;
- implementation-level wrapper details.

Material to omit from the main paper narrative:

- obsolete intermediate benchmark tables unless needed to explain why the sampled route was retired;
- broad experimental feasibility claims;
- many-element extrapolations not run with the final route.

## 10. Concrete next actions, prioritized

1. **Must-fix:** write a convention audit note reconciling `delta_g_W = (4/3) D_W` with the literature mass-derivative formula, including sign, units, and the definition of `Delta_E_W`.
2. **Must-fix:** translate the internal weak source into standard weak-charge notation or define clearly named internal coefficients that cannot be mistaken for standard `g_V` symbols.
3. **Must-fix:** verify the H-like `1s` state selection in the DKB spectrum across the final-route finite-difference calls.
4. **Should-fix:** run one independent magnetic-response or controlled-potential sign/prefactor check.
5. **Should-fix:** add small `rcavity`, order, and quadrature spot checks for one Pb and one U case.
6. **Should-fix:** either fix the exact-origin Fermi Coulomb branch or document that the routine is validated only for positive-radius quadrature use.
7. **Optional:** add realistic neutron-skin and diffuseness uncertainty inputs after the core convention checks are complete.

The project is in a strong internal state, but the next work should be theory/convention hardening, not more production tables.
6. Give compact validation: sampled route rejection, Coulomb validation, p/n recovery tests, finite-difference stability, and `npoints` convergence.
7. Present the Pb/U neutron-skin sensitivity table as the final internal benchmark.

Material that should move to appendices:

- detailed sampled-route history;
- long convergence tables;
- Fermi Coulomb analytic formula and numerical validation details;
- p/n common-density recovery tests;
- implementation-level wrapper details.

Material to omit from the main paper narrative:

- obsolete intermediate benchmark tables unless needed to explain why the sampled route was retired;
- broad experimental feasibility claims;
- many-element extrapolations not run with the final route.

## 10. Concrete next actions, prioritized

1. **Must-fix:** write a convention audit note reconciling `delta_g_W = (4/3) D_W` with the literature mass-derivative formula, including sign, units, and the definition of `Delta_E_W`.
2. **Must-fix:** translate the internal weak source into standard weak-charge notation or define clearly named internal coefficients that cannot be mistaken for standard `g_V` symbols.
3. **Must-fix:** verify the H-like `1s` state selection in the DKB spectrum across the final-route finite-difference calls.
4. **Should-fix:** run one independent magnetic-response or controlled-potential sign/prefactor check.
5. **Should-fix:** add small `rcavity`, order, and quadrature spot checks for one Pb and one U case.
6. **Should-fix:** either fix the exact-origin Fermi Coulomb branch or document that the routine is validated only for positive-radius quadrature use.
7. **Optional:** add realistic neutron-skin and diffuseness uncertainty inputs after the core convention checks are complete.

The project is in a strong internal state, but the next work should be theory/convention hardening, not more production tables.
