# Initial Physics Scope

## Observable

The subproject targets the specific-difference observable

\[
\delta_\xi g = g_{\mathrm{Li\text{-}like}} - \xi g_{\mathrm{H\text{-}like}},
\]

with `\xi` chosen to suppress the leading finite nuclear size contribution.

## First-Pass Physics Scope

- start with H-like `1s`
- start with Li-like `2s`
- prioritize Pb and U first
- compute absolute g-factor ingredients first
- then finite nuclear size sensitivity
- then determine `\xi`
- then build the specific-difference observable

## Numerical Scope

- reuse copied Dirac/B-spline infrastructure already present in `workspace_f2py/`
- reuse copied H-like and Li-like g-factor drivers before introducing any new weak-`g` numerics
- keep the first `\xi` definition at the Dirac+FNS level
- keep weak `g` work explicitly provisional until it is validated against local thesis or paper targets

## Out Of Scope For The First Pass

- boron-like ions
- full screened many-body uncertainty analysis
- publication-quality weak-`g` numbers for Li-like ions
- silent promotion from the spherical nuclear model to Fermi or neutron-skin variants

## Conventions To Preserve

- nuclear RMS radius input in fm
- existing `rms_fm -> rnucl` convention from the copied g-factor drivers
- existing `R_fm -> R_rel` convention from the copied weak-kernel path
- existing sign and normalization conventions from the reused copied source files
