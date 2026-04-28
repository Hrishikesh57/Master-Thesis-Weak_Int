# `gfactor_specific_difference`

Python-first, paper-oriented subproject for the specific-difference observable

\[
\delta_\xi g = g_{\mathrm{Li\text{-}like}} - \xi g_{\mathrm{H\text{-}like}} ,
\]

with the initial goal of suppressing the leading finite nuclear size (FNS) contribution and then testing whether a residual weak-interaction contribution survives that cancellation.

## Current Scope

- H-like `1s` and Li-like `2s`
- Pb and U first
- absolute g-factor ingredients first
- finite nuclear size sensitivity next
- initial `\xi` from Dirac+FNS only
- weak contribution treated as a provisional follow-on layer

## Current Non-Goals

- boron-like ions
- a full uncertainty budget for publication use
- a full screened-correction stack inside the first `\xi`
- silent reinterpretation of legacy lithium-like driver outputs beyond what is documented in the copied source tree

## Reuse Strategy

This subproject reuses the copied Fortran backends already present in `workspace_f2py/` instead of reimplementing B-spline or g-factor numerics:

- H-like `gfactor` drivers through `workspace_f2py/python_api/workspace_f2py_api/executables.py`
- Li-like `gfactor_lithium` suite through `workspace_f2py/programs/gfactor_lithium_f2py/`
- weak-potential normalization and sampled weak kernels through `workspace_f2py/programs/rad_f2py/`

The original `workspace/` tree remains untouched and is not part of this subproject.

## Layout

- `docs/`
  scope, roadmap, and code-reuse map
- `src/gfactor_specific_difference/`
  Python package for cases, backend interfaces, FNS response, `\xi`, and reporting
- `configs/`
  default Pb/U H-like and Li-like cases plus specific-difference run configs
- `scripts/`
  driver entrypoints for the first pass
- `tests/`
  smoke and convention checks
- `data/`
  project-local reference notes and future tabulated inputs
- `output/`
  generated summaries, tables, and plots

## First-Pass Status

- the project skeleton, docs, configs, and Python package are in place
- H-like Dirac baseline reuse is wired through the copied executable bridge
- Li-like copied-driver reuse is wired through the copied `f2py` suite
- finite-size response and initial `\xi` construction are implemented numerically
- weak `g` contribution remains explicitly provisional because the repository does not yet contain a validated legacy weak-`g` driver to wrap directly

## Example Entry Points

Placeholder-safe script mode:

```bash
python3 workspace_f2py/programs/gfactor_specific_difference/scripts/run_hlike_u.py --placeholder --json
python3 workspace_f2py/programs/gfactor_specific_difference/scripts/run_specific_difference.py --placeholder --json
```

Project tests:

```bash
python3 -m unittest discover -s workspace_f2py/programs/gfactor_specific_difference/tests
```

## Important Physics Notes

- The H-like baseline maps cleanly onto the copied `workspace_f2py/programs/gfactor/` driver family.
- The Li-like side reuses the copied `workspace_f2py/programs/gfactor_lithium/` driver family as-is.
  The final paper-level interpretation of the copied lithium `g_dirac` workflow inside the specific-difference stack still needs to be pinned down against the original thesis and paper context before publication claims are made.
- The default finite-size convention in this first pass matches the existing spherical finite-nucleus convention used by the copied g-factor drivers for reproducibility.
- Promotion to Fermi/skin studies is planned but not enabled by default.
