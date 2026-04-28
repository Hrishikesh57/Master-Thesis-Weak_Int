# Weak g-Factor Specific-Difference Plan

## Summary

Build a Python-first, paper-oriented subproject under `workspace_f2py/programs/gfactor_specific_difference/` that reuses copied Fortran backends already present in `workspace_f2py/`, keeps the original `workspace/` untouched, and stages the physics in this order:

1. reproduce absolute H-like and Li-like baseline `g`-factor ingredients,
2. extract finite nuclear size sensitivity numerically in a consistent Dirac/B-spline framework,
3. determine an initial Dirac+FNS cancellation coefficient `\xi`,
4. add a weak-interaction contribution on top of the same framework,
5. test whether a residual weak contribution survives in `\delta_\xi g`.

## Existing Code Findings

- H-like baseline already exists in `workspace_f2py/programs/gfactor/`.
- Li-like baseline already exists in `workspace_f2py/programs/gfactor_lithium/` and is already wrapped in `workspace_f2py/programs/gfactor_lithium_f2py/`.
- Nuclear and weak models already exist in `workspace_f2py/codes/bspqed/src/`.
- B-spline and DKB infrastructure already exist in `workspace_f2py/codes/bspqed/src/`.
- Weak radial kernels exist in `workspace_f2py/programs/rad_f2py/`, but there is no validated legacy weak-`g` driver to wrap directly.

## Implementation Defaults

- Subproject root: `workspace_f2py/programs/gfactor_specific_difference/`
- Primary style: Python-first
- Config format: TOML
- First systems: H-like Pb `1s`, H-like U `1s`, Li-like Pb `2s`, Li-like U `2s`
- Initial `\xi`: Dirac+FNS only
- Initial charge model: same spherical finite-size convention used by the existing copied g-factor drivers
- First weak target: H-like only

## First-Pass Deliverables

- clean subproject directory structure
- `README.md`
- `AGENTS.md`
- `docs/scope.md`
- `docs/plan.md`
- `docs/code_reuse_map.md`
- Python package under `src/gfactor_specific_difference/`
- default Pb/U config files
- stub and smoke driver scripts
- unit and backend-smoke tests

## Immediate Physics Ambiguities To Track

- exact weak `g`-factor formalism for H-like ions
- whether Li-like weak contribution requires a new screened weak kernel
- when to promote from the spherical FNS model to a Fermi or neutron-skin study
- when the larger correction stack should enter `\xi` beyond the initial Dirac+FNS definition
