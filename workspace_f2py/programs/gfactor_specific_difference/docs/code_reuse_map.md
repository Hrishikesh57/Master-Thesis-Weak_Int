# Code Reuse Map

## Reuse Directly

- `workspace_f2py/codes/bspqed/src/bspqed_dkb.f90`
  - `solve`
  - `radialgf`
- `workspace_f2py/codes/bspqed/src/bspqed_discretization.f90`
  - `nodwgt`
  - `storewf`
  - `storefunc`
  - `matelems`
- `workspace_f2py/codes/bspqed/src/bspqed_physconsts.f90`
  - `alpha`
  - `electron_reduced_compton`
  - weak constants
- `workspace_f2py/codes/bspqed/src/bspqed_potentials.f90`
  - `sphere_potential`
  - `fermi_potential`
  - `weak_potential`
- `workspace_f2py/python_api/workspace_f2py_api/executables.py`
  - H-like executable bridge
- `workspace_f2py/python_api/workspace_f2py_api/direct_wrappers.py`
  - Li-like direct wrappers
- `workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90`
  - active Li-like wrapped backend
- `workspace_f2py/programs/rad_f2py/src/rad_kernel_f2py.f90`
  - weak-potential normalization
  - sampled weak kernels

## Wrap

- H-like `gfactor` drivers, because they are still executable-style rather than direct reusable Python calls
- any new weak-`g` matrix-element routines needed for M1 and M2
- any direct nucleus-model evaluators needed for controlled FNS sweeps

## Copy Carefully

- wrapper design pattern from `workspace_f2py/programs/gfactor_lithium_f2py/`
- validation-script pattern from `workspace_f2py/programs/*_f2py/validation/`
- only if required, minimal H-like wrapper files under this subproject or a sibling wrapper path

## Avoid

- editing the original `workspace/`
- depending on interactive `stdin` workflows as the long-term scientific API
- using `workspace_f2py/programs/hartree/dhfc`
- bulk-copying large library trees when a narrow wrapper or executable bridge suffices
- using `workspace_f2py/programs/rad/src/something.f90` itself as the long-term API

## Reuse Caveats

- the copied H-like `gfactor` family maps cleanly onto H-like baseline ingredients
- the copied Li-like driver family is reused as-is, but the final paper-level interpretation of its `g_dirac` workflow inside the specific-difference stack remains a documented ambiguity to resolve
- the weak path currently provides normalization and sampled matrix-element infrastructure, not a validated legacy weak-`g` driver
