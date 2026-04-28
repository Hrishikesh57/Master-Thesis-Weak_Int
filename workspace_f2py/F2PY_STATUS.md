## `workspace_f2py` Direct `f2py` Status

`workspace_f2py/` is a full copy of the protected `workspace/` tree.

Current direct Python wrapper coverage:

- `programs/gfactor_lithium_f2py/` provides the copied lithium g-factor suite through `gfactor_lithium_suite`.
- `programs/selfen_f2py/` provides the copied self-energy suite through `selfen_suite`.
- `programs/hartree_f2py/` provides copied Hartree wrappers through `hartree_suite`.
- `programs/rad_f2py/` provides copied sampled weak/dipole kernels through `rad_kernel_suite`.
- Wrapped copied workflows:
  - `g_dirac`
  - `g_delbrueck`
  - `g_propagator_loop`
  - `g_uehling`
  - `g_uehling_all`
  - `g_wichmann_kroll`
  - `zero`
  - `one`
  - `many`
  - `dhfcb`
  - `shf`
  - sampled first-order weak radial integral
  - sampled second-order weak energy kernel
  - sampled discrete Stark/weak amplitude kernel

Current Python-first copied-workspace coverage beyond direct `f2py`:

- `python_api/workspace_f2py_api/` provides Python helpers that build and run copied Fortran executables from:
  - `programs/gfactor/`
  - `programs/gfactor_boron/`
  - `programs/gfactor_lithium/`
  - `programs/selfen/`
  - `programs/hartree/`
  - `programs/rad/`
- `python_api/workspace_f2py_api/` also provides one unified direct-wrapper import surface for the currently exposed copied `f2py` modules:
  - `compute_lithium_dirac_tables(...)`
  - `compute_lithium_delbrueck_tables(...)`
  - `compute_lithium_propagator_loop_feynman_table(...)`
  - `compute_lithium_uehling_tables(...)`
  - `compute_lithium_uehling_all_tables(...)`
  - `compute_lithium_wichmann_kroll_tables(...)`
  - `compute_selfen_zero(...)`
  - `compute_selfen_one(...)`
  - `compute_selfen_many(...)`
  - `compute_hartree_dhfcb(...)`
  - `compute_hartree_shf(...)`
  - `compute_rad_weak_potential_summary(...)`
  - `compute_rad_first_order_weak_shift(...)`
  - `compute_rad_second_order_shift(...)`
  - `compute_rad_discrete_amplitude(...)`
- `build_python_backend.sh` builds the copied executable families plus the copied lithium `f2py` module from `workspace_f2py/` only.

Not yet turned into direct `f2py` APIs:

- the rest of the copied program trees such as `programs/gfactor/`, `programs/gfactor_boron/`, and the remaining unwrapped interactive/file-driven paths in `programs/rad/` and `programs/hartree/`
- the bulk of the copied library modules under `codes/`
- the copied third-party SLATEC source tree under `build/slatec/`

Reason for the staged approach:

- the copied lithium g-factor suite, self-energy suite, consistent Hartree drivers, and sampled `rad` kernels are coherent, high-value workflows with trusted Fortran kernels
- wrapping every Fortran file in the copied workspace blindly would mix unrelated physics workflows, third-party numerical libraries, and driver programs that are not yet interface-clean
- several copied drivers still depend on formatted stdin/stdout workflows, file generation, OpenMP loops, or nested callback structures that are better preserved first through Python executable bridges
- the current path lets you work from Python across the copied workspace immediately, while exposing the scientifically clean lithium, self-energy, consistent Hartree, and sampled `rad` suites as direct `f2py` modules

Known copied-workspace build caveats:

- `programs/hartree/dhfc` currently fails to compile in the copied tree because the existing source passes `type(sphere_potential)` into interfaces declared as `class(function_type)`. That is an upstream copied-source issue, not a Python wrapper issue.
- `programs/hartree_f2py/` intentionally wraps `dhfcb` and `shf` only. The copied `dhfc` path remains blocked until that potential/interface mismatch is resolved explicitly.
- `programs/rad_f2py/` intentionally wraps the sampled kernels below `something.f90`, not the full interactive amplitude/file-discovery workflow or the basis-spectrum derived-type path.
