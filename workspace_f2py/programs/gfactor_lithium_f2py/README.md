# Lithium-Like g-Factor `f2py` Bridge

This directory is an isolated Python-interface path for the lithium-like bound-electron g-factor code.
The protected original workflow remains in `workspace/programs/gfactor_lithium`.
All `f2py` work for this bridge now lives inside the full copied workspace tree `workspace_f2py/`, so the original `workspace/` is no longer part of the wrapper build or validation path.

## What Was Copied

- `src/g_dirac_reference.f90`
  Exact copy of the original `workspace/programs/gfactor_lithium/src/g_dirac.f90` driver for traceability.
- `src/copied/bspqed_function.f90`
- `src/copied/bspqed_physconsts.f90`
- `src/copied/bspqed_quadrature.f90`
- `src/copied/bspqed_dkb.f90`
- `src/copied/bspqed_integrals.f90`
- `src/copied/dbspvd.f`
- `src/copied/dbspvn.f`
- `src/copied/dbvalu.f`
- `src/copied/dintrv.f`
- `src/copied/bspgf_lithium_reference.f90`
  Minimal upstream source snapshot set needed to build the first wrapped path without touching the protected original libraries.

## What Was Added

- `src/g_dirac_f2py.f90`
  A wrapper-oriented Fortran module that preserves the original `g_dirac` physics and calls the same `fbza` and `cbza` kernels from `bspgf_lithium`.
- `src/g_dirac_f2py_support.f90`
  Local support module holding only the sphere/linear potentials and knot-generation helpers needed by the first wrapper.
- `src/bspgf_lithium_dirac_kernels.f90`
  Reduced copy of the lithium kernel module containing only the trusted `fbza` and `cbza` paths plus the local helper routines they require.
- `src/xermsg_stub.f90`
  Minimal SLATEC error hook used only inside the copied wrapper build.
- `build_f2py.sh`
  Builds the Python extension module in place with `numpy.f2py`.
- `python/example_g_dirac.py`
  Minimal Python example for running the wrapped Dirac path.
- `validation/build_reference_g_dirac.sh`
  Builds the untouched `g_dirac` driver from the copied workspace into `validation/bin/` for comparison.
- `validation/validate_g_dirac.py`
  Compares the wrapped results against the original executable on the same input.

## Interface Scope

The first wrapped path is the lithium-like Dirac contribution (`g_dirac` logic).
This is the narrowest useful bridge because:

- the physical inputs are already scalar (`Z`, nuclear RMS radius, spline order),
- the trusted numerical core lives in existing library kernels (`fbza`, `cbza`),
- the original driver-side loop is thin and deterministic,
- no heavy file I/O or interactive control flow is required beyond the original stdin reads.

## Why The Build Is Local

The copied `workspace_f2py/include/lib/*.a` archives preserve the original executable path, but they were not compiled with `-fPIC`, so they cannot be linked into a Python extension shared object.

To preserve the original workflow exactly as a failsafe, the `f2py` build recompiles only the minimal copied source set inside this directory and does not modify the protected original archives or programs under `workspace/`.

## Build Commands

Build the Python module:

```bash
cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
./workspace_f2py/programs/gfactor_lithium_f2py/build_f2py.sh
```

Build the untouched original reference executable used for validation:

```bash
cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
./workspace_f2py/programs/gfactor_lithium_f2py/validation/build_reference_g_dirac.sh
```

## Run Commands

Run the Python example:

```bash
cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
python3 workspace_f2py/programs/gfactor_lithium_f2py/python/example_g_dirac.py
```

Run the validation comparison:

```bash
cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
python3 workspace_f2py/programs/gfactor_lithium_f2py/validation/validate_g_dirac.py
```

## Notes On Refactoring

- The protected original driver reads `Z` and `rms` from stdin and writes the convergence table to stdout.
- The wrapper version keeps the same conversions and kernel calls, but returns arrays for:
  `npoints`, `irr`, `red`, and `tot` in each gauge.
- No equations, scaling factors, normalization conventions, units, or B-spline settings were changed for valid inputs.
- Input checks were added only to prevent invalid Python calls from triggering lower-level `stop` statements for obviously nonphysical inputs.
- In the copied `bspqed_integrals.f90`, the four half-integer Bessel helper wrappers at the bottom were replaced by analytically equivalent spherical-Bessel recurrences so the first wrapped path can build without pulling in a much larger SLATEC special-function branch. The main integral formulas and QED matrix-element logic were left unchanged.

## Recommended Next Wrap Candidates

1. `g_propagator_loop`
   It has the same driver shape as `g_dirac` and returns the same `(irr, red, tot)` table structure.
2. `g_delbrueck`
   Also structurally simple, but introduces the `delbrueck_potential` type and its extra overall scaling factor.
3. `g_uehling` and `g_wichmann_kroll`
   Useful next targets once the simpler two-channel table path is stable, since they return six component channels plus a total.
```bash
cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
python3 workspace_f2py/programs/gfactor_lithium_f2py/python/example_g_dirac.py
```

Run the selectable suite example:

```bash
cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
python3 workspace_f2py/programs/gfactor_lithium_f2py/python/example_gfactor_lithium_suite.py --workflow uehling_all
```

Run the validation comparison:

```bash
cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
python3 workspace_f2py/programs/gfactor_lithium_f2py/validation/validate_g_dirac.py
```

Run the copied-workspace suite validation:

```bash
cd /home/hrishikeshsriman-narayan/Documents/Masters/Master-Thesis--Weak_Int
python3 -u workspace_f2py/programs/gfactor_lithium_f2py/validation/validate_suite.py
```

## Notes On Refactoring

- The protected original driver reads `Z` and `rms` from stdin and writes the convergence table to stdout.
- The wrapper version keeps the same conversions and kernel calls, but returns explicit arrays to Python for the same convergence tables.
- No equations, scaling factors, normalization conventions, units, or B-spline settings were changed for valid inputs.
- Input checks were added only to prevent invalid Python calls from triggering lower-level `stop` statements for obviously nonphysical inputs.
- In the copied `bspqed_integrals.f90`, the four half-integer Bessel helper wrappers at the bottom were replaced by analytically equivalent spherical-Bessel recurrences so the first wrapped path can build without pulling in a much larger SLATEC special-function branch. The main integral formulas and QED matrix-element logic were left unchanged.
- In the copied-workspace `bspqed_uehling.f90` and copied `bspgf_lithium_reference.f90`, nested quadrature callback procedures were refactored into saved-state module procedures so the shared object no longer requires an executable stack when imported from Python. The quadrature formulas themselves were not changed.

## Recommended Next Wrap Candidates

1. `g_propagator_loop`
   It has the same driver shape as `g_dirac` and returns the same `(irr, red, tot)` table structure.
2. `g_delbrueck`
   Also structurally simple, but introduces the `delbrueck_potential` type and its extra overall scaling factor.
3. `g_uehling` and `g_wichmann_kroll`
   Useful next targets once the simpler two-channel table path is stable, since they return six component channels plus a total.
