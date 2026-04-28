## `workspace_f2py/python_api`

This directory provides a Python-first interface to the copied `workspace_f2py/` tree.

There are two interface modes:

- Direct `f2py`:
  - `workspace_f2py_api.load_lithium_suite_module()` loads `programs/gfactor_lithium_f2py/gfactor_lithium_suite`.
- Python-to-executable bridges:
  - `workspace_f2py_api.run_gfactor_hydrogen_program(...)`
  - `workspace_f2py_api.run_gfactor_boron_program(...)`
  - `workspace_f2py_api.run_gfactor_lithium_program(...)`
  - `workspace_f2py_api.run_selfen_zero(...)`
  - `workspace_f2py_api.run_selfen_one(...)`
  - `workspace_f2py_api.run_selfen_many(...)`
  - `workspace_f2py_api.run_hartree_program(...)`
  - `workspace_f2py_api.run_rad_program(...)`

Build helpers:

- `workspace_f2py_api.build_lithium_f2py()` builds the copied lithium `f2py` module.
- `workspace_f2py_api.build_program_target(program_name, *targets)` builds selected copied executables.
- `workspace_f2py_api.build_python_backend()` runs `workspace_f2py/build_python_backend.sh`.

Example:

```bash
python3 workspace_f2py/python_api/examples/run_workspace_from_python.py
```

Current validation status:

- copied `gfactor` drivers build through their own copied Makefiles
- copied `selfen/zero` builds and runs through the Python bridge
- copied lithium `f2py` module imports through the Python package
- copied `hartree/dhfc` still fails at compile time due an existing source-level type mismatch in the copied Fortran
  - `workspace_f2py_api.run_rad_program(...)`

Build helpers:

- `workspace_f2py_api.build_lithium_f2py()` builds the copied lithium `f2py` module.
- `workspace_f2py_api.build_selfen_f2py()` builds the copied self-energy `f2py` module.
- `workspace_f2py_api.build_hartree_f2py()` builds the copied Hartree `f2py` module.
- `workspace_f2py_api.build_rad_f2py()` builds the copied `rad` kernel `f2py` module.
- `workspace_f2py_api.build_f2py_target(name)` builds one copied `f2py` target by name.
- `workspace_f2py_api.build_all_f2py()` builds all copied direct `f2py` modules.
- `workspace_f2py_api.build_program_target(program_name, *targets)` builds selected copied executables.
- `workspace_f2py_api.build_python_backend()` runs `workspace_f2py/build_python_backend.sh`.

Example:

```bash
python3 workspace_f2py/python_api/examples/run_workspace_from_python.py
```

Current validation status:

- copied `gfactor` drivers build through their own copied Makefiles
- copied `selfen/zero` builds and runs through the Python bridge
- copied lithium, self-energy, Hartree, and sampled-`rad` `f2py` modules import through the Python package
- copied `hartree/dhfc` still fails at compile time due an existing source-level type mismatch in the copied Fortran
