# `selfen` `f2py` Bridge

This directory provides a direct Python interface to the copied `workspace_f2py/programs/selfen` workflow without modifying that reference driver family in place.

Exposed workflows:

- `compute_zero_potential_self_energy`
- `compute_one_potential_self_energy`
- `compute_many_potential_self_energy`

Build:

```bash
bash workspace_f2py/programs/selfen_f2py/build_f2py.sh
```

Example:

```bash
python3 workspace_f2py/programs/selfen_f2py/python/example_selfen_suite.py
```
