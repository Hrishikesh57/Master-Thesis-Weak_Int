# `hartree_f2py`

Copied `f2py` wrapper path for the Hartree programs in `workspace_f2py/programs/hartree`.

What is wrapped:

- `hartree_dhfcb_f2py.compute_dhfcb_energies`
  - copied from the trusted `dhfcb` workflow
  - keeps the original Fermi-charge setup: `order=9`, `npoints=100`, `rcavity=50`
- `hartree_shf_f2py.compute_shf_energies`
  - copied from the trusted `shf` workflow
  - keeps the original Schrödinger Hartree-Fock setup: `order=6`, `npoints=36`, `rcavity=40`

What is not wrapped yet:

- `dhfc`
  - the copied `dhfc` source still mixes a `bspqed` sphere potential with `bspwf` interfaces expecting a different function hierarchy
  - that needs an explicit physics-aware resolution before it is safe to expose as a direct Python module

Build:

```bash
bash workspace_f2py/programs/hartree_f2py/build_f2py.sh
```

Example:

```bash
python3 workspace_f2py/programs/hartree_f2py/python/example_hartree_suite.py
```

Direct Python use:

```python
import sys
sys.path.insert(0, "workspace_f2py/programs/hartree_f2py")
import hartree_suite

core_e, valence_e, ierr = hartree_suite.hartree_dhfcb_f2py.compute_dhfcb_energies(
    10, 3.0, [1], [-1], [2], [-1]
)

core_e, valence_e, total_e, ierr = hartree_suite.hartree_shf_f2py.compute_shf_energies(
    2, [1], [0], [2], [0]
)
```
