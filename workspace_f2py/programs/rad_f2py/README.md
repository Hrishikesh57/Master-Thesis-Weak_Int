# `rad_f2py`

Copied `f2py` wrapper path for the sampled-kernel part of `workspace_f2py/programs/rad`.

What is wrapped:

- `rad_kernel_f2py.weak_potential_summary`
  - builds the weak source with the same `R_fm -> R_rel` conversion and weak normalization used by `something.f90`
- `rad_kernel_f2py.first_order_weak_shift`
  - sampled one-channel weak energy integral from `run_first_order_energy_correction`
- `rad_kernel_f2py.second_order_shift_sampled`
  - sampled two-channel weak matrix element and second-order energy correction from `run_second_order_energy_correction`
- `rad_kernel_f2py.discrete_amp_sampled`
  - one-intermediate sampled discrete amplitude contribution from the `amplitudes` path in `something.f90`
  - keeps the original Stark and weak amplitude finalization, including the `A_W / (3 sqrt(3))` normalization

What is not wrapped yet:

- the interactive file-discovery and prompt-driven top-level workflow in `something.f90`
- the basis-spectrum amplitude path with `basis_spectrum_data`
- the SCF/single-electron wavefunction generation path in `practise.f90`

Build:

```bash
bash workspace_f2py/programs/rad_f2py/build_f2py.sh
```

Example and validation:

```bash
python3 workspace_f2py/programs/rad_f2py/python/example_rad_kernel.py
```

The example:

- loads the copied sampled radial tables from `workspace_f2py/programs/rad`
- compares the wrapped first-order weak shift against the copied `bin/something` executable
- compares a wrapped `55: 6s -> 7s` single-intermediate `6p` amplitude contribution against the copied executable
