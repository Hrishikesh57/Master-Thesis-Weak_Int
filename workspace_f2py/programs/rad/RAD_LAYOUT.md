`rad` is organized so reference inputs live under `data/` and generated artifacts live under `outputs/`.

Top-level directories:

- `src/`: Fortran sources and helper batch scripts.
- `bin/`: built executables.
- `data/nuclear/`: nuclear reference inputs such as `radii_ref.txt`.
- `data/reference/radial_functions/`: analytic or converted reference radial tables used by `something`.
- `data/reference/energy_states/hydrogen_like/`: reference energy-state tables.
- `data/reference/wavefunctions/single_electron_legacy/`: older single-electron wavefunction tables kept as legacy references.
- `data/reference/wavefunctions/scf_legacy/`: older SCF wavefunction tables kept as legacy references.
- `data/grasp/raw/`: raw GRASP exports.
- `data/grasp/rel/`: GRASP data already converted to relativistic units.
- `outputs/wavefunctions/single_electron/`: new single-electron tables written by `practise`.
- `outputs/wavefunctions/scf/`: new SCF tables written by `practise`.
- `outputs/plots/`: generated diagnostic and comparison plots.
- `outputs/results/energy/hydrogen_like/`: generated second-order energy correction tables.
- `outputs/results/amplitudes/`: generated amplitude reports.
- `outputs/analysis/`: loose analysis tables and conversion outputs that are not source data.
- `outputs/figures/legacy/`: legacy figure files that were previously loose in the repo root.

Current lookup behavior:

- `src/practise.f90` writes fresh wavefunctions into `outputs/wavefunctions/...`.
- `src/something.f90` looks first in `outputs/wavefunctions/...`, then in `data/reference/...`, then in the old legacy paths as fallbacks.
- `plot_radial_table.py` resolves filenames from both `outputs/` and `data/` locations and saves plots into `outputs/plots/`.
