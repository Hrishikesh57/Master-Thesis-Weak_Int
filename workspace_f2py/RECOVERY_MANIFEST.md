# workspace_f2py recovery manifest

Recovered on 2026-04-28T02:08:36.

Recovered files were rebuilt from local Codex logs because `workspace_f2py` source files were not tracked by Git.

## Written files
- `workspace_f2py/Benchmark_ions.py` (278 lines): strict command-log source
- `workspace_f2py/Direct_dkb_convergence.py` (330 lines): strict command-log source
- `workspace_f2py/Direct_dkb_vs_sampled.py` (263 lines): strict command-log source
- `workspace_f2py/Fermi_route_convergence.py` (430 lines): strict command-log source
- `workspace_f2py/Final_pb_u_neutron_skin_table.py` (303 lines): strict command-log source
- `workspace_f2py/First_order.py` (266 lines): strict command-log source
- `workspace_f2py/Mass_derivative.py` (540 lines): strict command-log source
- `workspace_f2py/Mass_derivative_direct_dkb.py` (1365 lines): strict command-log source
- `workspace_f2py/Validate_fermi_coulomb.py` (269 lines): strict command-log source
- `workspace_f2py/Weak_source_diagnostic.py` (165 lines): strict command-log source
- `workspace_f2py/Mass_derivative_sign_check.py` (252 lines): exact Add File payload from Codex log
- `workspace_f2py/compute_weak_energy_corrections_basis.py` (790 lines): exact Add File payload from Codex log
- `workspace_f2py/programs/gfactor_specific_difference/configs/hlike_pb.toml` (29 lines): exact Add File payload from Codex log
- `workspace_f2py/programs/gfactor_specific_difference/configs/hlike_u.toml` (29 lines): exact Add File payload from Codex log
- `workspace_f2py/programs/gfactor_specific_difference/configs/lilike_pb.toml` (29 lines): exact Add File payload from Codex log
- `workspace_f2py/programs/gfactor_specific_difference/configs/lilike_u.toml` (29 lines): exact Add File payload from Codex log
- `workspace_f2py/programs/gfactor_specific_difference/configs/specific_difference_pb.toml` (9 lines): exact Add File payload from Codex log
- `workspace_f2py/programs/gfactor_specific_difference/configs/specific_difference_u.toml` (9 lines): exact Add File payload from Codex log
- `workspace_f2py/notes/paper_working_report.md` (2962 lines): strict report spans plus Section 42 from current prompt
- `workspace_f2py/outputs/weak_g_benchmark/mass_derivative_sign_check.csv` (5 lines): exact cat output from Codex log

## Known paper report gaps
- original report lines 1173-1269
- original report lines 2366-2419
- original report lines 2759-2759
## Verification
- `python -m py_compile` passed for every recovered `.py` file under `workspace_f2py`.
- `python workspace_f2py/Mass_derivative_sign_check.py` ran successfully and regenerated `outputs/weak_g_benchmark/mass_derivative_sign_check.{json,csv}`.

