# Validation record — 5 November 2026

## Native environment

A fresh virtual environment was created without global site packages. The installed versions were CPython 3.14.8, NumPy 2.4.6, SciPy 1.16.2, Matplotlib 3.10.8 and pytest 9.0.3 on Linux x86_64. The pinned requirements include transitive Python dependencies. Use the installation and run commands in the repository README; choose a fresh output directory for each run.

10 physics tests passed after a source-only build. The convergence driver passed all 58 comparisons across 16 Pb-208/U-238 calculations. It varied breakpoint count, spline order, quadrature and cavity radius; the independent mass-derivative and magnetic-response comparison also passed.

## Container environment

The companion Dockerfile was built and run with rootless Podman, using CPython 3.14.10 and the same Python dependency pins. Build and runtime both returned exit code zero. The thesis container rebuilt the Fortran/C extension from the included source and repeated all 58 numerical comparisons. Compiler flags and source/binary hashes are recorded in the generated build manifest.

The local verification image identifier is `98781ca02c4fbacda42920cfcd21626ce49254ecc863e0a830f86aeaefbc2adf`; this image has not been published to a registry. The demo's recorded observables agreed with the native example at the saved precision. The generated figures were visually inspected.

## Reproducibility limits

Native Docker Engine was unavailable, so the tested container engine was Podman with Dockerfile-compatible builds. Native Windows/macOS execution was not tested. The active GitHub Actions definition is at the repository root. The recorded native/container checks preceded the public commit; the live Actions result is recorded on GitHub. The pure-Python workflow's additional Python 3.12 job is a future check, not a claimed result.

The base image tag and native BLAS/compiler packages can change. The Python requirements are pinned, but bitwise-identical builds across platforms are not promised. The reference data are calculations of the stated model, not experimental observations. Scientific assumptions and remaining interpretation limits are described in the model report.

Portable reference JSON retains source hashes and scientific results but omits the original machine-specific compiler command paths. Rebuilding writes the full new command list to `workspace/programs/weak_spherical/build_manifest.json`.
