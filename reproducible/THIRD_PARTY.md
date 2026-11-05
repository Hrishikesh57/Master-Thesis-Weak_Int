# Source provenance and attribution

## Reused BspQed solver

The six modules in `workspace/codes/bspqed/src/` are reused from this repository's existing master’s-thesis source. The companion adds the documented local research corrections and a narrow source-built interface; it does not claim authorship of the whole radial DKB solver. The historical source is preserved at the repository root. The local headers do not provide a separate redistribution license, and this update does not assign a blanket license or invent an authorship credit. Existing source notices are retained.

## SLATEC dependency subset

The 37 routines in `workspace/build/slatec/src/` form the dependency closure needed by the included build. Their original authorship and reference headers are retained. Section 4 of the [Netlib Guide to the SLATEC Common Mathematical Library](https://netlib.org/slatec/guide) explicitly describes SLATEC as public-domain software. See [the accompanying notice](workspace/build/slatec/NOTICE.md). These routines are reused numerical infrastructure, not original thesis code.

## Scientific interfaces and corrections

`workspace/programs/weak_spherical/` and `workspace/programs/rad/src/` come from the later local research workspace. The September 2026 correction record distinguishes this baseline from the original submitted thesis. This portfolio update preserves the scientific implementation and adds a build-path fix: F2PY receives a relative kernel source filename so paths containing spaces work. The `scripts/` directory adds example orchestration and a convergence plot.

No private correspondence, unpublished manuscript draft, raw third-party paper PDF, retained compiled extension, or cached environment is added by this update. Public availability of this repository is not a grant of a new project-wide open-source license.
