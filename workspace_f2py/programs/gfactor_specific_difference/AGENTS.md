# AGENTS.md

## Role

Work as a conservative numerical atomic-physics contributor inside `workspace_f2py/programs/gfactor_specific_difference/`.

## Local Rules

- Reuse copied `workspace_f2py` backends first.
- Do not edit the protected original `workspace/` tree from this subproject.
- Keep units, sign conventions, normalization choices, and nuclear-model definitions explicit.
- Do not silently reinterpret any reused driver output.
- If logic is reused, cite the source path in code comments, docstrings, or documentation.

## Preferred Sources

- `workspace_f2py/programs/gfactor/`
- `workspace_f2py/programs/gfactor_lithium/`
- `workspace_f2py/programs/gfactor_lithium_f2py/`
- `workspace_f2py/programs/rad_f2py/`
- `workspace_f2py/codes/bspqed/src/`
- `Reading/Thesisvf2.pdf`
- `Reading/Cakir_PhD_Thesis.pdf`

## Forbidden Shortcuts

- no silent unit conversions
- no silent gauge switching
- no silent changes to `\xi` definition
- no bulk copying of Fortran library trees when a narrow wrapper or executable bridge suffices
- no dependence on prompt-driven `stdin` programs as the long-term scientific API without wrapping or documenting the interface

## Validation Expectation

When adding or changing numerics here:

- compare against the reused copied backend where possible
- keep the source path of the reused backend visible
- document any provisional physics assumptions immediately
