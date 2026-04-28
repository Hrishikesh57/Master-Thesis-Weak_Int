from __future__ import annotations

from functools import lru_cache
from typing import Any

from ._paths import ensure_python_api_on_path
from .cases import IonCase
from .gfactor_evaluation import GFactorComponents, WorkflowTable, make_dirac_components, make_workflow_table


ensure_python_api_on_path()


def _workspace_api():
    import workspace_f2py_api

    return workspace_f2py_api


def _normalize_rows(rows: list[tuple[float, ...]]) -> tuple[tuple[float, ...], ...]:
    return tuple(tuple(float(item) for item in row) for row in rows)


@lru_cache(maxsize=None)
def _cached_hlike_rows(workflow: str, znucl: int, rms_fm: float, kappa: int, nmain: int) -> tuple[tuple[float, ...], ...]:
    api = _workspace_api()
    result = api.run_gfactor_hydrogen_program(workflow, znucl, rms_fm, kappa, nmain)
    return _normalize_rows(result["rows"])


def _build_lithium_f2py_if_needed() -> None:
    api = _workspace_api()
    try:
        api.load_lithium_suite_module()
    except ModuleNotFoundError:
        api.build_lithium_f2py()


@lru_cache(maxsize=None)
def _cached_lilike_result(workflow: str, znucl: int, rms_fm: float, order: int) -> tuple[Any, ...]:
    api = _workspace_api()
    _build_lithium_f2py_if_needed()
    if workflow == "dirac":
        return tuple(api.compute_lithium_dirac_tables(znucl, rms_fm, order))
    if workflow == "delbrueck":
        return tuple(api.compute_lithium_delbrueck_tables(znucl, rms_fm, order))
    if workflow == "propagator_loop":
        return tuple(api.compute_lithium_propagator_loop_feynman_table(znucl, rms_fm, order))
    if workflow == "uehling":
        return tuple(api.compute_lithium_uehling_tables(znucl, rms_fm, order))
    if workflow == "uehling_all":
        return tuple(api.compute_lithium_uehling_all_tables(znucl, rms_fm, order))
    if workflow == "wichmann_kroll":
        return tuple(api.compute_lithium_wichmann_kroll_tables(znucl, rms_fm, order))
    raise ValueError(f"Unsupported lithium-like workflow: {workflow}")


def _require_hlike(case: IonCase) -> None:
    if not case.is_hlike:
        raise ValueError(f"Expected an H-like case, received {case.charge_state!r}.")


def _require_lilike(case: IonCase) -> None:
    if not case.is_lilike:
        raise ValueError(f"Expected a Li-like case, received {case.charge_state!r}.")


def run_hlike_correction(case: IonCase, workflow: str) -> WorkflowTable:
    _require_hlike(case)
    rows = _cached_hlike_rows(workflow, case.znucl, case.rms_fm, case.kappa, case.nmain)
    if not rows:
        raise RuntimeError(f"No numeric rows were parsed for H-like workflow {workflow!r}.")

    if workflow == "uehling_all":
        column_labels = ("npoints", "g_base", "g_with_uehling", "delta")
        selected_value = rows[-1][-1]
    else:
        column_labels = ("npoints", "value")
        selected_value = rows[-1][-1]

    source_map = {
        "dirac": "workspace_f2py/programs/gfactor/src/g_dirac.f90",
        "delbrueck": "workspace_f2py/programs/gfactor/src/g_delbrueck.f90",
        "uehling": "workspace_f2py/programs/gfactor/src/g_uehling.f90",
        "uehling_all": "workspace_f2py/programs/gfactor/src/g_uehling_all.f90",
        "wichmann_kroll": "workspace_f2py/programs/gfactor/src/g_wichmann_kroll.f90",
    }
    return make_workflow_table(
        workflow=workflow,
        rows=rows,
        column_labels=column_labels,
        selected_value=selected_value,
        source_path=source_map[workflow],
        notes=("Parsed through workspace_f2py_api executable reuse path.",),
    )


def run_hlike_dirac(case: IonCase) -> GFactorComponents:
    table = run_hlike_correction(case, "dirac")
    return make_dirac_components(
        case,
        backend="workspace_f2py_api.run_gfactor_hydrogen_program",
        dirac_table=table,
        notes=(
            "H-like baseline is reused through the copied executable bridge.",
            "The RMS-radius to finite-sphere radius conversion follows the copied gfactor driver convention.",
        ),
        source_paths=(
            "workspace_f2py/python_api/workspace_f2py_api/executables.py",
            table.source_path,
        ),
    )


def _sequence_to_rows(npoints: Any, *columns: Any) -> tuple[tuple[float, ...], ...]:
    npts = [int(item) for item in npoints]
    column_data = [[float(value) for value in column] for column in columns]
    rows = []
    for idx, npoint in enumerate(npts):
        row = [float(npoint)]
        for column in column_data:
            row.append(column[idx])
        rows.append(tuple(row))
    return tuple(rows)


def run_lilike_correction(case: IonCase, workflow: str) -> dict[str, WorkflowTable]:
    _require_lilike(case)
    result = _cached_lilike_result(workflow, case.znucl, case.rms_fm, case.order)
    source_map = {
        "dirac": "workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90",
        "delbrueck": "workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90",
        "propagator_loop": "workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90",
        "uehling": "workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90",
        "uehling_all": "workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90",
        "wichmann_kroll": "workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90",
    }
    tables: dict[str, WorkflowTable] = {}

    if workflow in {"dirac", "delbrueck"}:
        npoints, irr_f, red_f, tot_f, irr_c, red_c, tot_c, *_ = result
        tables["feynman"] = make_workflow_table(
            workflow,
            _sequence_to_rows(npoints, irr_f, red_f, tot_f),
            ("npoints", "irr", "red", "tot"),
            float(tot_f[-1]),
            source_map[workflow],
            gauge="feynman",
        )
        tables["coulomb"] = make_workflow_table(
            workflow,
            _sequence_to_rows(npoints, irr_c, red_c, tot_c),
            ("npoints", "irr", "red", "tot"),
            float(tot_c[-1]),
            source_map[workflow],
            gauge="coulomb",
        )
        return tables

    if workflow == "propagator_loop":
        npoints, irr, red, tot, *_ = result
        tables["feynman"] = make_workflow_table(
            workflow,
            _sequence_to_rows(npoints, irr, red, tot),
            ("npoints", "irr", "red", "tot"),
            float(tot[-1]),
            source_map[workflow],
            gauge="feynman",
        )
        return tables

    if workflow in {"uehling", "wichmann_kroll"}:
        npoints = result[0]
        f_columns = result[1:8]
        c_columns = result[8:15]
        tables["feynman"] = make_workflow_table(
            workflow,
            _sequence_to_rows(npoints, *f_columns),
            ("npoints", "a", "b", "e", "g1", "g2", "g3", "tot"),
            float(f_columns[-1][-1]),
            source_map[workflow],
            gauge="feynman",
        )
        tables["coulomb"] = make_workflow_table(
            workflow,
            _sequence_to_rows(npoints, *c_columns),
            ("npoints", "a", "b", "e", "g1", "g2", "g3", "tot"),
            float(c_columns[-1][-1]),
            source_map[workflow],
            gauge="coulomb",
        )
        return tables

    if workflow == "uehling_all":
        npoints, g_base_f, g_with_f, delta_f, g_base_c, g_with_c, delta_c, *_ = result
        tables["feynman"] = make_workflow_table(
            workflow,
            _sequence_to_rows(npoints, g_base_f, g_with_f, delta_f),
            ("npoints", "g_base", "g_with_uehling", "delta"),
            float(delta_f[-1]),
            source_map[workflow],
            gauge="feynman",
        )
        tables["coulomb"] = make_workflow_table(
            workflow,
            _sequence_to_rows(npoints, g_base_c, g_with_c, delta_c),
            ("npoints", "g_base", "g_with_uehling", "delta"),
            float(delta_c[-1]),
            source_map[workflow],
            gauge="coulomb",
        )
        return tables

    raise ValueError(f"Unsupported lithium-like workflow: {workflow}")


def run_lilike_dirac(case: IonCase) -> GFactorComponents:
    tables = run_lilike_correction(case, "dirac")
    selected_gauge = case.backend.lithium_selected_gauge
    if selected_gauge not in tables:
        raise ValueError(
            f"Requested lithium gauge {selected_gauge!r} is not available for workflow 'dirac'."
        )
    table = tables[selected_gauge]
    return make_dirac_components(
        case,
        backend="workspace_f2py/programs/gfactor_lithium_f2py",
        dirac_table=table,
        selected_gauge=selected_gauge,
        notes=(
            "Li-like baseline is reused through the copied lithium f2py suite.",
            "The copied lithium driver family is used as-is; the final paper-level interpretation of its 'g_dirac' workflow inside the specific-difference stack remains to be validated.",
        ),
        source_paths=(
            "workspace_f2py/programs/gfactor_lithium_f2py/src/gfactor_lithium_suite_f2py.f90",
            table.source_path,
        ),
    )
