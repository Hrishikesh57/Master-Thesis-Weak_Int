from __future__ import annotations

from dataclasses import dataclass, field

from .cases import IonCase


@dataclass
class WorkflowTable:
    workflow: str
    rows: tuple[tuple[float, ...], ...]
    column_labels: tuple[str, ...]
    selected_value: float
    gauge: str | None
    source_path: str
    notes: tuple[str, ...] = ()

    @property
    def npoints(self) -> tuple[int, ...]:
        if not self.rows:
            return ()
        return tuple(int(row[0]) for row in self.rows)


@dataclass
class GFactorComponents:
    case: IonCase
    backend: str
    selected_gauge: str | None
    dirac_value: float
    corrections: dict[str, float] = field(default_factory=dict)
    workflow_tables: dict[str, WorkflowTable] = field(default_factory=dict)
    source_paths: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    @property
    def total_value(self) -> float:
        return self.dirac_value + sum(self.corrections.values())

    @property
    def selected_value(self) -> float:
        return self.total_value


def make_workflow_table(
    workflow: str,
    rows: tuple[tuple[float, ...], ...],
    column_labels: tuple[str, ...],
    selected_value: float,
    source_path: str,
    *,
    gauge: str | None = None,
    notes: tuple[str, ...] = (),
) -> WorkflowTable:
    return WorkflowTable(
        workflow=workflow,
        rows=rows,
        column_labels=column_labels,
        selected_value=selected_value,
        gauge=gauge,
        source_path=source_path,
        notes=notes,
    )


def make_dirac_components(
    case: IonCase,
    backend: str,
    dirac_table: WorkflowTable,
    *,
    selected_gauge: str | None = None,
    notes: tuple[str, ...] = (),
    source_paths: tuple[str, ...] = (),
) -> GFactorComponents:
    return GFactorComponents(
        case=case,
        backend=backend,
        selected_gauge=selected_gauge,
        dirac_value=dirac_table.selected_value,
        workflow_tables={"dirac": dirac_table},
        source_paths=source_paths or (dirac_table.source_path,),
        notes=notes,
    )


def add_correction(
    components: GFactorComponents,
    name: str,
    value: float,
    table: WorkflowTable,
) -> None:
    components.corrections[name] = value
    components.workflow_tables[name] = table


def summarize_components(components: GFactorComponents) -> dict[str, object]:
    return {
        "label": components.case.label,
        "charge_state": components.case.charge_state,
        "isotope": components.case.isotope,
        "znucl": components.case.znucl,
        "nneutrons": components.case.nneutrons,
        "state_label": components.case.state_label,
        "rms_fm": components.case.rms_fm,
        "selected_gauge": components.selected_gauge,
        "dirac_value": components.dirac_value,
        "corrections": components.corrections,
        "total_value": components.total_value,
        "workflow_tables": {
            name: {
                "workflow": table.workflow,
                "gauge": table.gauge,
                "selected_value": table.selected_value,
                "npoints": list(table.npoints),
                "source_path": table.source_path,
            }
            for name, table in components.workflow_tables.items()
        },
        "source_paths": list(components.source_paths),
        "notes": list(components.notes),
    }
