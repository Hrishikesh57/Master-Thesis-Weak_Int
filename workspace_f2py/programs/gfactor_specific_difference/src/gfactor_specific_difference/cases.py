from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping
import tomllib


@dataclass(frozen=True)
class NuclearModelConfig:
    charge_model: str = "sphere"
    weak_model: str = "uniform_sphere"
    charge_radius_convention: str = "rms_fm"
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class BackendConfig:
    hlike_mode: str = "executable_bridge"
    lilike_mode: str = "direct_f2py"
    lithium_selected_gauge: str = "feynman"
    order: int = 9
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class IonCase:
    label: str
    charge_state: str
    isotope: str
    znucl: int
    nneutrons: int
    rms_fm: float
    kappa: int
    nmain: int
    state_label: str
    order: int = 9
    nuclear_model: NuclearModelConfig = field(default_factory=NuclearModelConfig)
    backend: BackendConfig = field(default_factory=BackendConfig)
    source_path: str | None = None
    notes: tuple[str, ...] = ()

    def with_rms_fm(self, rms_fm: float) -> "IonCase":
        return replace(self, rms_fm=rms_fm)

    @property
    def normalized_charge_state(self) -> str:
        return self.charge_state.strip().lower().replace("-", "").replace(" ", "")

    @property
    def is_hlike(self) -> bool:
        return self.normalized_charge_state == "hlike"

    @property
    def is_lilike(self) -> bool:
        return self.normalized_charge_state == "lilike"


def load_toml_document(path: str | Path) -> dict[str, Any]:
    path_obj = Path(path)
    with path_obj.open("rb") as handle:
        return tomllib.load(handle)


def load_ion_case(path: str | Path) -> IonCase:
    path_obj = Path(path)
    data = load_toml_document(path_obj)
    ion = data.get("ion", {})
    nuclear = data.get("nuclear_model", {})
    backend = data.get("backend", {})

    nuclear_cfg = NuclearModelConfig(
        charge_model=str(nuclear.get("charge_model", "sphere")),
        weak_model=str(nuclear.get("weak_model", "uniform_sphere")),
        charge_radius_convention=str(nuclear.get("charge_radius_convention", "rms_fm")),
        notes=tuple(str(item) for item in nuclear.get("notes", [])),
    )
    backend_cfg = BackendConfig(
        hlike_mode=str(backend.get("hlike_mode", "executable_bridge")),
        lilike_mode=str(backend.get("lilike_mode", "direct_f2py")),
        lithium_selected_gauge=str(backend.get("lithium_selected_gauge", "feynman")),
        order=int(backend.get("order", ion.get("order", 9))),
        notes=tuple(str(item) for item in backend.get("notes", [])),
    )
    return IonCase(
        label=str(ion["label"]),
        charge_state=str(ion["charge_state"]),
        isotope=str(ion["isotope"]),
        znucl=int(ion["znucl"]),
        nneutrons=int(ion["nneutrons"]),
        rms_fm=float(ion["rms_fm"]),
        kappa=int(ion["kappa"]),
        nmain=int(ion["nmain"]),
        state_label=str(ion["state_label"]),
        order=int(ion.get("order", backend_cfg.order)),
        nuclear_model=nuclear_cfg,
        backend=backend_cfg,
        source_path=str(path_obj),
        notes=tuple(str(item) for item in ion.get("notes", [])),
    )


def load_named_run_config(path: str | Path) -> Mapping[str, Any]:
    data = load_toml_document(path)
    return data.get("run", {})
