from __future__ import annotations

from pathlib import Path
import subprocess
from typing import Iterable


WORKSPACE_F2PY_ROOT = Path(__file__).resolve().parents[2]


def run_command(command: list[str], *, cwd: Path | None = None, stdin_text: str | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        input=stdin_text,
        text=True,
        capture_output=True,
        check=True,
    )


def format_stdin(values: Iterable[object]) -> str:
    return "".join(f"{value}\n" for value in values)


def numeric_rows(stdout: str, expected_columns: int) -> list[tuple[float, ...]]:
    rows: list[tuple[float, ...]] = []
    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != expected_columns:
            continue
        try:
            rows.append(tuple(float(item) for item in parts))
        except ValueError:
            continue
    return rows


def split_gauge_sections(stdout: str) -> tuple[str, str]:
    marker_f = "# Results in Feynman gauge"
    marker_c = "# Results in Coulomb gauge"
    if marker_f not in stdout or marker_c not in stdout:
        raise ValueError("Gauge markers were not found in the program output.")
    after_f = stdout.split(marker_f, 1)[1]
    feynman_text, coulomb_text = after_f.split(marker_c, 1)
    return feynman_text, coulomb_text
