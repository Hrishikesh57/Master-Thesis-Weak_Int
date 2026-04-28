from __future__ import annotations

from pathlib import Path
import sys


PACKAGE_ROOT = Path(__file__).resolve().parent
SRC_ROOT = PACKAGE_ROOT.parent
PROJECT_ROOT = SRC_ROOT.parent
PROGRAMS_ROOT = PROJECT_ROOT.parent
WORKSPACE_F2PY_ROOT = PROGRAMS_ROOT.parent
REPO_ROOT = WORKSPACE_F2PY_ROOT.parent
PYTHON_API_DIR = WORKSPACE_F2PY_ROOT / "python_api"
CONFIG_DIR = PROJECT_ROOT / "configs"
DOCS_DIR = PROJECT_ROOT / "docs"
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"


def ensure_python_api_on_path() -> None:
    python_api = str(PYTHON_API_DIR)
    if python_api not in sys.path:
        sys.path.insert(0, python_api)
