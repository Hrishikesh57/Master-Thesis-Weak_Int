#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from _script_support import PROJECT_ROOT, emit_payload
from gfactor_specific_difference import load_ion_case, run_hlike_dirac, summarize_components


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the H-like U 1s baseline through the specific-difference project.")
    parser.add_argument("--config", default=str(PROJECT_ROOT / "configs" / "hlike_u.toml"))
    parser.add_argument("--placeholder", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config)
    case = load_ion_case(config_path)
    if args.placeholder:
        emit_payload(
            {
                "status": "placeholder",
                "script": "run_hlike_u.py",
                "config": str(config_path),
                "label": case.label,
            },
            as_json=args.json,
        )
        return 0

    result = run_hlike_dirac(case)
    emit_payload(summarize_components(result), as_json=args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
