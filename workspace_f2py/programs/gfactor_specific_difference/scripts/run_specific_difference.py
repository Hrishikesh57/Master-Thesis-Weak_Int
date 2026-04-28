#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from _script_support import PROJECT_ROOT, emit_payload
from gfactor_specific_difference import (
    build_specific_difference,
    determine_xi,
    estimate_fns_response,
    estimate_hlike_weak_contribution,
    load_ion_case,
    load_named_run_config,
    run_hlike_dirac,
    run_lilike_dirac,
    summarize_run,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the first-pass Dirac+FNS specific-difference workflow.")
    parser.add_argument("--config", default=str(PROJECT_ROOT / "configs" / "specific_difference_u.toml"))
    parser.add_argument("--placeholder", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config)
    run_cfg = load_named_run_config(config_path)
    hydrogen_case = PROJECT_ROOT / "configs" / str(run_cfg["hydrogen_case"])
    lithium_case = PROJECT_ROOT / "configs" / str(run_cfg["lithium_case"])

    if args.placeholder:
        emit_payload(
            {
                "status": "placeholder",
                "script": "run_specific_difference.py",
                "config": str(config_path),
                "hydrogen_case": str(hydrogen_case),
                "lithium_case": str(lithium_case),
            },
            as_json=args.json,
        )
        return 0

    h_case = load_ion_case(hydrogen_case)
    li_case = load_ion_case(lithium_case)
    h_components = run_hlike_dirac(h_case)
    li_components = run_lilike_dirac(li_case)
    fns_grid = [float(value) for value in run_cfg["fns_rms_scan_fm"]]
    fns_h = estimate_fns_response(h_case, "hlike", fns_grid)
    fns_li = estimate_fns_response(li_case, "lilike", fns_grid)
    xi = determine_xi(li_case, h_case, fns_li, fns_h)
    weak_h = estimate_hlike_weak_contribution(h_case)
    specific_difference = build_specific_difference(li_components, h_components, xi, None, weak_h)
    emit_payload(
        summarize_run(h_components, li_components, fns_h, fns_li, xi, specific_difference, weak_h=weak_h),
        as_json=args.json,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
