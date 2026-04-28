#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import selfen_suite  # noqa: E402


def main() -> None:
    wrapper = selfen_suite.selfen_suite_f2py

    zero = wrapper.compute_zero_potential_self_energy(92, 5.8571, -1, 1, 80, 1.0e-4, 1.0e-4)
    print("Zero-potential:", zero[0], zero[1], zero[-1])

    one = wrapper.compute_one_potential_self_energy(92, 5.8571, -1, 1, 80, 1.0e-4, 1.0e-4)
    print("One-potential:", one[0], one[1], one[-1])

    many = wrapper.compute_many_potential_self_energy(92, 5.8571, -1, 1, 60, 1.0e-4, 1.0e-4, 1, 2, 2)
    print("Many-potential rows:")
    for row in zip(many[0], many[1], many[2], many[3], many[4]):
        print(row)


if __name__ == "__main__":
    main()
