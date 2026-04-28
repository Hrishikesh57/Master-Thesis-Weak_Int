import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import hartree_suite


def main() -> None:
    shf = hartree_suite.hartree_shf_f2py.compute_shf_energies(
        2,
        [1],
        [0],
        [2],
        [0],
    )
    print("SHF core energies:", shf[0])
    print("SHF valence energies:", shf[1])
    print("SHF total core energy:", shf[2])
    print("SHF ierr:", shf[3])

    dhfcb = hartree_suite.hartree_dhfcb_f2py.compute_dhfcb_energies(
        10,
        3.0,
        [1],
        [-1],
        [2],
        [-1],
    )
    print("DHFCB core energies:", dhfcb[0])
    print("DHFCB valence energies:", dhfcb[1])
    print("DHFCB ierr:", dhfcb[2])


if __name__ == "__main__":
    main()
