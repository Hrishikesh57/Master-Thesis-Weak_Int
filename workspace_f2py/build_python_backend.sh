#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

make -C "${SCRIPT_DIR}/programs/gfactor" \
    g_dirac g_delbrueck g_uehling g_uehling_all g_wichmann_kroll

make -C "${SCRIPT_DIR}/programs/gfactor_boron" \
    g_dirac g_delbrueck g_propagator_loop g_uehling g_uehling_all g_wichmann_kroll

make -C "${SCRIPT_DIR}/programs/gfactor_lithium" \
    g_dirac g_delbrueck g_propagator_loop g_uehling g_uehling_all g_wichmann_kroll

make -C "${SCRIPT_DIR}/programs/selfen" zero one many
make -C "${SCRIPT_DIR}/programs/hartree" dhfcb dhfc shf
make -C "${SCRIPT_DIR}/programs/rad" practise something theoretical_radial_functions_1s_2s bspwf_h

"${SCRIPT_DIR}/programs/gfactor_lithium_f2py/build_f2py.sh"
