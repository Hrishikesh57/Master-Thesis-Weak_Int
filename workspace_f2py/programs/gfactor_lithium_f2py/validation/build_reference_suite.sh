#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_COPY_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BIN_DIR="${SCRIPT_DIR}/bin"
MOD_DIR="${WORKSPACE_COPY_DIR}/include/mod"
LIB_DIR="${WORKSPACE_COPY_DIR}/include/lib"
SRC_DIR="${WORKSPACE_COPY_DIR}/programs/gfactor_lithium/src"

mkdir -p "${BIN_DIR}"

for name in g_dirac g_delbrueck g_propagator_loop g_uehling g_uehling_all g_wichmann_kroll; do
    gfortran \
        -Wall \
        -O2 \
        -I"${MOD_DIR}" \
        -L"${LIB_DIR}" \
        -o "${BIN_DIR}/${name}" \
        "${SRC_DIR}/${name}.f90" \
        -lbspgf \
        -lbspqed \
        -lbspwf \
        -lslatec \
        -llapack
done
