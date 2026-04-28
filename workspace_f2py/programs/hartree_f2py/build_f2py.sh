#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_COPY_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SRC_DIR="${SCRIPT_DIR}/src"
LIB_DIR="${WORKSPACE_COPY_DIR}/include/lib"
MOD_DIR="${WORKSPACE_COPY_DIR}/include/mod"
BUILD_DIR="${SCRIPT_DIR}/build"
MODULE_NAME="hartree_suite"
EXT_SUFFIX="$(python3 - <<'PY'
import sysconfig
print(sysconfig.get_config_var('EXT_SUFFIX'))
PY
)"
PYTHON_INCLUDE="$(python3 - <<'PY'
import sysconfig
print(sysconfig.get_paths()['include'])
PY
)"
NUMPY_INCLUDE="$(python3 - <<'PY'
import numpy
print(numpy.get_include())
PY
)"
F2PY_INCLUDE="$(python3 - <<'PY'
import numpy.f2py
print(numpy.f2py.get_include())
PY
)"

cd "${SCRIPT_DIR}"

rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

python3 -m numpy.f2py \
    --build-dir "${BUILD_DIR}" \
    -m "${MODULE_NAME}" \
    "${SRC_DIR}/dhfcb_f2py.f90" \
    "${SRC_DIR}/shf_f2py.f90"

gcc \
    -fPIC \
    -c \
    "${BUILD_DIR}/${MODULE_NAME}module.c" \
    -I"${PYTHON_INCLUDE}" \
    -I"${NUMPY_INCLUDE}" \
    -I"${F2PY_INCLUDE}" \
    -o "${BUILD_DIR}/${MODULE_NAME}module.o"

gcc \
    -fPIC \
    -c \
    "${F2PY_INCLUDE}/fortranobject.c" \
    -I"${PYTHON_INCLUDE}" \
    -I"${NUMPY_INCLUDE}" \
    -I"${F2PY_INCLUDE}" \
    -o "${BUILD_DIR}/fortranobject.o"

gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" -I"${MOD_DIR}" \
    -c "${SRC_DIR}/dhfcb_f2py.f90" -o "${BUILD_DIR}/dhfcb_f2py.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" -I"${MOD_DIR}" \
    -c "${SRC_DIR}/shf_f2py.f90" -o "${BUILD_DIR}/shf_f2py.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BUILD_DIR}/${MODULE_NAME}-f2pywrappers2.f90" -o "${BUILD_DIR}/${MODULE_NAME}-f2pywrappers2.o"

gfortran -shared -o "${SCRIPT_DIR}/${MODULE_NAME}${EXT_SUFFIX}" \
    "${BUILD_DIR}/${MODULE_NAME}module.o" \
    "${BUILD_DIR}/fortranobject.o" \
    "${BUILD_DIR}/dhfcb_f2py.o" \
    "${BUILD_DIR}/shf_f2py.o" \
    "${BUILD_DIR}/${MODULE_NAME}-f2pywrappers2.o" \
    -L"${LIB_DIR}" \
    -lbspwf \
    -lslatec \
    -llapack \
    $(python3-config --ldflags)
