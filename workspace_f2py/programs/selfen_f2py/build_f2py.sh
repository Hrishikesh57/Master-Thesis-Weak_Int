#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_COPY_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SRC_DIR="${SCRIPT_DIR}/src"
LI_F2PY_SRC_COPY_DIR="${WORKSPACE_COPY_DIR}/programs/gfactor_lithium_f2py/src/copied"
BSPQED_SRC_DIR="${WORKSPACE_COPY_DIR}/codes/bspqed/src"
BSPSE_SRC_DIR="${WORKSPACE_COPY_DIR}/codes/bspse/src"
LIB_DIR="${WORKSPACE_COPY_DIR}/include/lib"
BUILD_DIR="${SCRIPT_DIR}/build"
MODULE_NAME="selfen_suite"
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
    "${SRC_DIR}/selfen_suite_f2py.f90"

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

gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${LI_F2PY_SRC_COPY_DIR}/bspqed_function.f90" -o "${BUILD_DIR}/bspqed_function.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${LI_F2PY_SRC_COPY_DIR}/bspqed_physconsts.f90" -o "${BUILD_DIR}/bspqed_physconsts.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${LI_F2PY_SRC_COPY_DIR}/bspqed_quadrature.f90" -o "${BUILD_DIR}/bspqed_quadrature.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPQED_SRC_DIR}/bspqed_uehling.f90" -o "${BUILD_DIR}/bspqed_uehling.o"
gfortran -Wall -O2 -fPIC \
    -c "${WORKSPACE_COPY_DIR}/programs/gfactor_lithium_f2py/src/xermsg_stub.f90" -o "${BUILD_DIR}/xermsg_stub.o"
gfortran -Wall -O2 -fPIC \
    -c "${LI_F2PY_SRC_COPY_DIR}/dintrv.f" -o "${BUILD_DIR}/dintrv.o"
gfortran -Wall -O2 -fPIC \
    -c "${LI_F2PY_SRC_COPY_DIR}/dbspvn.f" -o "${BUILD_DIR}/dbspvn.o"
gfortran -Wall -O2 -fPIC \
    -c "${LI_F2PY_SRC_COPY_DIR}/dbspvd.f" -o "${BUILD_DIR}/dbspvd.o"
gfortran -Wall -O2 -fPIC \
    -c "${LI_F2PY_SRC_COPY_DIR}/dbvalu.f" -o "${BUILD_DIR}/dbvalu.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${LI_F2PY_SRC_COPY_DIR}/bspqed_dkb.f90" -o "${BUILD_DIR}/bspqed_dkb.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${LI_F2PY_SRC_COPY_DIR}/bspqed_integrals.f90" -o "${BUILD_DIR}/bspqed_integrals.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPQED_SRC_DIR}/bspqed_potentials.f90" -o "${BUILD_DIR}/bspqed_potentials.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPQED_SRC_DIR}/bspqed_utilities.f90" -o "${BUILD_DIR}/bspqed_utilities.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPQED_SRC_DIR}/bspqed_coefficients.f90" -o "${BUILD_DIR}/bspqed_coefficients.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPQED_SRC_DIR}/bspqed_fourier_density.f90" -o "${BUILD_DIR}/bspqed_fourier_density.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPQED_SRC_DIR}/bspqed_momentum.f90" -o "${BUILD_DIR}/bspqed_momentum.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPQED_SRC_DIR}/bspqed_discretization.f90" -o "${BUILD_DIR}/bspqed_discretization.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPSE_SRC_DIR}/bspse_zeropotential.f90" -o "${BUILD_DIR}/bspse_zeropotential.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPSE_SRC_DIR}/bspse_onepotential.f90" -o "${BUILD_DIR}/bspse_onepotential.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BSPSE_SRC_DIR}/bspse_manypotential.f90" -o "${BUILD_DIR}/bspse_manypotential.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${SRC_DIR}/selfen_suite_f2py.f90" -o "${BUILD_DIR}/selfen_suite_f2py.o"
gfortran -Wall -O2 -fPIC -J"${BUILD_DIR}" -I"${BUILD_DIR}" \
    -c "${BUILD_DIR}/${MODULE_NAME}-f2pywrappers2.f90" -o "${BUILD_DIR}/${MODULE_NAME}-f2pywrappers2.o"

gfortran -shared -o "${SCRIPT_DIR}/${MODULE_NAME}${EXT_SUFFIX}" \
    "${BUILD_DIR}/${MODULE_NAME}module.o" \
    "${BUILD_DIR}/fortranobject.o" \
    "${BUILD_DIR}/bspqed_function.o" \
    "${BUILD_DIR}/bspqed_physconsts.o" \
    "${BUILD_DIR}/bspqed_quadrature.o" \
    "${BUILD_DIR}/bspqed_uehling.o" \
    "${BUILD_DIR}/xermsg_stub.o" \
    "${BUILD_DIR}/dintrv.o" \
    "${BUILD_DIR}/dbspvn.o" \
    "${BUILD_DIR}/dbspvd.o" \
    "${BUILD_DIR}/dbvalu.o" \
    "${BUILD_DIR}/bspqed_dkb.o" \
    "${BUILD_DIR}/bspqed_integrals.o" \
    "${BUILD_DIR}/bspqed_potentials.o" \
    "${BUILD_DIR}/bspqed_utilities.o" \
    "${BUILD_DIR}/bspqed_coefficients.o" \
    "${BUILD_DIR}/bspqed_fourier_density.o" \
    "${BUILD_DIR}/bspqed_momentum.o" \
    "${BUILD_DIR}/bspqed_discretization.o" \
    "${BUILD_DIR}/bspse_zeropotential.o" \
    "${BUILD_DIR}/bspse_onepotential.o" \
    "${BUILD_DIR}/bspse_manypotential.o" \
    "${BUILD_DIR}/selfen_suite_f2py.o" \
    "${BUILD_DIR}/${MODULE_NAME}-f2pywrappers2.o" \
    -L"${LIB_DIR}" \
    -lslatec \
    -llapack \
    -lblas \
    $(python3-config --ldflags)
