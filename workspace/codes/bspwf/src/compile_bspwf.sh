#!/bin/bash

# Set the module name for the Python extension
MODULE_NAME="bspwf"

# List of all Fortran source files needed for the compilation
SOURCE_FILES="bspwf_constants.f90 bspwf_function.f90 bspwf_dirac.f90"

# Compile the Fortran code into a Python module using f2py
f2py3 -c -m $MODULE_NAME $SOURCE_FILES


