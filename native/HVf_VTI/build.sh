#!/bin/sh
# Build source only; no scientific forward model is run here.
set -eu
cd "$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
mkdir -p bin build
cd build
"${FC:-gfortran}" -O2 -ffree-line-length-none -I.. ../HV.f90 -o ../bin/HVf_VTI
cd ..
shasum -a 256 ./*.f90 ./*.f03 > bin/source.sha256
shasum -a 256 bin/HVf_VTI > bin/executable.sha256
