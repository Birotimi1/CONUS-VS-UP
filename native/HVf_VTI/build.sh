#!/bin/sh
# Build source only; no scientific forward model is run here.
set -eu
cd "$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
mkdir -p bin build
cd build
"${FC:-gfortran}" -O2 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -I.. ../HV.f90 -o ../bin/HVf_VTI
cd ..
# Linux/Expanse normally provides sha256sum; macOS provides shasum.
if command -v sha256sum >/dev/null 2>&1; then
    sha256sum ./*.f90 ./*.f03 build.sh > bin/source.sha256
    sha256sum bin/HVf_VTI > bin/executable.sha256
else
    shasum -a 256 ./*.f90 ./*.f03 build.sh > bin/source.sha256
    shasum -a 256 bin/HVf_VTI > bin/executable.sha256
fi
