#! /bin/bash

cmake \
  --no-warn-unused-cli \
  -S . \
  -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE:STRING=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE \
  -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ \
  -DCOCOA_ENABLE_MED=OFF \
  -DCOCOA_ENABLE_PROXPDE=OFF \
  -DCOCOA_ENABLE_OFORG=ON \
  $*

#   -DProXPDE_FETCH=False \
#   -DProXPDE_DIR:PATH=~/software/proxpde/lib/cmake/ProXPDE \
#   -DPYBIND11_PYTHON_VERSION=3.11 \
