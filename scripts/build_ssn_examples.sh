#!/usr/bin/env bash
# SPDX-License-Identifier: MPL-2.0
set -euo pipefail
cd "$(dirname "$0")/.."

ssn_build="${DPSIM_SSN_BUILD:-build_ssn_examples}"
extra=()
if [[ -n "${DPSIM_DEPENDENCY_CACHE:-}" ]]; then
  extra+=("-DFETCHCONTENT_SOURCE_DIR_SUITESPARSE-MODULE=$DPSIM_DEPENDENCY_CACHE/suitesparse-module-src")
  extra+=("-DFETCHCONTENT_SOURCE_DIR_JSON=$DPSIM_DEPENDENCY_CACHE/json-src")
  extra+=("-DFETCHCONTENT_SOURCE_DIR_READERWRITERQUEUE-MODULE=$DPSIM_DEPENDENCY_CACHE/readerwriterqueue-module-src")
fi

cmake -S . -B "$ssn_build" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_EXPORT_NO_PACKAGE_REGISTRY=ON \
  -DFETCH_GRID_DATA=OFF -DDPSIM_BUILD_DOC=OFF -DWITH_PYBIND=OFF \
  -DWITH_CIM=OFF -DWITH_VILLAS=OFF -DWITH_SUNDIALS=OFF -DWITH_OPENMP=OFF \
  -DWITH_SPARSE=ON -DWITH_KLU=ON \
  -DWITH_MARCH_NATIVE=ON -DWITH_LTO=OFF -DWITH_PROFILING=OFF "${extra[@]}"

targets=("$@")
if [[ ${#targets[@]} -eq 0 ]]; then
  for grid in 9 39 118 300; do
    for variant in Accuracy Dense Sparse; do
      targets+=("IEEE${grid}_SSN_${variant}")
    done
  done
fi
cmake --build "$ssn_build" --target "${targets[@]}" -j "${BUILD_JOBS:-4}"
