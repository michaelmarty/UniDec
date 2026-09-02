#!/usr/bin/env bash

set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_type="${UNIDEC_BUILD_TYPE:-Release}"
build_dir="${UNIDEC_BUILD_DIR:-${repo_dir}/build/native}"

cmake_args=(
    -S "${repo_dir}/unidec/src"
    -B "${build_dir}"
    "-DCMAKE_BUILD_TYPE=${build_type}"
)

if [[ "$(uname -s)" == "Linux" ]]; then
    # Prefer the native system installation. This also prevents WSL from
    # discovering a Windows HDF5 package through its inherited PATH.
    cmake_args+=("-DHDF5_ROOT=${UNIDEC_HDF5_ROOT:-/usr}")
elif [[ "$(uname -s)" == "Darwin" ]] && command -v brew >/dev/null 2>&1; then
    hdf5_prefix="$(brew --prefix hdf5)"
    fftw_prefix="$(brew --prefix fftw)"
    libomp_prefix="$(brew --prefix libomp)"
    export PKG_CONFIG_PATH="${fftw_prefix}/lib/pkgconfig:${hdf5_prefix}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
    cmake_args+=(
        "-DCMAKE_PREFIX_PATH=${hdf5_prefix};${fftw_prefix};${libomp_prefix}"
        "-DOpenMP_C_FLAGS=-Xpreprocessor -fopenmp"
        "-DOpenMP_C_LIB_NAMES=omp"
        "-DOpenMP_omp_LIBRARY=${libomp_prefix}/lib/libomp.dylib"
        "-DOpenMP_C_INCLUDE_DIR=${libomp_prefix}/include"
    )
fi

cmake "${cmake_args[@]}"
cmake --build "${build_dir}" --config "${build_type}" --parallel

case "$(uname -s)" in
    Darwin) executable="${repo_dir}/unidec/bin/unidecmac" ;;
    Linux) executable="${repo_dir}/unidec/bin/unideclinux" ;;
    *)
        echo "Unsupported platform for build_native.sh: $(uname -s)" >&2
        exit 2
        ;;
esac

if [[ ! -x "${executable}" ]]; then
    echo "Native executable was not deployed: ${executable}" >&2
    exit 1
fi

echo "Built ${executable}"
