#!/usr/bin/env sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_dir=$(CDPATH= cd -- "$script_dir/../.." && pwd)
build_dir="$project_dir/build/native-release"
install_dir="$project_dir/bin"

cmake -S "$script_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
cmake --build "$build_dir" --parallel
cmake --install "$build_dir" --prefix "$install_dir" --component Runtime
