# UniDec Native Engine Guide

This file applies to `public/UniDec/unidec/`, including the C/C++ engine in
`src/` and the package-consumed artifacts in `bin/`.

## Native Code

- CMake requires version 3.25, HDF5 C/HL, and single-precision FFTW or MKL;
  OpenMP is optional.
- Keep declarations, definitions, command dispatch, and CMake source lists in
  sync. Maintain portable C11/C++ across Windows, Linux, macOS, x64, and ARM64,
  guarding platform-specific behavior explicitly.
- Use existing CMake options such as `UNIDEC_ENABLE_OPENMP` and
  `UNIDEC_REQUIRE_INTELLLVM`; never hard-code a local toolchain.
- Treat array-size, HDF5-type, signedness, and narrowing warnings as possible
  data-corruption defects. For numerical work, test empty inputs, singleton
  dimensions, invalid values, and non-finite data.
- Normal builds copy the executable and shared library into `bin/`. Treat these
  tracked files as build outputs and review their changes separately.
- Never reuse an unknown or stale CMake cache to judge a source change. Use a
  clearly named, task-specific build directory outside the source tree.

## Native Validation

From `public/UniDec`, build the packaged executable on Linux or macOS with:

```shell
bash ./build_native.sh
```

On Windows with IntelLLVM, run this from the repository root:

```shell
public\UniDec\build_native_intel.bat
```

The wrapper anchors its paths, initializes Visual Studio and oneAPI, uses a
dedicated package-local CMake directory, and deploys `UniDec.exe` and
`unideclib.dll` to `unidec/bin/`. Use it instead of configuring from `src/` or
reusing an IDE cache. It expects Intel oneAPI (compiler and MKL), Ninja, and
HDF5 in locations discovered by the CMake project.

For a direct cross-platform build from the repository root:

```shell
cmake -S public/UniDec/unidec/src -B build/unidec -DCMAKE_BUILD_TYPE=Release
cmake --build build/unidec --config Release --parallel
```

Supply HDF5/FFTW roots or a toolchain as needed.

## WSL Linux Validation

Use a dedicated Linux virtual environment; never reuse the Windows project
interpreter from the mounted drive. In Ubuntu under WSL2:

```shell
sudo apt-get update
sudo apt-get install --yes build-essential ninja-build pkg-config libhdf5-dev libfftw3-dev python3-venv
python3 -m venv ~/.venvs/unidecdev
~/.venvs/unidecdev/bin/python -m pip install --upgrade pip "cmake>=3.25"
```

WSL normally appends the Windows `PATH`, which can make Linux CMake select
Windows HDF5 DLLs and fail while linking libraries such as `shlwapi`. Activate
the virtual environment, replace the inherited path with Linux locations, and
point HDF5 discovery at WSL's installation:

```shell
source ~/.venvs/unidecdev/bin/activate
export PATH="${VIRTUAL_ENV}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
cmake -S public/UniDec/unidec/src -B /tmp/unidec-build -G Ninja -DCMAKE_BUILD_TYPE=Release -DHDF5_ROOT=/usr
cmake --build /tmp/unidec-build --parallel
```

Run these in a WSL shell. For compound validation launched from PowerShell, put
the Bash commands in a temporary script rather than nesting shell metacharacters
inside `wsl.exe ... sh -lc`. Because the project deploys artifacts to the source
tree's parent `bin`, use a temporary source copy when validation must not modify
tracked package binaries.
