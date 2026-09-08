# UniDecDev Agent Guide

This file applies to the entire repository. A more deeply nested `AGENTS.md`, if one is added later, overrides these instructions for its subtree.

## Repository Map

- `public/UniDec/` is the main UniDec Python package, wxPython GUI suite, tests, and native UniDec engine.
- `public/UniDec/unidec/src/` contains the C/C++ engine and its CMake project.
- `public/UniDec/unidec/bin/` contains platform executables and shared libraries consumed by Python. Treat these as build outputs even when they are tracked.
- `public/UniDec/tests/` contains UniDec unit, GUI-launch, workflow, and native-build tests.
- `public/IsoDec/` is the independently packaged IsoDec Python/native project.
- `public/IsoDec/extern/IsoGen/` is a nested Git repository/submodule. Do not modify it unless the task explicitly includes IsoGen.
- `.github/workflows/` contains integration CI for this combined development repository. Package-specific workflows also live below each package's `.github/workflows/` directory.
- `Scripts/` contains project utilities and examples; do not treat them as the primary package APIs.

See `SKILLS.md` for a task-to-workflow map and focused validation recipes.

## Working Agreements

1. Read the nearest README, package metadata, tests, and build files before changing behavior.
2. Inspect Git status before editing. Preserve user changes and keep unrelated modifications out of the patch.
3. Make the smallest coherent change. Avoid opportunistic rewrites, mass formatting, or generated-file churn.
4. Use repository-relative paths in code and documentation. Never commit a developer-specific interpreter, dependency, or build-directory path.
5. Keep Python 3.10 compatibility for UniDec and Python 3.9 compatibility for IsoDec unless package metadata is deliberately changed.
6. Do not hand-edit native executables, shared libraries, generated build trees, caches, or IDE metadata. Rebuild artifacts from source when an explicitly requested deliverable requires them.
7. Do not change release, publishing, signing, or deployment configuration unless the request explicitly covers it.

## Cross-Layer Changes

Many UniDec features cross several boundaries:

```text
wx control -> Python config -> config/HDF5 serialization -> Python engine dispatch
           -> C command dispatch -> native implementation -> output import/display
```

When a setting or workflow crosses this path:

- Search for every read, write, default, copy, and reset of the setting.
- Keep GUI-to-config and config-to-GUI conversion symmetric.
- Keep text configuration and HDF5 names/types consistent with the native reader.
- Update Python and C dispatch together. A callable native function is not integrated until the intended command path reaches it.
- Preserve backward-compatible defaults. A new opt-in numeric mode should normally leave existing behavior active at zero and select the new path only for an explicitly valid positive value.
- Validate both the legacy branch and the new branch.
- Verify output datasets, shapes, and import code together; do not infer the native contract from only one layer.

For UniChrom and MetaUniDec work, inspect at least the relevant controller/engine, shared GUI controls, configuration serialization, `h5io.c`, native command entry point, and workflow tests.

## Python and GUI Code

- Use the interpreter configured for the package being edited; do not assume `python` resolves to a particular machine installation.
- Install UniDec development dependencies from `public/UniDec` with `python -m pip install -e ".[dev]"` when a prepared environment is not already available.
- Keep GUI construction deterministic. Tests must not depend on a user's saved default file, recent-file list, or local instrument data.
- wxPython objects generally need a real `wx.App`, the correct parent lifetime, and cleanup after tests. Use subprocess isolation for launch smoke tests when global GUI state could leak.
- Headless Linux cannot run display-dependent wx tests without a display server. Use focused non-GUI tests there or provide an appropriate display environment.
- Put regression coverage in `public/UniDec/tests/` and follow its `unittest` style for UniDec.
- Put IsoDec regression coverage in `public/IsoDec/tests/` and follow its `pytest` style.

## Native Code

- UniDec's CMake project requires CMake 3.25, HDF5 C/HL, and single-precision FFTW or MKL. OpenMP is optional.
- Keep declarations, definitions, command dispatch, and CMake source lists synchronized.
- Maintain portable C11/C++ code across Windows, Linux, macOS, x64, and ARM64. Guard compiler- or platform-specific behavior explicitly.
- Use the CMake options already provided, including `UNIDEC_ENABLE_OPENMP` and `UNIDEC_REQUIRE_INTELLLVM`, rather than hard-coding a local toolchain.
- Normal UniDec CMake builds copy the executable and shared library into `public/UniDec/unidec/bin/`. Review those artifact changes separately from source changes.
- Treat warnings involving array sizes, HDF5 types, signedness, and narrowing as possible data-corruption defects, not cosmetic issues.
- For numerical changes, test edge cases such as empty inputs, singleton dimensions, invalid parameter values, and non-finite data.

## Validation

Run the narrowest relevant checks first, then broaden according to risk.

From `public/UniDec`:

```shell
python -m unittest tests.test_unichrom_workflow -v
python -m unittest tests.test_gui_launch -v
python -m unittest discover -s tests -p "test_*.py" -v
```

On Linux or macOS, build the packaged UniDec native executable with:

```shell
bash ./build_native.sh
```

For a direct cross-platform native build, configure a separate build directory and supply HDF5/FFTW roots or a toolchain as needed:

```shell
cmake -S public/UniDec/unidec/src -B build/unidec -DCMAKE_BUILD_TYPE=Release
cmake --build build/unidec --config Release --parallel
```

From `public/IsoDec`:

```shell
git submodule update --init --recursive
python -m pip install -e ".[test]"
python -m pytest
```

Do not install dependencies, initialize submodules, or perform a full native rebuild merely to validate a documentation-only or isolated Python change. If required dependencies, GUI display support, fixtures, or a compiler toolchain are unavailable, run the checks that remain meaningful and report the exact limitation.

## Review Checklist

- Behavior matches the request, including default and error paths.
- Public or serialized names remain compatible unless migration is part of the task.
- GUI, Python, HDF5/config, native dispatch, and output-import layers agree.
- Focused regression tests cover all affected branches.
- No unrelated source, binary, cache, IDE, or build-tree changes are included.
- Documentation and command examples match the repository's current files and CI.
