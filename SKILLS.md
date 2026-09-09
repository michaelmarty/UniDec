# UniDecDev Workflow Guide

This companion to `AGENTS.md` routes common tasks to the files and checks most likely to matter. It is repository documentation, not an installable Codex skill package.

## Agent Principles

Apply these principles when using any workflow in this guide:

- Think before coding. State assumptions, surface ambiguity and tradeoffs, and ask for clarification when they materially affect the result.
- Choose the simplest complete implementation. Do not add speculative features, one-use abstractions, or unrequested flexibility.
- Keep changes surgical. Match local style, avoid unrelated cleanup, and remove only orphans created by the current change.
- Define success in terms of checks. For multi-step work, state a short plan with a verification step for each part, then run the narrowest relevant checks and broaden them according to risk.

## Choose the Workflow

| Task | Inspect first | Minimum useful validation |
|---|---|---|
| UniDec Python engine behavior | `public/UniDec/unidec/engine.py`, relevant module, nearby tests | Focused `unittest` module |
| UniDec/UniChrom GUI behavior | launcher/controller, engine, `modules/gui_elements/` or `metaunidec/gui_elements/` | Focused workflow test plus `tests.test_gui_launch` when construction changes |
| Configuration option | config default, import/export, GUI bindings, HDF5/text serialization, native reader | Tests for default, round trip, and each dispatch branch |
| UniDec native command or algorithm | `unidec/src/CMakeLists.txt`, command entry point, header/implementation, Python caller | Native build plus focused Python/native workflow test |
| HDF5 input/output contract | Python HDF5 writer/reader, `unidec/src/h5io.c`, consumers of each dataset | Temporary-file round trip and shape/type assertions |
| IsoDec Python behavior | `public/IsoDec/isodec/`, package tests, `pyproject.toml` | Focused `pytest`, then full `python -m pytest` |
| IsoDec native behavior | `public/IsoDec/isodec/src/`, bindings, IsoGen boundary | Build/install package and run native API tests |
| Packaging or CI | relevant `pyproject.toml`, package workflow, root integration workflow | Reproduce the changed workflow command locally where practical |

## Focused UniDec Test Recipes

Run these from `public/UniDec`:

```shell
# GUI-independent UniChrom routing and configuration behavior
python -m unittest tests.test_unichrom_workflow -v

# Major window construction and control availability; requires GUI support
python -m unittest tests.test_gui_launch -v

# Standard and MetaUniDec workflow behavior
python -m unittest tests.test_unidec_workflows tests.test_metaunidec_workflows -v

# Ion-mobility engine separation
python -m unittest tests.test_im_engine_split -v

# Native build/command smoke coverage
python -m unittest tests.test_native_build -v

# Full local UniDec suite
python -m unittest discover -s tests -p "test_*.py" -v
```

Some GUI/workflow tests use importer data. CI provides it through `UNIDEC_IMPORTER_TEST_DATA` and can require it with `UNIDEC_REQUIRE_IMPORTER_TEST_DATA=1`. Do not silently substitute unrelated local data when reproducing those tests.

## Configuration-to-Native Checklist

Use this checklist when adding or changing a field such as a processing width, threshold, or mode selector:

1. Define a safe Python default.
2. Import the current value into every relevant GUI.
3. Export and validate user input back into the config.
4. Persist and reload the value through each supported config format.
5. Write the HDF5 key with the type and location expected by C.
6. Read it in native code with an explicit fallback for older files.
7. Dispatch only when its activation condition is met.
8. Retain the previous processing path when the option is disabled.
9. Import native outputs and confirm expected dimensions before plotting.
10. Test enabled, disabled/default, malformed, and missing-value cases.

## Native Build Selection

- Linux/macOS UniDec: prefer `public/UniDec/build_native.sh`; it handles package deployment and platform output names.
- Windows UniDec: mirror `.github/workflows/test-unidec-gui-workflows.yml`. Use an out-of-tree CMake build with HDF5 and FFTW supplied through vcpkg or explicit roots.
- Windows IntelLLVM: enable `UNIDEC_REQUIRE_INTELLLVM=ON` when the compiler choice is part of the requirement.
- ARM64 without OpenMP support: configure with `UNIDEC_ENABLE_OPENMP=OFF`.
- IsoDec: build through its Python package (`pip install` or `python -m build`) so scikit-build-core exercises the same integration used for wheels.

Never reuse an unknown or stale CMake cache to judge a source change. Prefer a clearly named, task-specific build directory outside the source tree.

## Completion Notes

When handing work back, state:

- what behavior changed;
- which source and contract boundaries were updated;
- which focused and broad checks passed;
- which checks were skipped and why;
- whether tracked native artifacts changed.
