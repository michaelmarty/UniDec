# UniDec Package Guide

This file applies to `public/UniDec/`, the UniDec Python package, wxPython GUI
suite, tests, public scripts, and package build/release files. The native engine
has additional instructions in `unidec/AGENTS.md`.

## Package Boundaries

- Keep Python 3.10 compatibility unless package metadata is deliberately
  changed.
- `tests/` contains unit, GUI-launch, workflow, and native-build tests. Put
  UniDec regression coverage there and follow its `unittest` style.
- `PublicScripts/` may import reusable package APIs, but package code must not
  import from or depend on `PublicScripts/` or the repository-level `Scripts/`
  tree. Tests importing either script tree are script tests and must remain
  outside automatic public-package testing.
- Use the package's configured interpreter. Install development dependencies
  from this directory with `python -m pip install -e ".[dev]"` only when a
  prepared environment is unavailable. For `PublicScripts/` work, use
  `python -m pip install -e ".[scripts]"`, or `.[dev,scripts]` when both are
  needed; specialized scripts may still require vendor runtimes.

## Cross-Layer Changes

Many settings and workflows cross this path:

```text
wx control -> Python config -> config/HDF5 serialization -> Python engine dispatch
           -> C command dispatch -> native implementation -> output import/display
```

- Search for every read, write, default, copy, and reset of a changed setting.
- Keep GUI-to-config and config-to-GUI conversion symmetric.
- Keep text configuration and HDF5 names/types consistent with the native
  reader; update Python and C dispatch together.
- Preserve backward-compatible defaults. A new opt-in numeric mode should
  normally keep existing behavior at zero and activate only for an explicitly
  valid positive value.
- Validate legacy and new branches. Check output datasets, shapes, and import
  code together rather than inferring the native contract from one layer.
- For UniChrom and MetaUniDec, inspect the relevant controller/engine, shared
  GUI controls, configuration serialization, `h5io.c`, native command entry
  point, and workflow tests.

## Python and GUI Tests

- Keep GUI construction deterministic; tests must not depend on saved defaults,
  recent files, or local instrument data.
- wxPython objects generally require a real `wx.App`, correct parent lifetime,
  and cleanup. Use subprocess isolation for launch smoke tests when global GUI
  state could leak.
- Headless Linux cannot run display-dependent wx tests without a display server.
  Use focused non-GUI tests or provide a suitable display environment.
- `tests/run_headless_tests.py` is the authoritative headless suite. It blocks
  wxPython and GUI-package imports even when installed. Keep GUI-independent
  calculations, configuration, serialization, and engine behavior outside GUI
  entry points and `gui_elements`; never weaken the guard for misplaced code.
- `tests/run_automatic_tests.py` is the broader GUI-capable suite and excludes
  tests targeting `Scripts/` and `PublicScripts/`.

Run the narrowest relevant checks first from `public/UniDec`:

```shell
python -m unittest tests.test_unichrom_workflow -v
python tests/run_headless_tests.py
python -m unittest tests.test_gui_launch -v
python tests/run_automatic_tests.py
```

On managed Windows hosts, Python may create a temporary directory but later
deny NumPy, HDF5, or ordinary writes with `PermissionError`; redirecting `TEMP`
and `TMP` into the repository may not help. When the traceback clearly matches
this sandbox restriction, rerun the unchanged test with sandbox escalation.
Do not alter application code, tests, or temp paths as a workaround, and report
both the original limitation and escalated result.
