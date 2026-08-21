# Building and publishing IsoDec

IsoDec uses scikit-build-core. Every wheel compiles `isodeclib` and the pinned
IsoGen submodule, then installs both native libraries beside one another in
`isodec/bin`.

## Prepare a release

1. Update `isodec/_version.py` and `CITATION.cff`.
2. Update release notes and verify the intended IsoGen submodule commit.
3. Run `git submodule update --init --recursive` and the full test suite.
4. Review `THIRD_PARTY_NOTICES.md` and the bundled license files.
5. Push the release commit and run the **Build and publish** workflow.

The workflow builds Windows, Linux, macOS Intel, and macOS Apple Silicon
wheels plus a source distribution. A manual run creates the matching GitHub
release and can publish through PyPI trusted publishing.

## Local build

```shell
git submodule update --init --recursive
python -m pip install --upgrade build twine
python -m build
python -m twine check dist/*
```

Source builds require CMake, a C/C++ compiler, and FFTW. Windows uses the FFTW
import library in IsoGen. Linux distributions typically provide FFTW through
`libfftw3-dev`; on macOS use `brew install fftw`.

Test the wheel in a clean environment from outside the repository:

```shell
python -m venv wheel-test
wheel-test/Scripts/python -m pip install dist/isodec-*.whl
wheel-test/Scripts/python -m pytest tests --import-mode=importlib
```

Use `wheel-test/bin/python` on Linux and macOS. Linux wheels are repaired with
auditwheel and macOS wheels with delocate so FFTW is bundled.
