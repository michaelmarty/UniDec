# Development

Clone with submodules and install an editable build:

```shell
git clone --recurse-submodules https://github.com/michaelmarty/IsoDec.git
cd IsoDec
python -m pip install -e ".[test]"
python -m pytest
```

The CMake project is rooted at `isodec/src`. It adds
`extern/IsoGen/src` as a subdirectory and links `isodeclib` to the `isogen`
target. Do not copy IsoGen source or prebuilt native libraries into `isodec`.

Build documentation with:

```shell
python -m pip install -e ".[docs]"
python -m mkdocs build --strict
```

See `PUBLISHING.md` for platform wheel and release procedures.
