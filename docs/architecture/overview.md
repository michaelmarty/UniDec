# UniDec internal architecture

This contributor guide traces the Python and C implementation. The [Python API map](../api/index.md) is the starting point for choosing a supported workflow. The [internal module reference](../api/internals/index.md) indexes implementation symbols.

```text
wx launcher and view -> mode presenter -> Python engine
Python session -----------------------> Python engine
                                            |
                              configuration and prepared data
                                            |
                               UniDec executable in unidec/bin
                                            |
                         output files or updated HDF5 datasets
                                            |
                          Python result imports -> caller or view
```

- **Entry points.** The launcher selects the presenter and the presenter owns the engine and view. The view exchanges settings with the engine through the presenter. Python scripts and sessions can create engines directly. Sessions add path handling, managed workspaces, stage freshness, and bounded results. See the [GUI extension guide](../guides/gui-extensions.md).
- **Python engines.** `unidec.engine.UniDec` owns the conventional workflow. `modules.IMEng.UniDecIM` specializes mobility input, preprocessing, and result import. `modules.CDEng.UniDecCD` builds a charge-by-m/z histogram from events before deconvolution. `metaunidec.mudeng.MetaUniDec` manages projects, with `modules.ChromEng.ChromEngine` deriving project spectra from chromatogram windows.
- **Native implementation.** `unidec/src/UniDec.c` loads the configuration and chooses the command path. Mode-specific implementations read the inputs and write outputs. The Python importer translates those outputs into arrays and peak objects suited to the mode.
- **Build products.** `unidec/src/CMakeLists.txt` builds the executable and `unideclib` shared library and deploys them to `unidec/bin`. The normal workflows below run the executable in a subprocess. The separate `unidec.modules.unidecwrapper` module loads the shared library through ctypes. Its ABI and memory ownership are a separate integration surface.

```{toctree}
:maxdepth: 1

engine-state
native-boundary
```
