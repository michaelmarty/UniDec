# Sessions and MCP

**Tier 1.** Sessions provide a headless API for callers that need a dedicated workspace and structured, bounded results. Pass an explicit `working_directory` for outputs and `allowed_input_roots` for source files. Inspect `ok`, `error`, and state availability after each operation. Native and importer requirements depend on the [mode](modes.md).

## Single spectra

Import `UniDecSession` from `unidec.session`. It defaults to conventional MS. For IM-MS, pass an instance of `unidec.modules.IMEng.UniDecIM` as `engine`. For CD-MS, pass an instance of `unidec.modules.CDEng.UniDecCD` and an `engine_adapter` constructed with `unidec.session.UniDecCDSessionEngineAdapter(engine)` using that same instance.

Call `open_spectrum`, inspect or apply settings with `get_parameters` and `set_parameters`, then call `run("all")`. Use `get_results` for peak summaries and `get_plot` to write a managed plot. `save_state` creates a ZIP checkpoint. `load_state` restores a checkpoint from an allowed location.

The [conventional-MS example](../guides/conventional-ms.md) shows the full session lifecycle from opening data through rerunning a restored checkpoint. The [IM-MS example](../guides/im-ms.md) inspects mobility and CCS grids. The [CD-MS example](../guides/cd-ms.md) inspects the charge histogram and fit.

## MetaUniDec

Import `MetaUniDecCollectionSession` from `unidec.meta_session`. Call `open_collection` to copy a source HDF5 project into the workspace, then call `run` and `get_results`. Use `save_checkpoint` and `load_checkpoint` for HDF5 snapshots.

The [MetaUniDec workflow](../guides/metaunidec.md) builds a collection, inspects its mass axis and peak extracts, and reruns a restored checkpoint.

(unichrom)=
## UniChrom

Import `UniChromSession` from `unidec.unichrom_session`. Call `open_chromatogram` to create the managed HDF5 project from a raw source, then call `create_time_windows`, `run`, and `get_results`.

The [UniChrom workflow](../guides/unichrom.md) imports a small mzML, creates two time windows, and inspects the generated HDF5 project.

## MCP clients

Install the `agent` extra and launch `unidec-mcp` with the acquisition mode, input roots, and working directory. The [MCP guide](../guides/agent-interface.md) lists tool names, the supported plotting and checkpoint operations, and the optional live GUI bridge.
