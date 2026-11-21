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

## Supported Python reference

| Entry point | Supported operations |
| --- | --- |
| {any}`UniDecSession <unidec.session.UniDecSession>` | {any}`open_spectrum <unidec.session.UniDecSession.open_spectrum>`, {any}`get_state <unidec.session.UniDecSession.get_state>`, {any}`get_parameters <unidec.session.UniDecSession.get_parameters>`, {any}`set_parameters <unidec.session.UniDecSession.set_parameters>`, {any}`run <unidec.session.UniDecSession.run>`, {any}`get_results <unidec.session.UniDecSession.get_results>`, {any}`get_plot <unidec.session.UniDecSession.get_plot>`, {any}`save_state <unidec.session.UniDecSession.save_state>`, {any}`load_state <unidec.session.UniDecSession.load_state>`, {any}`get_run_manifest <unidec.session.UniDecSession.get_run_manifest>`, {any}`compare_runs <unidec.session.UniDecSession.compare_runs>`, {any}`export_regression_fixture <unidec.session.UniDecSession.export_regression_fixture>` |
| {any}`UniDecCDSessionEngineAdapter <unidec.session.UniDecCDSessionEngineAdapter>` | Use with a CD-MS engine when constructing a single-spectrum session. |
| {any}`MetaUniDecCollectionSession <unidec.meta_session.MetaUniDecCollectionSession>` | {any}`open_collection <unidec.meta_session.MetaUniDecCollectionSession.open_collection>`, {any}`get_state <unidec.meta_session.MetaUniDecCollectionSession.get_state>`, {any}`get_parameters <unidec.meta_session.MetaUniDecCollectionSession.get_parameters>`, {any}`set_parameters <unidec.meta_session.MetaUniDecCollectionSession.set_parameters>`, {any}`run <unidec.meta_session.MetaUniDecCollectionSession.run>`, {any}`get_results <unidec.meta_session.MetaUniDecCollectionSession.get_results>`, {any}`save_checkpoint <unidec.meta_session.MetaUniDecCollectionSession.save_checkpoint>`, {any}`load_checkpoint <unidec.meta_session.MetaUniDecCollectionSession.load_checkpoint>`, {any}`get_run_manifest <unidec.meta_session.MetaUniDecCollectionSession.get_run_manifest>`, {any}`compare_runs <unidec.meta_session.MetaUniDecCollectionSession.compare_runs>`, {any}`export_regression_fixture <unidec.meta_session.MetaUniDecCollectionSession.export_regression_fixture>` |
| {any}`UniChromSession <unidec.unichrom_session.UniChromSession>` | {any}`open_chromatogram <unidec.unichrom_session.UniChromSession.open_chromatogram>`, {any}`get_state <unidec.unichrom_session.UniChromSession.get_state>`, {any}`get_parameters <unidec.unichrom_session.UniChromSession.get_parameters>`, {any}`set_parameters <unidec.unichrom_session.UniChromSession.set_parameters>`, {any}`create_time_windows <unidec.unichrom_session.UniChromSession.create_time_windows>`, {any}`run <unidec.unichrom_session.UniChromSession.run>`, {any}`get_results <unidec.unichrom_session.UniChromSession.get_results>`, {any}`get_run_manifest <unidec.unichrom_session.UniChromSession.get_run_manifest>`, {any}`compare_runs <unidec.unichrom_session.UniChromSession.compare_runs>`, {any}`export_regression_fixture <unidec.unichrom_session.UniChromSession.export_regression_fixture>` |
| {any}`unidec-mcp entry point <unidec.mcp_server.main>` | Starts the STDIO server using the supplied command-line settings. |
| {any}`create_server <unidec.mcp_server.create_server>` | Registers the tool set on an MCP server for a supplied session. |
