# UniDec MCP Agent Interface

The optional MCP adapter exposes one headless conventional-MS, IM-MS, CD-MS, MetaUniDec, or UniChrom session through a bounded set of tools. An optional local bridge can also route workflow commands to a separately running conventional UniDec GUI. It does not provide shell access, arbitrary Python or filesystem operations, native command selection, DoubleDec, or SMAS.

## Install and start

From `public/UniDec`, install UniDec and the optional MCP dependency into the Python environment that Codex will use:

```shell
python -m pip install -e ".[agent]"
```

The adapter currently targets the MCP 1.x `FastMCP` API, so the optional dependency is constrained to `mcp>=1.9,<2`. If MCP 2.x was installed previously, rerun the command above after updating the checkout, or downgrade it directly:

```shell
python -m pip install --upgrade "mcp>=1.9,<2"
```

Verify the selected API before registering the server:

```shell
python -c "from mcp.server.fastmcp import FastMCP; print('MCP 1.x FastMCP import OK')"
```

The server requires at least one readable input root and a separate working directory. Repeat `--input-root` when spectra or checkpoints live under multiple directories.

```shell
python -m unidec.mcp_server --input-root C:/Data/UniDec --working-directory C:/Data/UniDecAgent
```

Conventional MS is the default. Start a separate IM-MS server with `--mode imms`; one server owns exactly one mode and session:

```shell
python -m unidec.mcp_server --mode imms --input-root C:/Data/UniDec --working-directory C:/Data/UniDecIMAgent
```

IM-MS exposes the reviewed common parameters plus arrival-time bounds, arrival-time smoothing/background subtraction, arrival-time peak width, and CCS bounds/bin size. Instrument calibration controls remain outside the public parameter registry. In IM-MS mode, `get_plot` additionally accepts `ccs_distribution`; `processed_spectrum` is the m/z projection of the processed 2D grid.

Start CD-MS with `--mode cdms`. CD-MS adds reviewed charge-calibration slope, charge binning/resolution, charge peak width/smoothing, and calibration subtype parameters. Its event cache and all processing/deconvolution outputs remain below the managed working directory. CD-MS checkpoints restore through the versioned `rawdata.npz` event cache and intentionally mark computed results stale until rerun.

Start a separate MetaUniDec collection server with `--mode metaunidec`:

```shell
python -m unidec.mcp_server --mode metaunidec --input-root C:/Data/UniDec --working-directory C:/Data/MetaUniDecAgent
```

MetaUniDec uses a separate eleven-tool collection contract described below. Opening a project first copies its HDF5 file into a unique operation directory, and all configuration writes and native execution use that managed copy. The source project is retained only as provenance and is not processed in place.

Start UniChrom with `--mode unichrom`:

```shell
python -m unidec.mcp_server --mode unichrom --input-root C:/Data/Chromatograms --working-directory C:/Data/UniChromAgent
```

UniChrom accepts a raw chromatogram file or vendor directory below an input root. Its derived HDF5 project and native outputs are routed into managed storage; no adjacent project or sidecar is created beside the raw source. The initial contract supports regular time windows, not regular scan windows, detected chromatographic-peak windows, or manual selections.

The installed `unidec-mcp` entry point accepts the same arguments.

The server uses STDIO for MCP messages. Engine progress is sent to stderr so it does not corrupt the protocol stream. Do not use the command as an interactive console program.

## Configure Codex

Codex supports local STDIO MCP servers, and its CLI and IDE extension share MCP configuration. The default user configuration is `~/.codex/config.toml`; a trusted project can instead use `.codex/config.toml`. See the [official Codex MCP documentation](https://learn.chatgpt.com/docs/extend/mcp) for the current configuration and approval options.

Add the server from a terminal:

```shell
codex mcp add unidec -- python -m unidec.mcp_server --input-root C:/Data/UniDec --working-directory C:/Data/UniDecAgent
codex mcp list
```

Alternatively, configure it directly. Use paths valid on the host; forward slashes make Windows TOML examples easier to read.

```toml
[mcp_servers.unidec]
command = "python"
args = [
  "-m", "unidec.mcp_server",
  "--input-root", "C:/Data/UniDec",
  "--working-directory", "C:/Data/UniDecAgent",
]
cwd = "<path-to-UniDecDev>/public/UniDec"
startup_timeout_sec = 30
tool_timeout_sec = 600
default_tools_approval_mode = "prompt"
enabled = true
```

Restart the IDE extension after changing its MCP configuration. Keep approval mode at `prompt` until the local roots and expected tool behavior have been reviewed. A host timeout can stop waiting for a tool response, but it cannot safely interrupt a native UniDec calculation already in progress.

## Connect a live GUI

Choose a bridge descriptor below the MCP working directory, then set it before starting the conventional GUI from PowerShell:

```powershell
$env:UNIDEC_GUI_BRIDGE_FILE = "C:/Data/UniDecAgent/gui-bridge.json"
gunidec --unidec
```

Add the same path to the MCP server arguments:

```toml
args = [
  "-m", "unidec.mcp_server",
  "--input-root", "C:/Data/UniDec",
  "--working-directory", "C:/Data/UniDecAgent",
  "--gui-bridge-file", "C:/Data/UniDecAgent/gui-bridge.json",
]
```

With this option configured, the server exposes three additional tools:

1. `get_gui_parameters(names)` reads values from the controls currently visible in the live GUI and returns the conventional-MS parameter schema.
2. `set_gui_parameters(values)` validates the complete visible public parameter set and the proposed update on an isolated configuration copy. It commits only a valid result and then refreshes the visible controls.
3. `run_gui(stage)` accepts `process`, `deconvolve`, `pick_peaks`, or `all` and routes the workflow through the normal GUI presenter.

All three tools are serialized and execute GUI access on the wx main thread. They report `busy` for an overlapping agent command and `gui_unavailable` when the descriptor is missing, stale, or unreachable. Parameter responses include a GUI revision that increments after each successful agent update.

The bridge is bound to `127.0.0.1` and authenticates requests with a random token stored in the descriptor. The GUI removes its descriptor when it closes. The MCP process rereads the descriptor for every call, so the GUI can be closed and restarted without restarting Codex. The live GUI and the twelve single-spectrum headless tools are separate sessions; `run_gui` acts only on the file and settings already visible in the GUI.

## Tool workflow

Use the tools in this order:

1. `open_spectrum(path)` opens a file below a configured input root.
2. `get_state()` and `get_parameters(names)` inspect the current revision and supported parameter schema.
3. `set_parameters(values)` validates and applies a parameter update atomically.
4. `run(stage)` executes `process`, `deconvolve`, `pick_peaks`, or `all`.
5. `get_results(max_peaks)` returns bounded structured results; `get_plot(name)` writes a managed PNG artifact.
6. `get_run_manifest(operation_id)` returns the latest run record when the ID is omitted, or one retained run by operation ID. `compare_runs(before_operation_id, after_operation_id)` reports parameter changes and deltas for R-squared, peak count, and dominant peak mass.
7. `export_regression_fixture(operation_id, path)` writes a versioned JSON description for a retained failed run. The destination must be a new `.json` file below the working directory.
8. `save_state(path)` writes a new ZIP below the working directory; `load_state(path)` restores a ZIP from an input root or the working directory.

Every mutation has an operation ID and invalidates stale downstream results where appropriate. Mutations are serialized and can return `busy`. Check `ok`, `error`, `state.revision`, and `state.available` before continuing. `save_state` refuses to overwrite an existing checkpoint.

Successful and failed run attempts retain a bounded in-memory manifest. Each manifest records the UniDec version, mode, input identity, public parameters and changes from the preceding run, workload counts, requested stage, UTC start time, duration, status, error, and bounded result metrics. Single-spectrum and MetaUniDec files up to 256 MiB also receive a SHA-256 hash; larger files explicitly report that hashing was omitted to avoid delaying normal data loading. A session retains the latest 100 manifests and at most 25 height-ranked peaks per manifest. These records are session-local and are not restored from a checkpoint.

Regression-fixture descriptions contain only a schema identifier, version and mode, source run/revision, input identity, public parameters, bounded workload counts, requested stage, and expected error. Error messages are limited to 2,000 characters. The export does not copy the input, include internal output paths or session identifiers, or generate executable test code. Successful runs are rejected, and existing files are never overwritten.

### MetaUniDec collection workflow

The `metaunidec` mode exposes collection-specific names rather than treating an HDF5 project as one spectrum:

1. `open_collection(path)` copies and opens an HDF5 project from an allowed input root.
2. `get_collection_state()` returns at most 100 spectrum summaries plus aggregate readiness; it does not return complete per-spectrum arrays.
3. `get_parameters(names)` and `set_parameters(values)` read or atomically update reviewed collection-wide parameters. Spectrum variables and inclusion are not parameter aliases and are not exposed yet.
4. `run_collection()` performs preprocessing, native deconvolution, aggregate grid generation, and peak extraction.
5. `get_collection_results(max_peaks)` returns at most 100 height-ranked aggregate peaks and bounded mass metadata.
6. `get_collection_run_manifest(operation_id)` and `compare_collection_runs(before_operation_id, after_operation_id)` return bounded run provenance and parameter/aggregate metric deltas. `export_collection_regression_fixture(operation_id, path)` writes a new bounded JSON description for a retained failed run.
7. `save_collection_checkpoint(path)` writes a new HDF5 snapshot below the working directory; `load_collection_checkpoint(path)` opens a permitted snapshot through the same managed-copy boundary. Loaded results are stale until `run_collection()` succeeds.

MetaUniDec checkpoints never overwrite an existing file. To load a checkpoint created under a different server working directory, that directory must also be configured as an input root.

### UniChrom workflow

The `unichrom` mode exposes ten tools:

1. `open_chromatogram(path)` opens a permitted raw file or vendor directory and creates an empty managed HDF5 project.
2. `get_chromatogram_state()` returns bounded TIC/BPC ranges and counts, scan metadata, window count, and result freshness; complete chromatogram arrays are not returned.
3. `get_parameters(names)` and `set_parameters(values)` expose reviewed common deconvolution settings plus chromatogram selection, time/scan window widths, optional time bounds, chromatographic peak width, and coupled-deconvolution width.
4. `create_time_windows(width, start, end)` replaces the managed collection with regular time windows. `start` and `end` must either both be omitted or form an increasing finite range.
5. `run_unichrom()` preprocesses and deconvolves the managed windows and performs aggregate peak extraction.
6. `get_unichrom_results(max_peaks)` returns at most 100 height-ranked aggregate peaks, mass bounds, window count, and UniScore when available.
7. `get_unichrom_run_manifest(operation_id)` and `compare_unichrom_runs(before_operation_id, after_operation_id)` return bounded run provenance and parameter/aggregate metric deltas. `export_unichrom_regression_fixture(operation_id, path)` writes a new bounded JSON description for a retained failed run.

UniChrom checkpoint tools are intentionally absent. The current derived HDF5 project is not a portable raw-chromatogram checkpoint because regenerating windows still requires the original importer-supported source.

## Security boundary

- Input paths are resolved and must remain below an explicitly configured root.
- The server mode is fixed at startup as `ms`, `imms`, `cdms`, `metaunidec`, or `unichrom`; changing modes is not a tool operation.
- Plots, operation files, and saved checkpoints remain below the working directory.
- Checkpoint archives with absolute or parent-traversal members are rejected before extraction.
- The adapter declares twelve single-spectrum tools, eleven MetaUniDec collection tools, or ten UniChrom tools and, only for the single-spectrum server when configured, three allowlisted GUI tools; it has no generic command, import, or file tool.
- Run the server with a dedicated working directory and only the minimum input roots needed for the task.
