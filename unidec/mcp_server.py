"""Optional MCP adapter for the headless :class:`UniDecSession` API."""

import argparse
from contextlib import redirect_stdout
import sys

from unidec.session import UniDecCDSessionEngineAdapter, UniDecSession
from unidec.gui_agent_bridge import GuiAgentBridgeClient


SERVER_INSTRUCTIONS = """Operate only through the declared UniDec tools. Open inputs only below an administrator-configured input root; all generated artifacts and checkpoints stay below the managed working directory. Never invent paths or assume stale results are current. Use this sequence: open_spectrum, inspect or set parameters, run, then get_results or get_plot. Use get_run_manifest and compare_runs for bounded provenance and objective before/after summaries. export_regression_fixture may write a non-overwriting JSON description for a retained failed run. Mutations are serialized and may return busy. Check every response's ok/error and state fields before continuing. save_state never overwrites an existing file. If the GUI tools are available, they control a separate live conventional-MS GUI session and do not share state with the headless tools. Inspect that session with get_gui_parameters, update it with set_gui_parameters, and run it with run_gui. The headless server exposes the configured conventional MS, IM-MS, or CD-MS mode only; it does not expose a shell, arbitrary Python, generic file access, native command selection, MetaUniDec, UniChrom, DoubleDec, or SMAS."""

COLLECTION_SERVER_INSTRUCTIONS = """Operate only through the declared MetaUniDec collection tools. Open HDF5 projects only below an administrator-configured input root. The server copies every source project into its managed working directory before opening or changing it; the source remains provenance and is never processed in place. Use this sequence: open_collection, inspect or atomically set collection-wide parameters, run_collection, then get_collection_results. Use get_collection_run_manifest and compare_collection_runs for bounded provenance and aggregate before/after summaries. export_collection_regression_fixture may write a non-overwriting JSON description for a retained failed run. Mutations are serialized and may return busy. Check every response's ok/error and state fields before continuing. Collection checkpoints are non-overwriting HDF5 snapshots below the managed working directory; a loaded checkpoint is stale until run_collection succeeds. Spectrum variables and inclusion are not parameters and are not exposed in this initial MCP contract. The server does not expose a shell, arbitrary Python, generic file access, native command selection, per-spectrum editing, or live MetaUniDec GUI control."""

UNICHROM_SERVER_INSTRUCTIONS = """Operate only through the declared UniChrom tools. Open raw chromatogram files or vendor directories only below an administrator-configured input root. Derived HDF5 projects and native outputs stay below the managed working directory; the source is never processed in place. Use this sequence: open_chromatogram, inspect or atomically set reviewed parameters, create_time_windows, run_unichrom, then get_unichrom_results. Use get_unichrom_run_manifest and compare_unichrom_runs for bounded provenance and aggregate before/after summaries. export_unichrom_regression_fixture may write a non-overwriting JSON description for a retained failed run. Mutations are serialized and may return busy. Check every response's ok/error and state fields before continuing. Only regular time-window creation is exposed initially; scan windows, detected-peak windows, manual windows, checkpoints, and live GUI control are not available. The server does not expose a shell, arbitrary Python, generic file access, or native command selection."""


def _session_call(method, *args, **kwargs):
    """Keep legacy engine progress messages off the MCP STDIO protocol stream."""
    with redirect_stdout(sys.stderr):
        return method(*args, **kwargs)


def register_tools(server, session, annotations_factory, gui_bridge=None):
    """Register the stable tool contract; separated for SDK-free contract tests."""

    def tool(name, description, *, read_only=False, idempotent=False):
        annotations = annotations_factory(
            readOnlyHint=read_only,
            destructiveHint=False,
            idempotentHint=idempotent,
            openWorldHint=False,
        )
        return server.tool(name=name, description=description, annotations=annotations)

    @tool("open_spectrum", "Open one spectrum from a configured input root.")
    def open_spectrum(path: str) -> dict:
        return _session_call(session.open_spectrum, path)

    @tool("get_state", "Return session identity, status, revision, and result freshness.", read_only=True,
          idempotent=True)
    def get_state() -> dict:
        return _session_call(session.get_state)

    @tool("get_parameters", "Return parameter values and their machine-readable schema.", read_only=True,
          idempotent=True)
    def get_parameters(names: list[str] | None = None) -> dict:
        return _session_call(session.get_parameters, names)

    @tool("set_parameters", "Atomically validate and update supported UniDec parameters.")
    def set_parameters(values: dict) -> dict:
        return _session_call(session.set_parameters, values)

    @tool("run", "Run process, deconvolve, pick_peaks, or all for the current revision.")
    def run(stage: str = "all") -> dict:
        return _session_call(session.run, stage)

    @tool("get_results", "Return bounded structured peaks from the current completed revision.", read_only=True,
          idempotent=True)
    def get_results(max_peaks: int = 100) -> dict:
        return _session_call(session.get_results, max_peaks)

    @tool("get_plot", "Render a named plot to a PNG inside the managed working directory.", idempotent=True)
    def get_plot(name: str) -> dict:
        return _session_call(session.get_plot, name)

    @tool("get_run_manifest", "Return bounded provenance for one run operation.", read_only=True,
          idempotent=True)
    def get_run_manifest(operation_id: int | None = None) -> dict:
        return _session_call(session.get_run_manifest, operation_id)

    @tool("compare_runs", "Compare parameters and summary metrics for two recorded runs.", read_only=True,
          idempotent=True)
    def compare_runs(before_operation_id: int, after_operation_id: int) -> dict:
        return _session_call(session.compare_runs, before_operation_id, after_operation_id)

    @tool("export_regression_fixture", "Write a bounded JSON description for one retained failed run.")
    def export_regression_fixture(operation_id: int, path: str) -> dict:
        return _session_call(session.export_regression_fixture, operation_id, path)

    @tool("save_state", "Save a new ZIP checkpoint inside the managed working directory.")
    def save_state(path: str) -> dict:
        return _session_call(session.save_state, path)

    @tool("load_state", "Load a ZIP checkpoint from a configured input root or the working directory.")
    def load_state(path: str) -> dict:
        return _session_call(session.load_state, path)

    if gui_bridge is not None:
        @tool("get_gui_parameters", "Return visible parameter values and schema from the separate live GUI.",
              read_only=True, idempotent=True)
        def get_gui_parameters(names: list[str] | None = None) -> dict:
            return _session_call(gui_bridge.get_parameters, names)

        @tool("set_gui_parameters", "Atomically validate and update parameters in the separate live GUI.")
        def set_gui_parameters(values: dict) -> dict:
            return _session_call(gui_bridge.set_parameters, values)

        @tool("run_gui", "Run process, deconvolve, pick_peaks, or all in the separate live UniDec GUI.")
        def run_gui(stage: str = "all") -> dict:
            return _session_call(gui_bridge.run, stage)

    return server


def register_collection_tools(server, session, annotations_factory):
    """Register the bounded MetaUniDec collection contract."""

    def tool(name, description, *, read_only=False, idempotent=False):
        annotations = annotations_factory(
            readOnlyHint=read_only,
            destructiveHint=False,
            idempotentHint=idempotent,
            openWorldHint=False,
        )
        return server.tool(name=name, description=description, annotations=annotations)

    @tool("open_collection", "Copy and open one MetaUniDec HDF5 project from a configured input root.")
    def open_collection(path: str) -> dict:
        return _session_call(session.open_collection, path)

    @tool("get_collection_state", "Return bounded collection identity, spectra, status, and freshness.",
          read_only=True, idempotent=True)
    def get_collection_state() -> dict:
        return _session_call(session.get_state)

    @tool("get_parameters", "Return collection-wide parameter values and schema.", read_only=True,
          idempotent=True)
    def get_parameters(names: list[str] | None = None) -> dict:
        return _session_call(session.get_parameters, names)

    @tool("set_parameters", "Atomically validate and update collection-wide UniDec parameters.")
    def set_parameters(values: dict) -> dict:
        return _session_call(session.set_parameters, values)

    @tool("run_collection", "Process, deconvolve, and pick aggregate peaks for the collection.")
    def run_collection() -> dict:
        return _session_call(session.run, "all")

    @tool("get_collection_results", "Return bounded aggregate peaks for the current revision.",
          read_only=True, idempotent=True)
    def get_collection_results(max_peaks: int = 100) -> dict:
        return _session_call(session.get_results, max_peaks)

    @tool("get_collection_run_manifest", "Return bounded provenance for one collection run.",
          read_only=True, idempotent=True)
    def get_collection_run_manifest(operation_id: int | None = None) -> dict:
        return _session_call(session.get_run_manifest, operation_id)

    @tool("compare_collection_runs", "Compare parameters and aggregate metrics for two collection runs.",
          read_only=True, idempotent=True)
    def compare_collection_runs(before_operation_id: int, after_operation_id: int) -> dict:
        return _session_call(session.compare_runs, before_operation_id, after_operation_id)

    @tool("export_collection_regression_fixture",
          "Write a bounded JSON description for one retained failed collection run.")
    def export_collection_regression_fixture(operation_id: int, path: str) -> dict:
        return _session_call(session.export_regression_fixture, operation_id, path)

    @tool("save_collection_checkpoint", "Save a new HDF5 collection checkpoint in managed storage.")
    def save_collection_checkpoint(path: str) -> dict:
        return _session_call(session.save_checkpoint, path)

    @tool("load_collection_checkpoint", "Load a permitted HDF5 checkpoint as a managed collection copy.")
    def load_collection_checkpoint(path: str) -> dict:
        return _session_call(session.load_checkpoint, path)

    return server


def register_unichrom_tools(server, session, annotations_factory):
    """Register the bounded UniChrom workflow contract."""

    def tool(name, description, *, read_only=False, idempotent=False):
        annotations = annotations_factory(
            readOnlyHint=read_only,
            destructiveHint=False,
            idempotentHint=idempotent,
            openWorldHint=False,
        )
        return server.tool(name=name, description=description, annotations=annotations)

    @tool("open_chromatogram", "Open a permitted chromatogram into a managed HDF5 project.")
    def open_chromatogram(path: str) -> dict:
        return _session_call(session.open_chromatogram, path)

    @tool("get_chromatogram_state", "Return bounded TIC/BPC, scan, window, status, and freshness metadata.",
          read_only=True, idempotent=True)
    def get_chromatogram_state() -> dict:
        return _session_call(session.get_state)

    @tool("get_parameters", "Return reviewed UniChrom and deconvolution parameter metadata.",
          read_only=True, idempotent=True)
    def get_parameters(names: list[str] | None = None) -> dict:
        return _session_call(session.get_parameters, names)

    @tool("set_parameters", "Atomically validate and update reviewed UniChrom parameters.")
    def set_parameters(values: dict) -> dict:
        return _session_call(session.set_parameters, values)

    @tool("create_time_windows", "Replace the managed collection with regular time windows.")
    def create_time_windows(width: float, start: float | None = None,
                            end: float | None = None) -> dict:
        return _session_call(session.create_time_windows, width, start, end)

    @tool("run_unichrom", "Process, deconvolve, and pick aggregate peaks for the managed windows.")
    def run_unichrom() -> dict:
        return _session_call(session.run, "all")

    @tool("get_unichrom_results", "Return bounded aggregate results for the current windows.",
          read_only=True, idempotent=True)
    def get_unichrom_results(max_peaks: int = 100) -> dict:
        return _session_call(session.get_results, max_peaks)

    @tool("get_unichrom_run_manifest", "Return bounded provenance for one UniChrom run.",
          read_only=True, idempotent=True)
    def get_unichrom_run_manifest(operation_id: int | None = None) -> dict:
        return _session_call(session.get_run_manifest, operation_id)

    @tool("compare_unichrom_runs", "Compare parameters and aggregate metrics for two UniChrom runs.",
          read_only=True, idempotent=True)
    def compare_unichrom_runs(before_operation_id: int, after_operation_id: int) -> dict:
        return _session_call(session.compare_runs, before_operation_id, after_operation_id)

    @tool("export_unichrom_regression_fixture",
          "Write a bounded JSON description for one retained failed UniChrom run.")
    def export_unichrom_regression_fixture(operation_id: int, path: str) -> dict:
        return _session_call(session.export_regression_fixture, operation_id, path)

    return server


def create_server(session, gui_bridge=None, collection=False, unichrom=False):
    """Create the FastMCP server, importing the optional dependency only when used."""
    try:
        from mcp.server.fastmcp import FastMCP
        from mcp.types import ToolAnnotations
    except ImportError as error:
        raise RuntimeError(
            'The MCP adapter requires the optional agent dependencies. '
            'Install them with: python -m pip install -e ".[agent]"'
        ) from error

    if collection:
        instructions = COLLECTION_SERVER_INSTRUCTIONS
    elif unichrom:
        instructions = UNICHROM_SERVER_INSTRUCTIONS
    else:
        instructions = SERVER_INSTRUCTIONS
    server = FastMCP("UniDec", instructions=instructions)
    if collection:
        return register_collection_tools(server, session, ToolAnnotations)
    if unichrom:
        return register_unichrom_tools(server, session, ToolAnnotations)
    return register_tools(server, session, ToolAnnotations, gui_bridge=gui_bridge)


def build_parser():
    parser = argparse.ArgumentParser(description="Serve one headless UniDec session over MCP STDIO.")
    parser.add_argument(
        "--mode",
        choices=("ms", "imms", "cdms", "metaunidec", "unichrom"),
        default="ms",
        help="UniDec engine mode for this server (default: ms).",
    )
    parser.add_argument(
        "--input-root",
        action="append",
        required=True,
        dest="input_roots",
        help="Directory that the server may read (repeat for additional roots).",
    )
    parser.add_argument(
        "--working-directory",
        required=True,
        help="Directory reserved for session outputs, plots, and checkpoints.",
    )
    parser.add_argument(
        "--gui-bridge-file",
        help="Descriptor written by a live conventional UniDec GUI agent bridge.",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    with redirect_stdout(sys.stderr):
        if args.mode == "metaunidec":
            from unidec.meta_session import MetaUniDecCollectionSession

            session = MetaUniDecCollectionSession(
                working_directory=args.working_directory,
                allowed_input_roots=args.input_roots,
            )
        elif args.mode == "unichrom":
            from unidec.unichrom_session import UniChromSession

            session = UniChromSession(
                working_directory=args.working_directory,
                allowed_input_roots=args.input_roots,
            )
        else:
            engine = None
            engine_adapter = None
            if args.mode == "imms":
                from unidec.modules.IMEng import UniDecIM

                engine = UniDecIM(ignore_args=True, silent=True)
            elif args.mode == "cdms":
                from unidec.modules.CDEng import UniDecCD

                engine = UniDecCD()
                engine_adapter = UniDecCDSessionEngineAdapter(engine)
            session = UniDecSession(
                engine=engine,
                engine_adapter=engine_adapter,
                working_directory=args.working_directory,
                allowed_input_roots=args.input_roots,
            )
    gui_bridge = GuiAgentBridgeClient(args.gui_bridge_file) if args.gui_bridge_file else None
    server = create_server(
        session,
        gui_bridge=gui_bridge,
        collection=args.mode == "metaunidec",
        unichrom=args.mode == "unichrom",
    )
    server.run(transport="stdio")


if __name__ == "__main__":
    main()
