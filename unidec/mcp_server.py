"""Optional MCP adapter for the headless :class:`UniDecSession` API."""

import argparse
from contextlib import redirect_stdout
import sys

from unidec.session import UniDecCDSessionEngineAdapter, UniDecSession
from unidec.gui_agent_bridge import GuiAgentBridgeClient


SERVER_INSTRUCTIONS = """Operate only through the declared UniDec tools. Open inputs only below an administrator-configured input root; all generated artifacts and checkpoints stay below the managed working directory. Never invent paths or assume stale results are current. Use this sequence: open_spectrum, inspect or set parameters, run, then get_results or get_plot. Use get_run_manifest and compare_runs for bounded provenance and objective before/after summaries. export_regression_fixture may write a non-overwriting JSON description for a retained failed run. Mutations are serialized and may return busy. Check every response's ok/error and state fields before continuing. save_state never overwrites an existing file. If the GUI tools are available, they control a separate live conventional-MS GUI session and do not share state with the headless tools. Inspect that session with get_gui_parameters, update it with set_gui_parameters, and run it with run_gui. The headless server exposes the configured conventional MS, IM-MS, or CD-MS mode only; it does not expose a shell, arbitrary Python, generic file access, native command selection, MetaUniDec, UniChrom, DoubleDec, or SMAS."""


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


def create_server(session, gui_bridge=None):
    """Create the FastMCP server, importing the optional dependency only when used."""
    try:
        from mcp.server.fastmcp import FastMCP
        from mcp.types import ToolAnnotations
    except ImportError as error:
        raise RuntimeError(
            'The MCP adapter requires the optional agent dependencies. '
            'Install them with: python -m pip install -e ".[agent]"'
        ) from error

    server = FastMCP("UniDec", instructions=SERVER_INSTRUCTIONS)
    return register_tools(server, session, ToolAnnotations, gui_bridge=gui_bridge)


def build_parser():
    parser = argparse.ArgumentParser(description="Serve one headless UniDec session over MCP STDIO.")
    parser.add_argument(
        "--mode",
        choices=("ms", "imms", "cdms"),
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
    server = create_server(session, gui_bridge=gui_bridge)
    server.run(transport="stdio")


if __name__ == "__main__":
    main()
