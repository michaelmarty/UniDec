import unittest

from unidec.mcp_server import build_parser, register_tools


EXPECTED_TOOLS = {
    "open_spectrum",
    "get_state",
    "get_parameters",
    "set_parameters",
    "run",
    "get_results",
    "get_plot",
    "get_run_manifest",
    "compare_runs",
    "export_regression_fixture",
    "save_state",
    "load_state",
}
GUI_TOOLS = {"get_gui_parameters", "set_gui_parameters", "run_gui"}


class FakeServer:
    def __init__(self):
        self.tools = {}

    def tool(self, **metadata):
        def decorate(function):
            self.tools[metadata["name"]] = (function, metadata)
            return function
        return decorate


class FakeSession:
    def __init__(self):
        self.calls = []

    def __getattr__(self, name):
        def call(*args):
            self.calls.append((name, args))
            return {"ok": True, "method": name, "args": list(args)}
        return call


class TestMCPAdapter(unittest.TestCase):
    def setUp(self):
        self.server = FakeServer()
        self.session = FakeSession()
        register_tools(self.server, self.session, lambda **values: values)

    def test_registers_only_the_bounded_tool_contract_with_annotations(self):
        self.assertEqual(set(self.server.tools), EXPECTED_TOOLS)
        read_only = {"get_state", "get_parameters", "get_results", "get_run_manifest", "compare_runs"}
        idempotent = read_only | {"get_plot"}
        for name, (_, metadata) in self.server.tools.items():
            annotations = metadata["annotations"]
            self.assertEqual(annotations["readOnlyHint"], name in read_only)
            self.assertEqual(annotations["idempotentHint"], name in idempotent)
            self.assertFalse(annotations["destructiveHint"])
            self.assertFalse(annotations["openWorldHint"])

    def test_every_tool_delegates_to_the_session(self):
        calls = {
            "open_spectrum": ("input.txt",),
            "get_state": (),
            "get_parameters": (["startz"],),
            "set_parameters": ({"startz": 5},),
            "run": ("process",),
            "get_results": (4,),
            "get_plot": ("mass_distribution",),
            "get_run_manifest": (3,),
            "compare_runs": (3, 6),
            "export_regression_fixture": (3, "failure.json"),
            "save_state": ("checkpoint.zip",),
            "load_state": ("checkpoint.zip",),
        }
        for name, args in calls.items():
            function = self.server.tools[name][0]
            response = function(*args)
            self.assertTrue(response["ok"])
            self.assertEqual(response["method"], name)
        self.assertEqual([name for name, _ in self.session.calls], list(calls))

    def test_optional_gui_tools_delegate_to_the_separate_bridge(self):
        bridge = FakeSession()
        server = FakeServer()
        register_tools(server, self.session, lambda **values: values, gui_bridge=bridge)

        self.assertEqual(set(server.tools), EXPECTED_TOOLS | GUI_TOOLS)
        calls = {
            "get_gui_parameters": (["startz"],),
            "set_gui_parameters": ({"startz": 5},),
            "run_gui": ("all",),
        }
        for name, args in calls.items():
            response = server.tools[name][0](*args)
            self.assertTrue(response["ok"])
        self.assertEqual(
            bridge.calls,
            [
                ("get_parameters", (["startz"],)),
                ("set_parameters", ({"startz": 5},)),
                ("run", ("all",)),
            ],
        )
        get_annotations = server.tools["get_gui_parameters"][1]["annotations"]
        self.assertTrue(get_annotations["readOnlyHint"])
        self.assertTrue(get_annotations["idempotentHint"])

    def test_cli_requires_explicit_roots_and_supports_multiple_input_roots(self):
        arguments = build_parser().parse_args([
            "--input-root", "first", "--input-root", "second", "--working-directory", "work",
        ])
        self.assertEqual(arguments.input_roots, ["first", "second"])
        self.assertEqual(arguments.working_directory, "work")
        self.assertEqual(arguments.mode, "ms")

        arguments = build_parser().parse_args([
            "--input-root", "first", "--working-directory", "work",
            "--gui-bridge-file", "bridge.json",
        ])
        self.assertEqual(arguments.gui_bridge_file, "bridge.json")
        self.assertEqual(
            build_parser().parse_args([
                "--mode", "imms", "--input-root", "first", "--working-directory", "work",
            ]).mode,
            "imms",
        )
        self.assertEqual(
            build_parser().parse_args([
                "--mode", "cdms", "--input-root", "first", "--working-directory", "work",
            ]).mode,
            "cdms",
        )


if __name__ == "__main__":
    unittest.main()
