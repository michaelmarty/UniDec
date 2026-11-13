import unittest

from unidec.mcp_server import (
    build_parser,
    register_collection_tools,
    register_tools,
    register_unichrom_tools,
)


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
COLLECTION_TOOLS = {
    "open_collection",
    "get_collection_state",
    "get_parameters",
    "set_parameters",
    "run_collection",
    "get_collection_results",
    "get_collection_run_manifest",
    "compare_collection_runs",
    "export_collection_regression_fixture",
    "save_collection_checkpoint",
    "load_collection_checkpoint",
}
UNICHROM_TOOLS = {
    "open_chromatogram",
    "get_chromatogram_state",
    "get_parameters",
    "set_parameters",
    "create_time_windows",
    "run_unichrom",
    "get_unichrom_results",
    "get_unichrom_run_manifest",
    "compare_unichrom_runs",
    "export_unichrom_regression_fixture",
}


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

    def test_collection_tools_use_separate_names_and_delegate(self):
        server = FakeServer()
        register_collection_tools(server, self.session, lambda **values: values)
        self.assertEqual(set(server.tools), COLLECTION_TOOLS)
        calls = {
            "open_collection": ("project.hdf5",),
            "get_collection_state": (),
            "get_parameters": (["startz"],),
            "set_parameters": ({"startz": 5},),
            "run_collection": (),
            "get_collection_results": (4,),
            "get_collection_run_manifest": (3,),
            "compare_collection_runs": (3, 6),
            "export_collection_regression_fixture": (3, "failure.json"),
            "save_collection_checkpoint": ("checkpoint.hdf5",),
            "load_collection_checkpoint": ("checkpoint.hdf5",),
        }
        expected_methods = {
            "open_collection": "open_collection",
            "get_collection_state": "get_state",
            "get_parameters": "get_parameters",
            "set_parameters": "set_parameters",
            "run_collection": "run",
            "get_collection_results": "get_results",
            "get_collection_run_manifest": "get_run_manifest",
            "compare_collection_runs": "compare_runs",
            "export_collection_regression_fixture": "export_regression_fixture",
            "save_collection_checkpoint": "save_checkpoint",
            "load_collection_checkpoint": "load_checkpoint",
        }
        for name, args in calls.items():
            response = server.tools[name][0](*args)
            self.assertTrue(response["ok"])
            self.assertEqual(response["method"], expected_methods[name])
        read_only = {
            "get_collection_state", "get_parameters", "get_collection_results",
            "get_collection_run_manifest", "compare_collection_runs",
        }
        for name, (_, metadata) in server.tools.items():
            self.assertEqual(metadata["annotations"]["readOnlyHint"], name in read_only)

    def test_unichrom_tools_use_separate_names_and_delegate(self):
        server = FakeServer()
        register_unichrom_tools(server, self.session, lambda **values: values)
        self.assertEqual(set(server.tools), UNICHROM_TOOLS)
        calls = {
            "open_chromatogram": ("sample.mzML",),
            "get_chromatogram_state": (),
            "get_parameters": (["time_window"],),
            "set_parameters": ({"time_window": 1.0},),
            "create_time_windows": (1.0, 0.0, 5.0),
            "run_unichrom": (),
            "get_unichrom_results": (4,),
            "get_unichrom_run_manifest": (3,),
            "compare_unichrom_runs": (3, 6),
            "export_unichrom_regression_fixture": (3, "failure.json"),
        }
        expected_methods = {
            "open_chromatogram": "open_chromatogram",
            "get_chromatogram_state": "get_state",
            "get_parameters": "get_parameters",
            "set_parameters": "set_parameters",
            "create_time_windows": "create_time_windows",
            "run_unichrom": "run",
            "get_unichrom_results": "get_results",
            "get_unichrom_run_manifest": "get_run_manifest",
            "compare_unichrom_runs": "compare_runs",
            "export_unichrom_regression_fixture": "export_regression_fixture",
        }
        for name, args in calls.items():
            response = server.tools[name][0](*args)
            self.assertTrue(response["ok"])
            self.assertEqual(response["method"], expected_methods[name])

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
        self.assertEqual(
            build_parser().parse_args([
                "--mode", "metaunidec", "--input-root", "first", "--working-directory", "work",
            ]).mode,
            "metaunidec",
        )
        self.assertEqual(
            build_parser().parse_args([
                "--mode", "unichrom", "--input-root", "first", "--working-directory", "work",
            ]).mode,
            "unichrom",
        )


if __name__ == "__main__":
    unittest.main()
