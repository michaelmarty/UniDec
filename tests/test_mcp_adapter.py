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
    "save_state",
    "load_state",
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
        read_only = {"get_state", "get_parameters", "get_results"}
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
            "save_state": ("checkpoint.zip",),
            "load_state": ("checkpoint.zip",),
        }
        for name, args in calls.items():
            function = self.server.tools[name][0]
            response = function(*args)
            self.assertTrue(response["ok"])
            self.assertEqual(response["method"], name)
        self.assertEqual([name for name, _ in self.session.calls], list(calls))

    def test_cli_requires_explicit_roots_and_supports_multiple_input_roots(self):
        arguments = build_parser().parse_args([
            "--input-root", "first", "--input-root", "second", "--working-directory", "work",
        ])
        self.assertEqual(arguments.input_roots, ["first", "second"])
        self.assertEqual(arguments.working_directory, "work")


if __name__ == "__main__":
    unittest.main()
