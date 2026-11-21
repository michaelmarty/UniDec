"""Focused checks for the published Python API and synthetic examples."""

import importlib
import inspect
from pathlib import Path
import pkgutil
import re
import socket
import subprocess
import unittest
from unittest.mock import patch


DOCS = Path(__file__).resolve().parents[1] / "docs"

# Review this inventory with the supported workflow contracts. Keep inherited
# operations listed under the class that supplies their documented behavior.
TIER1 = {
    "engine.md": {
        "unidec.engine.UniDec": (
            "open_file",
            "pass_data_in",
            "process_data",
            "run_unidec",
            "pick_peaks",
            "save_state",
            "load_state",
        ),
        "unidec.modules.IMEng.UniDecIM": (
            "process_data",
            "load_state",
        ),
        "unidec.modules.CDEng.UniDecCD": (
            "open_cdms_file",
            "process_data",
            "run_deconvolution",
        ),
        "unidec.metaunidec.mudeng.MetaUniDec": (
            "open",
            "process_data",
            "run_unidec",
            "pick_peaks",
        ),
        "unidec.modules.ChromEng.ChromEngine": (
            "open_chrom",
            "run_unidec",
            "pick_peaks",
        ),
    },
    "configuration.md": {
        "unidec.modules.unidecstructure.ParameterSpec": (),
        "unidec.modules.unidecstructure.UniDecConfig": (
            "parameter_schema",
            "parameter_values",
            "validate_parameters",
            "update_parameters",
        ),
    },
    "data-and-peaks.md": {
        "unidec.modules.unidecstructure.DataContainer": (),
        "unidec.modules.peakstructure.Peak": (),
        "unidec.modules.peakstructure.Peaks": (
            "add_peaks",
        ),
    },
    "sessions.md": {
        "unidec.session.UniDecSession": (
            "open_spectrum",
            "get_state",
            "get_parameters",
            "set_parameters",
            "run",
            "get_results",
            "get_plot",
            "save_state",
            "load_state",
            "get_run_manifest",
            "compare_runs",
            "export_regression_fixture",
        ),
        "unidec.session.UniDecCDSessionEngineAdapter": (),
        "unidec.meta_session.MetaUniDecCollectionSession": (
            "open_collection",
            "get_state",
            "get_parameters",
            "set_parameters",
            "run",
            "get_results",
            "save_checkpoint",
            "load_checkpoint",
            "get_run_manifest",
            "compare_runs",
            "export_regression_fixture",
        ),
        "unidec.unichrom_session.UniChromSession": (
            "open_chromatogram",
            "get_state",
            "get_parameters",
            "set_parameters",
            "create_time_windows",
            "run",
            "get_results",
            "get_run_manifest",
            "compare_runs",
            "export_regression_fixture",
        ),
        "unidec.mcp_server.main": (),
        "unidec.mcp_server.create_server": (),
    },
}

HEADLESS_EXAMPLES = {
    "api/configuration.md": ("validate-parameters",),
    "api/data-and-peaks.md": ("inspect-peak-container",),
    "guides/scientific-analysis.md": (
        "preprocess-spectrum", "score-chromatographic-peak", "fit-peak-models",
        "match-candidate-masses", "map-mass-defects",
    ),
}

CURATED_REF = re.compile(r"\{any\}`[^`\n]+ <(unidec\.[\w.]+)>`")
NAMED_PYTHON_BLOCK = re.compile(
    r"^```\{code-block\} python\n:name: ([\w-]+)\n\n(.*?)^```$", re.M | re.S,
)
PLACEHOLDER = re.compile(r"^(?:todo|tbd|fixme|placeholder|stub|not implemented)\b", re.I)


class TestSupportedDocumentation(unittest.TestCase):
    def test_inventory_matches_curated_references(self):
        for page, groups in TIER1.items():
            with self.subTest(page=page):
                expected = set(groups)
                for name, members in groups.items():
                    expected.update(f"{name}.{member}" for member in members)
                actual = set(CURATED_REF.findall((DOCS / "api" / page).read_text(encoding="utf-8")))
                self.assertEqual(expected, actual)

    def test_supported_symbols_have_contract_docstrings(self):
        for groups in TIER1.values():
            for name, members in groups.items():
                with self.subTest(symbol=name):
                    self.check_docstring(pkgutil.resolve_name(name))
                for member in members:
                    with self.subTest(symbol=f"{name}.{member}"):
                        self.check_docstring(pkgutil.resolve_name(f"{name}.{member}"))

    def check_docstring(self, obj):
        docstring = inspect.cleandoc(obj.__doc__ or "")
        summary = docstring.split("\n", 1)[0]
        self.assertRegex(summary, r"\w+\s+\w+\s+\w+", "Describe the public behavior in the summary")
        self.assertNotRegex(summary, PLACEHOLDER, "Replace the placeholder with a contract")

    def test_named_headless_examples(self):
        from UniDecImporter.ImporterFactory import ImporterFactory

        # Matplotlib can invoke font discovery during import on Linux.
        importlib.import_module("matplotlib.font_manager")
        for page, names in HEADLESS_EXAMPLES.items():
            path = DOCS / page
            source = path.read_text(encoding="utf-8")
            for name in names:
                with self.subTest(page=page, example=name):
                    blocks = [m for m in NAMED_PYTHON_BLOCK.finditer(source) if m[1] == name]
                    self.assertEqual(len(blocks), 1, "Expected one named Python code block")
                    block = blocks[0]
                    line_offset = source.count("\n", 0, block.start(2))
                    code = compile("\n" * line_offset + block[2], str(path), "exec", optimize=0)
                    with (
                        patch.object(ImporterFactory, "create_importer",
                                     side_effect=AssertionError("Use synthetic input in checked examples")),
                        patch.object(subprocess, "Popen",
                                     side_effect=AssertionError("Checked examples must run in-process")),
                        patch.object(socket.socket, "connect",
                                     side_effect=AssertionError("Checked examples must run offline")),
                        patch.object(socket.socket, "connect_ex",
                                     side_effect=AssertionError("Checked examples must run offline")),
                    ):
                        exec(code, {"__name__": "__main__"})


if __name__ == "__main__":
    from tests.run_headless_tests import install_headless_import_guard

    install_headless_import_guard()
    unittest.main()
