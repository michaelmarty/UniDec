"""Run the GUI-independent test suite while rejecting GUI imports."""

import importlib.abc
import importlib.util
from pathlib import Path
import sys
import unittest


TEST_DIRECTORY = Path(__file__).resolve().parent
PACKAGE_DIRECTORY = TEST_DIRECTORY.parent
HEADLESS_TEST_MODULES = (
    "tests.test_batch_exports",
    "tests.test_batch_parameters",
    "tests.test_bundled_cd_processing",
    "tests.test_bundled_config_compatibility",
    "tests.test_bundled_im_processing",
    "tests.test_bundled_ms_processing",
    "tests.test_cd_engine",
    "tests.test_cd_scoring",
    "tests.test_engine_failures",
    "tests.test_engine_history",
    "tests.test_engine_state",
    "tests.test_ht_engine",
    "tests.test_im_engine_split",
    "tests.test_isodec_output_names",
    "tests.test_matchtools",
    "tests.test_mcp_adapter",
    "tests.test_meta_dataset_io",
    "tests.test_parameter_metadata",
    "tests.test_peak_quantitation",
    "tests.test_plotting_labels",
    "tests.test_tools_extraction",
    "tests.test_tools_peaks",
    "tests.test_tools_preprocessing",
    "tests.test_uccd_binary",
    "tests.test_unichrom_cd_stack_shapes",
    "tests.test_unichrom_workflow",
    "tests.test_unidec_engine_io",
    "tests.test_unidec_session",
)
BLOCKED_GUI_MODULES = (
    "wx",
    "unidec.modules.gui_elements",
    "unidec.metaunidec.gui_elements",
)


def is_gui_module(module_name):
    return any(module_name == root or module_name.startswith(root + ".") for root in BLOCKED_GUI_MODULES)


class HeadlessImportGuard(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if is_gui_module(fullname):
            raise ModuleNotFoundError(
                f"Headless tests cannot import GUI dependency {fullname!r}; "
                "move GUI-independent behavior into a headless module"
            )
        return None


def install_headless_import_guard():
    original_find_spec = importlib.util.find_spec

    def headless_find_spec(name, package=None):
        if is_gui_module(name):
            return None
        return original_find_spec(name, package)

    importlib.util.find_spec = headless_find_spec
    sys.meta_path.insert(0, HeadlessImportGuard())


def main():
    if str(PACKAGE_DIRECTORY) not in sys.path:
        sys.path.insert(0, str(PACKAGE_DIRECTORY))
    install_headless_import_guard()
    suite = unittest.defaultTestLoader.loadTestsFromNames(HEADLESS_TEST_MODULES)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()
