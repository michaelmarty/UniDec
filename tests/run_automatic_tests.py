"""Run package tests while excluding tests that target research scripts."""

import ast
from pathlib import Path
import sys
import unittest


TEST_DIRECTORY = Path(__file__).resolve().parent
PACKAGE_DIRECTORY = TEST_DIRECTORY.parent
SCRIPT_ROOTS = {"Scripts", "PublicScripts"}

if str(PACKAGE_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIRECTORY))


def targets_script_directories(test_path):
    """Return whether a test imports or constructs a path into a script tree."""
    tree = ast.parse(test_path.read_text(encoding="utf-8"), filename=str(test_path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            components = node.value.replace("\\", "/").split("/")
            if SCRIPT_ROOTS.intersection(components):
                return True
        elif isinstance(node, ast.Import):
            if any(alias.name.split(".", 1)[0] in SCRIPT_ROOTS for alias in node.names):
                return True
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.split(".", 1)[0] in SCRIPT_ROOTS:
                return True
    return False


def build_suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    excluded = []
    for test_path in sorted(TEST_DIRECTORY.glob("test_*.py")):
        if targets_script_directories(test_path):
            excluded.append(test_path.name)
            continue
        suite.addTests(loader.discover(
            start_dir=str(TEST_DIRECTORY),
            pattern=test_path.name,
            top_level_dir=str(TEST_DIRECTORY),
        ))
    return suite, excluded


def main():
    suite, excluded = build_suite()
    print(f"Excluded {len(excluded)} test modules that target Scripts/** or PublicScripts/**")
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()
