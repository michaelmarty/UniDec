import os
import shutil
import sys
from importlib.resources import as_file, files
from pathlib import Path


def has_gui_display():
    if sys.platform.startswith("win") or sys.platform == "darwin":
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


def copy_unidec_example(destination, *relative_parts):
    resource = files("unidec").joinpath("bin", "Example Data", *relative_parts)
    with as_file(resource) as source:
        if not source.is_file():
            raise FileNotFoundError(f"Packaged UniDec example not found: {source}")
        target = Path(destination, source.name)
        shutil.copy2(source, target)
    return target


def find_importer_test_data():
    """Locate the separately distributed UniDecImporter fixtures."""
    configured = os.environ.get("UNIDEC_IMPORTER_TEST_DATA")
    candidates = []
    if configured:
        candidates.append(Path(configured))

    # Common source layout: UniDecDev and UniDecImporter are sibling checkouts.
    candidates.append(Path(__file__).resolve().parents[4] / "UniDecImporter" / "TestData")

    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    return None


def importer_fixture(test_case, *relative_parts, required=True):
    """Return a real UniDecImporter fixture, rejecting missing Git LFS objects."""
    root = find_importer_test_data()
    path = None if root is None else root.joinpath(*relative_parts)
    available = path is not None and path.exists()
    if available and path.is_file():
        with path.open("rb") as fixture:
            available = not fixture.read(64).startswith(b"version https://git-lfs.github.com/spec/")
    if available:
        return path

    relative = os.path.join(*relative_parts)
    message = f"UniDecImporter fixture is missing or is only a Git LFS pointer: {relative}"
    if required and os.environ.get("UNIDEC_REQUIRE_IMPORTER_TEST_DATA") == "1":
        test_case.fail(message)
    test_case.skipTest(message)


def copy_workflow_spectra(destination):
    """Prefer importer fixtures, falling back to examples in the UniDec wheel."""
    importer_data = find_importer_test_data()
    if importer_data is not None:
        single_scan = importer_data / "SingleScan"
        sources = [single_scan / "test_txt.txt", single_scan / "test_csv.csv"]
        if all(source.is_file() for source in sources):
            targets = []
            for source in sources:
                target = Path(destination, source.name)
                shutil.copy2(source, target)
                targets.append(target)
            return targets

    return [copy_unidec_example(destination, name) for name in ("ADH.txt", "BSA.txt")]
