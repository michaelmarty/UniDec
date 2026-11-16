"""Generate Sphinx API pages while excluding import-unsafe scripts."""

import argparse
from pathlib import Path

from sphinx.ext.apidoc import main as apidoc_main


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIRECTORY = REPOSITORY_ROOT / "docs" / "api" / "generated"
PACKAGE_DIRECTORY = REPOSITORY_ROOT / "unidec"

# These are executable training, test, data-generation, or native-library
# wrappers rather than import-safe public API modules.
EXCLUDE_PATTERNS = [
    PACKAGE_DIRECTORY / "IsoDec" / "*.so",
    PACKAGE_DIRECTORY / "IsoDec" / "train.py",
    PACKAGE_DIRECTORY / "IsoDec" / "isogenwrapper.py",
    PACKAGE_DIRECTORY / "IsoDec" / "IsoGen" / "isogenc.py",
    PACKAGE_DIRECTORY / "IsoDec" / "IsoGen" / "*training*.py",
    PACKAGE_DIRECTORY / "metaunidec" / "test_MUD.py",
    PACKAGE_DIRECTORY / "modules" / "unidecwrapper.py",
]


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=SOURCE_DIRECTORY)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    output_directory = arguments.output_dir.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    apidoc_main(
        [
            "--force",
            "--remove-old",
            "--module-first",
            "--output-dir",
            str(output_directory),
            str(PACKAGE_DIRECTORY),
            *(str(pattern) for pattern in EXCLUDE_PATTERNS),
        ]
    )
