"""Generate Sphinx API pages while excluding import-unsafe scripts."""

import argparse
import csv
from io import StringIO
from pathlib import Path
import sys

from sphinx.ext.apidoc import main as apidoc_main


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIRECTORY = REPOSITORY_ROOT / "docs" / "api" / "generated"
PACKAGE_DIRECTORY = REPOSITORY_ROOT / "unidec"

# These are executable training, test, data-generation, or native-library
# wrappers rather than import-safe public API modules.
EXCLUDE_PATTERNS = [
    PACKAGE_DIRECTORY / "__main__.py",
    PACKAGE_DIRECTORY / "gui_self_test.py",
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


def write_parameter_reference(output_directory: Path) -> None:
    """Render the public registry into a Sphinx table."""
    from unidec.modules.unidecstructure import _PARAMETER_SPECS

    rows = StringIO()
    writer = csv.writer(rows, lineterminator="\n", quoting=csv.QUOTE_ALL)
    for spec in _PARAMETER_SPECS:
        if not spec.public:
            continue
        bounds = []
        if spec.minimum is not None:
            bounds.append(f"≥ {spec.minimum:g}")
        if spec.maximum is not None:
            bounds.append(f"≤ {spec.maximum:g}")
        if spec.choices is not None:
            bounds.append("One of " + ", ".join(map(str, spec.choices)))
        writer.writerow([
            f"``{spec.attribute}``", spec.dtype.__name__, spec.description,
            spec.units or "—", "; ".join(bounds) or "—", ", ".join(spec.modes),
            "Yes" if spec.optional else "No", spec.affects,
        ])

    title = "Public parameter reference"
    page = [
        ".. Generated from ParameterSpec by generate_api_docs.py.",
        "",
        title,
        "=" * len(title),
        "",
        "Bounds are inclusive. Optional parameters accept ``None`` to clear a value.",
        "The stage column is registry metadata. The conventional Python session",
        "currently clears freshness for all stages after a successful update.",
        "",
        ".. csv-table::",
        '   :header: "Name", "Type", "Description", "Units", "Bounds or choices", "Modes", "Optional", "Stage"',
        "",
        *(f"   {row}" for row in rows.getvalue().splitlines()),
        "",
    ]
    (output_directory / "parameters.rst").write_text("\n".join(page), encoding="utf-8")


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
    sys.path.insert(0, str(REPOSITORY_ROOT))
    write_parameter_reference(output_directory)
