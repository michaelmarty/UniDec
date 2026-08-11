"""Print the processed CD-MS ion count for every NPZ file in a folder."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# Allow this deeply nested script to be run directly from PowerShell or PyCharm.
SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = next(
    parent
    for parent in SCRIPT_PATH.parents
    if (parent / "unidec" / "modules" / "CDEng.py").is_file()
)
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from unidec.modules.CDEng import UniDecCD


# Change this path to choose the folder used when no command-line path is given.
DEFAULT_FOLDER = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput\Amgen Rabbit IgG 21 Injections Set 3\merged"
)
REPORT_FILENAME = "CDMS Ion Counts.csv"


def natural_sort_key(path: Path) -> list[object]:
    """Sort numbered filenames as 1, 2, ..., 10 instead of 1, 10, ..., 2."""
    return [
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", path.name)
    ]


def find_npz_files(folder: Path) -> list[Path]:
    """Find NPZ files directly in folder, excluding cached subfolder files."""
    return sorted(
        (
            path
            for path in folder.iterdir()
            if path.is_file() and path.suffix.casefold() == ".npz"
        ),
        key=natural_sort_key,
    )


def process_file(path: Path) -> int:
    """Load and process one CD-MS NPZ file and return its filtered ion count."""
    engine = UniDecCD()
    engine.open_file(str(path))
    if engine.process_data(transform=False) != 1:
        raise RuntimeError("CDEng returned an empty histogram")
    return len(engine.farray)


def format_report(
    folder: Path,
    file_count: int,
    results: list[tuple[str, int]],
    failures: list[tuple[str, str]],
) -> str:
    """Create the printable ion-count table saved after processing."""
    name_width = max([len("NPZ file"), *(len(filename) for filename, _ in results)])
    separator = "-" * (name_width + 15)
    lines = [
        "CD-MS Processed Ion Counts",
        f"Folder: {folder}",
        "",
        f"{'NPZ file':<{name_width}}  {'Ion count':>12}",
        separator,
    ]
    lines.extend(
        f"{filename:<{name_width}}  {ion_count:>12,}"
        for filename, ion_count in results
    )
    lines.extend(
        [
            separator,
            f"Successful files: {len(results)}/{file_count}",
            f"Sum of per-file ion counts: {sum(count for _, count in results):,}",
        ]
    )

    if failures:
        lines.extend(["", "Failed files"])
        lines.extend(f"{filename}: {message}" for filename, message in failures)

    return "\n".join(lines) + "\n"


def save_csv_report(report_path: Path, results: list[tuple[str, int]]) -> None:
    """Save successful file counts as an Excel-compatible CSV table."""
    with report_path.open("w", newline="", encoding="utf-8-sig") as report_file:
        writer = csv.writer(report_file)
        writer.writerow(["NPZ File", "Ion Count"])
        writer.writerows(results)


def process_folder(folder: Path) -> int:
    """Process top-level NPZ files and print and save their ion counts."""
    if not folder.is_dir():
        print(f"ERROR: Folder does not exist or is not accessible: {folder}")
        return 1

    npz_files = find_npz_files(folder)
    if not npz_files:
        print(f"ERROR: No NPZ files found directly in: {folder}")
        return 1

    print(f"Found {len(npz_files)} NPZ file(s) in: {folder}")

    results: list[tuple[str, int]] = []
    failures: list[tuple[str, str]] = []
    for index, path in enumerate(npz_files, start=1):
        print(f"\n[{index}/{len(npz_files)}] Processing: {path.name}")
        try:
            ion_count = process_file(path)
        except Exception as error:
            failures.append((path.name, str(error)))
            print(f"FAILED: {path.name}: {error}")
        else:
            results.append((path.name, ion_count))
            print(f"ION COUNT: {path.name}: {ion_count:,}")

    report = format_report(folder, len(npz_files), results, failures)
    print(f"\n{report}", end="")

    report_path = folder / REPORT_FILENAME
    try:
        save_csv_report(report_path, results)
    except OSError as error:
        print(f"ERROR: Could not save report to {report_path}: {error}")
        return 1

    print(f"Saved ion-count table: {report_path}")
    return 1 if failures else 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process CD-MS NPZ files with CDEng and print their ion counts."
    )
    parser.add_argument(
        "folder",
        nargs="?",
        type=Path,
        default=DEFAULT_FOLDER,
        help="Folder containing NPZ files (defaults to DEFAULT_FOLDER in this script).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(process_folder(parse_args().folder))
