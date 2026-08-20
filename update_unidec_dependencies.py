"""Upgrade UniDec's declared pip dependencies in the active Python environment."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


PROJECT_DIRECTORY = Path(__file__).resolve().parent / "public" / "UniDec"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Upgrade UniDec dependencies using public/UniDec/pyproject.toml."
    )
    parser.add_argument(
        "--extras",
        metavar="EXTRAS",
        help="Comma-separated optional dependency groups to include (for example: dev,img).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show the upgrades pip would perform without installing them.",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if not PROJECT_DIRECTORY.is_dir():
        raise FileNotFoundError(f"UniDec project directory not found: {PROJECT_DIRECTORY}")

    package_source = str(PROJECT_DIRECTORY)
    if arguments.extras:
        package_source = f"{package_source}[{arguments.extras}]"

    command = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "--upgrade",
        "--upgrade-strategy",
        "eager",
    ]
    if arguments.dry_run:
        command.append("--dry-run")
    command.append(package_source)

    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
