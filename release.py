"""Local pre-release validation and artifact build."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


def run(command):
    print("+", *command)
    subprocess.run(command, check=True)


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-tests", action="store_true")
    args = parser.parse_args(argv)

    if subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True, check=True).stdout.strip():
        raise SystemExit("Commit or stash working-tree changes before releasing")
    run(["git", "submodule", "update", "--init", "--recursive"])
    if not args.skip_tests:
        run([sys.executable, "-m", "pytest"])
    for path in (Path("build"), Path("dist"), Path("wheelhouse")):
        shutil.rmtree(path, ignore_errors=True)
    run([sys.executable, "-m", "build"])
    artifacts = [str(path) for path in Path("dist").iterdir()]
    run([sys.executable, "-m", "twine", "check", *artifacts])


if __name__ == "__main__":
    main()
