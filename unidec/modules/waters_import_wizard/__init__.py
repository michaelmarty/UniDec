__author__ = 'michael.marty'

import os
import platform
from pathlib import Path


def configure_masslynx_raw_dll():
    """Point the Waters importer at UniDec's bundled MassLynx SDK DLL."""
    if platform.system() != "Windows" or "MASSLYNX_RAW_DLL" in os.environ:
        return

    dll_path = Path(__file__).resolve().parents[2] / "bin" / "MassLynxRaw.dll"
    if dll_path.is_file():
        os.environ["MASSLYNX_RAW_DLL"] = str(dll_path)


configure_masslynx_raw_dll()
