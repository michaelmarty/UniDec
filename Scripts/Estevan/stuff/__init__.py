"""Compatibility namespace; new scripts live in the named category folders."""
from pathlib import Path
_root = Path(__file__).resolve().parents[1]
__path__.extend(str(_root / name) for name in ("confocal", "sds_page", "sec", "mass_spec", "fluorescence_uvvis", "lipids", "dls", "bca", "radioactivity", "shared"))
