from __future__ import annotations
from typing import Optional, Tuple
import os
import numpy as np
import matplotlib.pyplot as plt

# Use the previously provided thermo_raw_utils logic inline (simplified import here)
# If you already have these helpers in a separate file, you can import them;
# otherwise include a minimal wrapper that calls plot_raw_chrom_and_spectrum from your utils.

# To avoid circular complexity, we provide a thin wrapper that expects raw_readers installed.
try:
    from thermo_raw_utils import plot_raw_chrom_and_spectrum as _plot_raw
except Exception:
    _plot_raw = None

def plot_raw_chrom_and_spectrum(
    raw_path: str,
    ms_level: int = 1,
    max_scans: Optional[int] = 200,
    spectrum_mode: str = "apex",
    spectrum_index: Optional[int] = None,
    normalize_spectrum: bool = True,
    figsize: Tuple[float, float] = (10, 6),
    dpi: int = 120,
    font_name: str = "Arial",
    font_size: int = 10,
    save: Optional[str] = None,
    show: bool = False,
    prefer_native: bool = True,
    msconvert_path: Optional[str] = None,
):
    """
    Delegate to thermo_raw_utils.plot_raw_chrom_and_spectrum if available.
    """
    if _plot_raw is None:
        raise RuntimeError("Thermo RAW utilities not available. Ensure thermo_raw_utils.py is present and imports succeed.")
    return _plot_raw(
        raw_path=raw_path,
        ms_level=ms_level,
        max_scans=max_scans,
        spectrum_mode=spectrum_mode,
        spectrum_index=spectrum_index,
        normalize_spectrum=normalize_spectrum,
        figsize=figsize,
        dpi=dpi,
        font_name=font_name,
        font_size=font_size,
        save=save,
        show=show,
        prefer_native=prefer_native,
        msconvert_path=msconvert_path,
    )