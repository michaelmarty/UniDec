# plot_helpers.py
from __future__ import annotations
from contextlib import contextmanager
from typing import Sequence, Optional, Tuple, Dict, Iterable
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# Constants
DEFAULT_WIDTH_IN = 10.0
ROW_HEIGHT_IN = 1.35
MIN_HEIGHT_IN = 3.8
DISPLAY_DPI = 110

DEFAULT_CHROM_COLORS = ["#333333", "#00A9A5", "#C76BA0", "#2C7FB8", "#7A7A7A", "#00A9A5"]
DEFAULT_CDMS_COLORS  = ["#CC79A7", "#2C7FB8", "#7A7A7A", "#00A9A5"]
DEFAULT_MWD_COLORS   = ["#0F0C08", "#3CAEA3"]
DEFAULT_MZ_COLORS    = ["#2C7FB8", "#CC79A7", "#00A9A5", "#7A7A7A", "#F0E442"]

def rc_for_style(font_name: str, font_size: int, tick_length: float=5.0, tick_width: float=1.0) -> Dict:
    return {
        "font.family": font_name,
        "font.size": font_size,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size,
        "ytick.labelsize": font_size,
        "legend.fontsize": font_size,
        "axes.linewidth": 1.2,
        "xtick.major.size": tick_length,
        "ytick.major.size": tick_length,
        "xtick.major.width": tick_width,
        "ytick.major.width": tick_width,
        "legend.frameon": False,
        "legend.borderpad": 0.2,
        "legend.handlelength": 2.0,
        "legend.handletextpad": 0.6,
    }

@contextmanager
def style_context(font_name: str="Arial", font_size: int=10, tick_length: float=5.0, tick_width: float=1.0):
    rc = rc_for_style(font_name, font_size, tick_length, tick_width)
    with plt.rc_context(rc):
        yield

def ensure_colors(colors: Optional[Iterable[str]], n: int, base: Sequence[str]) -> list:
    if colors is None or len(list(colors)) < n:
        k = (n + len(base) - 1) // len(base)
        return (list(base) * k)[:n]
    return list(colors)[:n]

def compute_target_height(n_rows: int, row_height: float=ROW_HEIGHT_IN, min_height: float=MIN_HEIGHT_IN, padding: float=0.8) -> float:
    return max(min_height, row_height * max(1, n_rows) + padding)

def place_second_label(ax, x: np.ndarray, y: np.ndarray, label: Optional[str], color: str, x_range: Tuple[float,float], font_size: int):
    if not label:
        return
    if x is None or y is None or x.size == 0 or y.size == 0:
        return
    mask = (x >= x_range[0]) & (x <= x_range[1])
    if not np.any(mask):
        mask = np.ones_like(x, dtype=bool)
    x_sel = x[mask]
    y_sel = y[mask]
    if x_sel.size == 0 or y_sel.size == 0:
        return
    i_max = int(np.argmax(y_sel))
    peak_x = float(x_sel[i_max]); peak_y = float(y_sel[i_max])
    x0, x1 = x_range
    xr = max(x1 - x0, 1e-9)
    xf = float(np.clip((peak_x - x0) / xr + 0.05, 0.02, 0.98))
    handle = Line2D([0], [0], color=color, lw=1.5)
    leg = ax.legend([handle], [label], loc="center", bbox_to_anchor=(xf, (peak_y - ax.get_ylim()[0]) / max(1e-9, ax.get_ylim()[1]-ax.get_ylim()[0])), frameon=False, fontsize=font_size)
    ax.add_artist(leg)