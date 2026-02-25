# mz_spectrum_plot_utils.py
from __future__ import annotations
from typing import Sequence, Tuple, Optional, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# Consistent sizing for stacked plots
DEFAULT_WIDTH_IN = 10.0     # keep width constant; height scales with rows
ROW_HEIGHT_IN     = 1.35     # per subplot
MIN_HEIGHT_IN     = 3.8      # minimum total height


def read_mz_spectrum_file(file_path: str) -> pd.DataFrame:
    """
    Read an m/z spectrum text file.
    Expected format: two columns (mz, intensity), with or without header.
    The parser is resilient to whitespace- or comma-separated data.

    Returns DataFrame with columns ['mz', 'intensity'] and drops NaNs.
    """
    # Try pandas with common separators
    for sep, header in [(None, None), (r"\s+", None), (",", None), (",", 0), (r"\s+", 0)]:
        try:
            df = pd.read_csv(file_path, sep=sep, header=header, engine="python")
            # Reduce to first two numeric columns
            num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            if len(num_cols) >= 2:
                df = df[num_cols[:2]].copy()
            else:
                df = df.iloc[:, :2].copy()
            df.columns = ["mz", "intensity"]
            df = df.dropna()
            df["mz"] = pd.to_numeric(df["mz"], errors="coerce")
            df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
            df = df.dropna()
            if len(df) > 2:
                return df
        except Exception:
            continue

    # Final fallback: numpy loadtxt
    arr = np.loadtxt(file_path)
    if arr.ndim == 1:
        arr = np.column_stack([np.arange(len(arr)), arr])
    df = pd.DataFrame(arr[:, :2], columns=["mz", "intensity"]).dropna()
    return df


def _apply_style(font_name: str, font_size: int, tick_length: float, tick_width: float):
    mpl.rcParams.update({
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
    })


def _place_second_label_near_peak(
    ax: mpl.axes.Axes,
    x: np.ndarray,
    y: np.ndarray,
    label: Optional[str],
    color: str,
    x_range: Tuple[float, float],
    font_size: int,
):
    """
    Place 'second legend' text near the apex of the peak, to the right of max.
    """
    if not label:
        return
    mask = (x >= x_range[0]) & (x <= x_range[1])
    if not np.any(mask):
        mask = np.ones_like(x, dtype=bool)
    x_sel = x[mask]
    y_sel = y[mask]
    i_max = int(np.argmax(y_sel))
    x_peak = float(x_sel[i_max])
    y_peak = float(y_sel[i_max])

    dx = 0.04 * (x_range[1] - x_range[0])
    x_text = min(x_range[1] - 0.01 * (x_range[1] - x_range[0]), x_peak + dx)
    ax.text(
        x_text, y_peak,
        label,
        color=color,
        va="center", ha="left",
        fontsize=font_size,
    )


def plot_single_mz_spectrum(
    dataset: str,
    label: str,
    second_label: Optional[str] = None,
    x_range: Tuple[float, float] = (100.0, 2000.0),
    color: str = "#2C7FB8",
    normalize: bool = True,
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (10, 4),
    dpi: int = 300,
    tick_length: float = 5.0,
    tick_width: float = 1.0,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Single m/z spectrum with optional 'second label near peak'.
    When normalize=True, intensities are scaled to percent and y-ticks show
    $$[0, \%, 100]$$.
    """
    _apply_style(font_name, font_size, tick_length, tick_width)

    df = read_mz_spectrum_file(dataset)
    df = df[(df["mz"] >= x_range[0]) & (df["mz"] <= x_range[1])].copy()
    if len(df) == 0:
        df = read_mz_spectrum_file(dataset)

    x = df["mz"].astype(float).values
    y = df["intensity"].astype(float).values

    if normalize:
        m = float(np.nanmax(y))
        if m > 0:
            y = (y / m) * 100.0

    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

    ax.plot(x, y, color=color, linewidth=2.0, label=label)

    # Limits and margins
    ax.set_xlim(*x_range)
    ax.margins(y=0)

    # Clean spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Y ticks
    if normalize:
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 50, 100])
        ax.set_yticklabels(["0", "%", "100"])
        ax.set_ylabel("Relative Intensity (%)")
    else:
        ax.set_ylabel("Intensity (a.u.)")

    ax.set_xlabel("m/z")

    # Primary legend: frameless, left-center
    ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.5))

    # Second label near peak (colored text)
    _place_second_label_near_peak(
        ax=ax, x=x, y=y, label=second_label, color=color,
        x_range=x_range, font_size=font_size
    )

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()

    return fig, ax


def plot_mz_spectrum_offset_stack(
    datasets: Sequence[str],
    labels: Sequence[str],
    second_labels: Sequence[Optional[str]],
    x_range: Tuple[float, float] = (100.0, 2000.0),
    colors: Optional[Sequence[str]] = None,
    y_offset_step: float = 120.0,  # kept for API compatibility; not used in per-row stack
    y_tick_labels: Tuple[str, str, str] = ("0", "%", "100"),
    normalize: bool = True,
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (10, 8),
    dpi: int = 300,
    tick_length: float = 5.0,
    tick_width: float = 1.0,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Stacked m/z spectra (one axis per dataset), styled consistently.

    - Each subplot has a left-center legend (sample label).
    - Optional 'second label near peak' is placed as colored text to the right of the apex.
    - When normalize=True, y-ticks show $$[0, \%, 100]$$.
    - Figure height scales with number of rows to keep aspect ratio consistent.

    Note: y_offset_step is retained for compatibility but not needed in per-row stacking.
    """
    n = len(datasets)
    if len(labels) != n:
        raise ValueError("labels length must match datasets length.")
    if second_labels is None or len(second_labels) != n:
        second_labels = [None] * n

    if colors is None or len(colors) < n:
        base = ["#2C7FB8", "#CC79A7", "#00A9A5", "#7A7A7A", "#F0E442", "#56B4E9"]
        k = (n + len(base) - 1) // len(base)
        colors = (base * k)[:n]
    else:
        colors = list(colors)[:n]

    _apply_style(font_name, font_size, tick_length, tick_width)

    # Compute consistent height
    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = max(MIN_HEIGHT_IN, (ROW_HEIGHT_IN * max(1, n)) + 0.8)

    fig, axes = plt.subplots(n, 1, figsize=(width_in, target_h), sharex=True, dpi=dpi)
    if n == 1:
        axes = [axes]

    for i, (path, lab, sec_lab, color) in enumerate(zip(datasets, labels, second_labels, colors)):
        df = read_mz_spectrum_file(path)
        df = df[(df["mz"] >= x_range[0]) & (df["mz"] <= x_range[1])].copy()
        if len(df) == 0:
            df = read_mz_spectrum_file(path)

        x = df["mz"].astype(float).values
        y = df["intensity"].astype(float).values

        if normalize:
            m = float(np.nanmax(y))
            if m > 0:
                y = (y / m) * 100.0

        ax = axes[i]
        ax.plot(x, y, color=color, linewidth=2.0, label=lab)

        # Limits and margins
        ax.set_xlim(*x_range)
        if normalize:
            ax.set_ylim(0, 100)
        ax.margins(y=0)

        # Clean spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Y ticks and labels
        if normalize:
            ax.set_yticks([0, 50, 100])
            if y_tick_labels and len(y_tick_labels) == 3:
                ax.set_yticklabels(list(y_tick_labels))
        else:
            ymin = float(np.nanmin(y))
            ymax = float(np.nanmax(y))
            ax.set_yticks([ymin, ymax])

        # Tick styling
        ax.tick_params(axis="y", length=tick_length, width=tick_width)
        ax.tick_params(axis="x", length=tick_length, width=tick_width)

        # Primary legend on left
        ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.5))

        # Second label near peak
        _place_second_label_near_peak(
            ax=ax, x=x, y=y, label=sec_lab, color=color,
            x_range=x_range, font_size=font_size
        )

        # Hide x tick labels on intermediate rows
        if i < n - 1:
            ax.tick_params(labelbottom=False)

    # Shared labels
    fig.supxlabel("m/z")
    fig.supylabel("Relative Intensity (%)" if normalize else "Intensity (a.u.)")

    # Layout
    plt.tight_layout(rect=(0.08, 0.08, 0.98, 0.96), h_pad=0.6)
    plt.subplots_adjust(hspace=0.55)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()

    return fig, axes