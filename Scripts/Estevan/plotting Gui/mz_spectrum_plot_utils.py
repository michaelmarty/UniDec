# mz_spectrum_plot_utils.py
from __future__ import annotations
from typing import Sequence, Tuple, Optional, Union, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _detect_columns(df: pd.DataFrame) -> Tuple[str, str]:
    """
    Detect m/z and intensity columns (case-insensitive).
    Accepts possible m/z names: 'mz', 'm_z', 'm/z'.
    Accepts intensity names: 'intensity', 'Intensity', 'y', 'signal'.
    Fallback: first two columns.
    Returns (mz_col, intensity_col).
    """
    cols_lower = {str(c).lower(): c for c in df.columns}
    mz_candidates = ["mz", "m_z", "m/z"]
    y_candidates = ["intensity", "y", "signal"]

    mz_col = None
    y_col = None

    for key in mz_candidates:
        if key in cols_lower:
            mz_col = cols_lower[key]
            break
    for key in y_candidates:
        if key in cols_lower:
            y_col = cols_lower[key]
            break

    # Fallback: use the first two columns
    if mz_col is None or y_col is None:
        ordered = list(df.columns)
        if len(ordered) < 2:
            raise ValueError("DataFrame must have at least two columns for m/z and intensity.")
        mz_col = ordered[0]
        y_col = ordered[1]

    return mz_col, y_col


def _load_spectrum(source: Union[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Load or validate mass spectrum data with columns: m/z (or variants) and intensity.
    Adds 'mz' and 'Percent' columns (normalized intensity).
    - If source is a file path, expects whitespace-delimited text with at least two columns:
      first: m/z, second: intensity (headerless).
    - If source is a DataFrame, auto-detects m/z and intensity columns.
    """
    if isinstance(source, str):
        df = pd.read_csv(source, sep=r"\s+", header=None, comment="#", engine="python")
        if df.shape[1] < 2:
            raise ValueError("Input file must have at least two columns: m/z and intensity.")
        df = df.iloc[:, :2].copy()
        df.columns = ["mz", "Intensity"]
    elif isinstance(source, pd.DataFrame):
        df = source.copy()
        mz_col, y_col = _detect_columns(df)
        df = df[[mz_col, y_col]].copy()
        df.columns = ["mz", "Intensity"]
    else:
        raise TypeError("source must be a file path (str) or a pandas DataFrame.")

    # Ensure numeric
    df["mz"] = pd.to_numeric(df["mz"], errors="coerce")
    df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
    df = df.dropna().sort_values("mz").reset_index(drop=True)

    # Normalize to percent
    max_i = df["Intensity"].max()
    df["Percent"] = (df["Intensity"] / max_i) * 100.0 if max_i > 0 else 0.0
    return df


def plot_single_mz_spectrum(
    dataset: Union[str, pd.DataFrame],
    label: str,
    second_label: Optional[str] = None,
    x_range: Optional[Tuple[float, float]] = None,
    color: str = "#2C7FB8",
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Plot a single normalized m/z spectrum (Percent) with optional second legend near the peak.
    Returns (fig, ax).
    """
    plt.rcParams.update({"font.family": font_name, "font.size": font_size})
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

    df = _load_spectrum(dataset)
    ax.plot(df["mz"], df["Percent"], color=color, linewidth=1.8, label=label)

    # Axis ranges
    if x_range is not None:
        ax.set_xlim(*x_range)
    else:
        # Use data-driven x-range with small margins
        xmin, xmax = float(df["mz"].min()), float(df["mz"].max())
        xr = xmax - xmin if xmax > xmin else 1.0
        ax.set_xlim(xmin - 0.02 * xr, xmax + 0.02 * xr)

    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative Intensity (%)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper left")

    # Optional second legend near the main peak
    if second_label:
        peak_idx = int(df["Percent"].idxmax())
        peak_x = float(df["mz"].iloc[peak_idx])
        peak_y = float(df["Percent"].iloc[peak_idx])
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        xr = max(x1 - x0, 1e-9)
        yr = max(y1 - y0, 1e-9)
        xf = float(np.clip((peak_x - x0) / xr + 0.05, 0.02, 0.98))
        yf = float(np.clip((peak_y - y0) / yr, 0.02, 0.98))
        handle = Line2D([0], [0], color=color, lw=1.5)
        leg2 = ax.legend([handle], [second_label], loc="center", bbox_to_anchor=(xf, yf), frameon=False, fontsize=font_size)
        ax.add_artist(leg2)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax


def plot_mz_spectrum_offset_stack(
    datasets: Sequence[Union[str, pd.DataFrame]],
    labels: Sequence[str],
    second_labels: Optional[Sequence[Optional[str]]] = None,
    x_range: Optional[Tuple[float, float]] = None,
    colors: Optional[Sequence[str]] = None,
    y_offset_step: float = 120.0,
    y_tick_labels: Sequence[str] = ("0", "%", "100"),
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Plot multiple normalized m/z spectra with vertical offsets (Percent + offset),
    optional second legend near each peak. Returns (fig, axes).
    """
    n = len(datasets)
    if len(labels) != n:
        raise ValueError("labels length must match datasets length.")
    if second_labels is None:
        second_labels = [None] * n
    if len(second_labels) != n:
        raise ValueError("second_labels length must match datasets length (or None).")

    # Colors
    if colors is None or len(colors) < n:
        base = ["#2C7FB8", "#CC79A7", "#7A7A7A", "#00A9A5", "#F0E442"]
        k = (n + len(base) - 1) // len(base)
        colors = (base * k)[:n]
    else:
        colors = list(colors)[:n]

    # Global style
    plt.rcParams.update({"font.family": font_name, "font.size": font_size})
    fig, axes = plt.subplots(n, 1, sharex=True, figsize=figsize, dpi=dpi)
    if n == 1:
        axes = [axes]

    # Prepare all spectra for a global x-range if not provided
    spectra = [ _load_spectrum(src) for src in datasets ]
    if x_range is None:
        xmin = min(float(df["mz"].min()) for df in spectra)
        xmax = max(float(df["mz"].max()) for df in spectra)
        xr = xmax - xmin if xmax > xmin else 1.0
        x_range = (xmin - 0.02 * xr, xmax + 0.02 * xr)

    for i, (df, label, second, color) in enumerate(zip(spectra, labels, second_labels, colors)):
        ax = axes[i]
        offset = i * y_offset_step

        # Plot normalized spectrum with offset
        ax.plot(df["mz"], df["Percent"] + offset, color=color, linewidth=1.5, label=label)

        # Style
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(*x_range)
        ax.set_ylim(offset, offset + 101.5)
        ax.set_yticks([offset, offset + 50, offset + 100])
        ax.set_yticklabels(list(y_tick_labels))

        # First legend (upper-left)
        leg1 = ax.legend([ax.lines[-1]], [label], loc="upper left", frameon=False, fontsize=font_size)
        ax.add_artist(leg1)

        # Optional second legend near peak
        if second:
            peak_idx = int(df["Percent"].idxmax())
            peak_x = float(df["mz"].iloc[peak_idx])
            peak_y = float(df["Percent"].iloc[peak_idx] + offset)
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            xr = max(x1 - x0, 1e-9)
            yr = max(y1 - y0, 1e-9)
            xf = float(np.clip((peak_x - x0) / xr + 0.05, 0.02, 0.98))
            yf = float(np.clip((peak_y - y0) / yr, 0.02, 0.98))
            handle = Line2D([0], [0], color=color, lw=1.5)
            leg2 = ax.legend([handle], [second], loc="center", bbox_to_anchor=(xf, yf), frameon=False, fontsize=font_size)
            ax.add_artist(leg2)

        # Bottom x-label only on last axis
        if i < n - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("m/z")

    # Shared y label
    fig.text(0.05, 0.5, "Relative Intensity (%)", va="center", rotation="vertical", fontsize=font_size)

    # Layout
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95, hspace=0.6)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()

    return fig, axes