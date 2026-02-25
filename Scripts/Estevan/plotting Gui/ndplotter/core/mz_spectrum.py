from __future__ import annotations
from typing import Sequence, Tuple, Optional, Union, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def _detect_columns(df: pd.DataFrame) -> Tuple[str, str]:
    cols_lower = {str(c).lower(): c for c in df.columns}
    mz_candidates = ["mz", "m_z", "m/z"]
    y_candidates = ["intensity", "y", "signal"]
    mz_col = next((cols_lower[k] for k in mz_candidates if k in cols_lower), None)
    y_col = next((cols_lower[k] for k in y_candidates if k in cols_lower), None)
    if mz_col is None or y_col is None:
        ordered = list(df.columns)
        if len(ordered) < 2:
            raise ValueError("DataFrame must have at least two columns for m/z and intensity.")
        mz_col, y_col = ordered[0], ordered[1]
    return mz_col, y_col

def _load_spectrum(source: Union[str, pd.DataFrame]) -> pd.DataFrame:
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
    df["mz"] = pd.to_numeric(df["mz"], errors="coerce")
    df["Intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
    df = df.dropna().sort_values("mz").reset_index(drop=True)
    return df

def _y_values(df: pd.DataFrame, normalize: bool) -> Tuple[np.ndarray, str]:
    if normalize:
        max_i = df["Intensity"].max()
        y = (df["Intensity"] / max_i) * 100.0 if max_i > 0 else np.zeros_like(df["Intensity"].values, dtype=float)
        y_label = "Relative Intensity (%)"
    else:
        y = df["Intensity"].values.astype(float)
        y_label = "Intensity"
    return y, y_label

def plot_single_mz_spectrum(
    dataset: Union[str, pd.DataFrame],
    label: str,
    second_label: Optional[str] = None,
    x_range: Optional[Tuple[float, float]] = None,
    color: str = "#2C7FB8",
    normalize: bool = True,
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False,
):
    plt.rcParams.update({"font.family": font_name, "font.size": font_size})
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

    df = _load_spectrum(dataset)
    y, y_label = _y_values(df, normalize)
    ax.plot(df["mz"], y, color=color, linewidth=1.8, label=label)

    if x_range is not None:
        ax.set_xlim(*x_range)
    else:
        xmin, xmax = float(df["mz"].min()), float(df["mz"].max())
        xr = xmax - xmin if xmax > xmin else 1.0
        ax.set_xlim(xmin - 0.02 * xr, xmax + 0.02 * xr)

    if normalize:
        ax.set_ylim(0, 100)
    else:
        ymax = float(np.nanmax(y)) if y.size else 1.0
        ax.set_ylim(0, ymax * 1.05)

    ax.set_xlabel("m/z")
    ax.set_ylabel(y_label)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper left")

    if second_label:
        peak_idx = int(np.argmax(y))
        peak_x = float(df["mz"].iloc[peak_idx])
        peak_y = float(y[peak_idx])
        x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
        xr = max(x1 - x0, 1e-9); yr = max(y1 - y0, 1e-9)
        xf = float(np.clip((peak_x - x0) / xr + 0.05, 0.02, 0.98))
        yf = float(np.clip((peak_y - y0) / yr, 0.02, 0.98))
        handle = Line2D([0], [0], color=color, lw=1.5)
        leg2 = ax.legend([handle], [second_label], loc="center", bbox_to_anchor=(xf, yf), frameon=False, fontsize=font_size)
        ax.add_artist(leg2)

    if save: fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show: plt.show()
    return fig, ax

def plot_mz_spectrum_offset_stack(
    datasets: Sequence[Union[str, pd.DataFrame]],
    labels: Sequence[str],
    second_labels: Optional[Sequence[Optional[str]]] = None,
    x_range: Optional[Tuple[float, float]] = None,
    colors: Optional[Sequence[str]] = None,
    y_offset_step: float = 120.0,
    y_tick_labels: Optional[Sequence[str]] = None,
    normalize: bool = True,
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False,
):
    n = len(datasets)
    if len(labels) != n:
        raise ValueError("labels length must match datasets length.")
    if second_labels is None:
        second_labels = [None] * n
    if len(second_labels) != n:
        raise ValueError("second_labels length must match datasets length (or None).")

    if colors is None or len(colors) < n:
        base = ["#2C7FB8", "#CC79A7", "#00A9A5", "#7A7A7A", "#F0E442"]
        k = (n + len(base) - 1) // len(base)
        colors = (base * k)[:n]
    else:
        colors = list(colors)[:n]

    if y_tick_labels is None:
        y_tick_labels = ("0", "%", "100") if normalize else ("", "", "")

    plt.rcParams.update({"font.family": font_name, "font.size": font_size})
    fig, axes = plt.subplots(n, 1, sharex=True, figsize=figsize, dpi=dpi)
    if n == 1:
        axes = [axes]

    spectra = [_load_spectrum(src) for src in datasets]

    if x_range is None:
        xmin = min(float(df["mz"].min()) for df in spectra if len(df))
        xmax = max(float(df["mz"].max()) for df in spectra if len(df))
        xr = xmax - xmin if xmax > xmin else 1.0
        x_range = (xmin - 0.02 * xr, xmax + 0.02 * xr)

    for i, (df, label, second, color) in enumerate(zip(spectra, labels, second_labels, colors)):
        ax = axes[i]
        offset = i * y_offset_step
        y, y_label = _y_values(df, normalize)

        ax.plot(df["mz"], y + offset, color=color, linewidth=1.5, label=label)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(*x_range)

        if normalize:
            ax.set_ylim(offset, offset + 101.5)
            ax.set_yticks([offset, offset + 50, offset + 100])
        else:
            ymax = float(np.nanmax(y)) if y.size else 1.0
            ax.set_ylim(offset, offset + ymax * 1.05)
            ax.set_yticks([offset, offset + ymax / 2, offset + ymax])

        ax.set_yticklabels(list(y_tick_labels))
        leg1 = ax.legend([ax.lines[-1]], [label], loc="upper left", frameon=False, fontsize=font_size)
        ax.add_artist(leg1)

        if second:
            peak_idx = int(np.argmax(y))
            peak_x = float(df["mz"].iloc[peak_idx])
            peak_y = float(y[peak_idx] + offset)
            x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
            xr = max(x1 - x0, 1e-9); yr = max(y1 - y0, 1e-9)
            xf = float(np.clip((peak_x - x0) / xr + 0.05, 0.02, 0.98))
            yf = float(np.clip((peak_y - y0) / yr, 0.02, 0.98))
            handle = Line2D([0], [0], color=color, lw=1.5)
            leg2 = ax.legend([handle], [second], loc="center", bbox_to_anchor=(xf, yf), frameon=False, fontsize=font_size)
            ax.add_artist(leg2)

        if i < n - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("m/z")

    fig.text(0.05, 0.5, y_label, va="center", rotation="vertical", fontsize=font_size)
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95, hspace=0.6)

    if save: fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show: plt.show()
    return fig, axes