# cdms_mass_plot_utils.py
from __future__ import annotations
from typing import Sequence, Tuple, Optional, Union, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def _load_mass(source: Union[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Load or validate mass distribution data with columns: Mass_Da, Intensity.
    Adds Mass_kDa and Percent.
    """
    if isinstance(source, str):
        df = pd.read_csv(source, sep=r"\s+", header=None)
        df.columns = ["Mass_Da", "Intensity"]
    elif isinstance(source, pd.DataFrame):
        df = source.copy()
        if not {"Mass_Da", "Intensity"}.issubset(df.columns):
            raise ValueError("DataFrame must have 'Mass_Da' and 'Intensity'")
    else:
        raise TypeError("source must be str or pd.DataFrame.")

    df["Mass_kDa"] = df["Mass_Da"] / 1000.0
    max_i = df["Intensity"].max()
    df["Percent"] = (df["Intensity"] / max_i) * 100.0 if max_i > 0 else 0.0
    return df

def plot_single_cdms_mass(
    dataset: Union[str, pd.DataFrame],
    label: str,
    second_label: Optional[str] = None,
    x_range: Tuple[float, float] = (0, 500),
    color: str = "#2C7FB8",
    y_labels: Sequence[str] = ("0", "%", "100"),
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False
):
    """
    Plot a single mass distribution with optional second legend near the peak.
    """
    plt.rcParams.update({"font.family": font_name, "font.size": font_size})
    fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

    df = _load_mass(dataset)
    ax.plot(df["Mass_kDa"], df["Percent"], color=color, linewidth=1.8, label=label)
    ax.set_xlim(*x_range)
    ax.set_xlabel("Mass (kDa)")
    ax.set_ylabel("Relative Intensity (%)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper left")

    if second_label:
        peak_idx = int(df["Percent"].idxmax())
        peak_x = float(df["Mass_kDa"].iloc[peak_idx])
        peak_y = float(df["Percent"].iloc[peak_idx])
        # place second legend near peak
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        xr = max(x1 - x0, 1e-9)
        yr = max(y1 - y0, 1e-9)
        xf = (peak_x - x0) / xr
        yf = (peak_y - y0) / yr
        xf = float(np.clip(xf + 0.05, 0.02, 0.98))
        yf = float(np.clip(yf, 0.02, 0.98))
        handle = Line2D([0], [0], color=color, lw=1.5)
        leg2 = ax.legend([handle], [second_label], loc="center", bbox_to_anchor=(xf, yf), frameon=False, fontsize=font_size)
        ax.add_artist(leg2)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax

def plot_cdms_mass_offset_stack(
    datasets: Sequence[Union[str, pd.DataFrame]],
    labels: Sequence[str],
    second_labels: Optional[Sequence[Optional[str]]] = None,
    x_range: Tuple[float, float] = (0, 500),
    colors: Optional[Sequence[str]] = None,
    y_offset_step: float = 120.0,
    y_labels: Sequence[str] = ("0", "%", "100"),
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False
):
    """
    Stacked mass distributions with vertical offsets and optional second legend near peaks.
    """
    n = len(datasets)
    if len(labels) != n:
        raise ValueError("labels length must match datasets length.")
    if second_labels is None:
        second_labels = [None] * n
    if len(second_labels) != n:
        raise ValueError("second_labels length must match datasets length (or None).")

    if colors is None or len(colors) < n:
        base = ["#CC79A7", "#CC79A7", "#2C7FB8", "#7A7A7A", "#00A9A5", "#F0E442"]
        k = (n + len(base) - 1) // len(base)
        colors = (base * k)[:n]
    else:
        colors = list(colors)[:n]

    plt.rcParams.update({"font.family": font_name, "font.size": font_size})
    fig, axes = plt.subplots(n, 1, sharex=True, figsize=figsize, dpi=dpi)
    if n == 1:
        axes = [axes]

    for i, (src, label, second) in enumerate(zip(datasets, labels, second_labels)):
        df = _load_mass(src)
        ax = axes[i]
        offset = i * y_offset_step
        color = colors[i]
        ax.plot(df["Mass_kDa"], df["Percent"] + offset, color=color, linewidth=1.5, label=label)
        # style
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(*x_range)
        ax.set_ylim(offset, offset + 101.5)
        ax.set_yticks([offset, offset + 50, offset + 100])
        ax.set_yticklabels(list(y_labels))
        # legends
        leg1 = ax.legend([ax.lines[-1]], [label], loc="upper left", frameon=False, fontsize=font_size)
        ax.add_artist(leg1)

        if second:
            peak_idx = int(df["Percent"].idxmax())
            peak_x = float(df["Mass_kDa"].iloc[peak_idx])
            peak_y = float(df["Percent"].iloc[peak_idx] + offset)
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            xr = max(x1 - x0, 1e-9)
            yr = max(y1 - y0, 1e-9)
            xf = (peak_x - x0) / xr
            yf = (peak_y - y0) / yr
            xf = float(np.clip(xf + 0.05, 0.02, 0.98))
            yf = float(np.clip(yf, 0.02, 0.98))
            handle = Line2D([0], [0], color=color, lw=1.5)
            leg2 = ax.legend([handle], [second], loc="center", bbox_to_anchor=(xf, yf), frameon=False, fontsize=font_size)
            ax.add_artist(leg2)

        if i < n - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Mass (kDa)")

    fig.text(0.05, 0.5, "Relative Intensity (%)", va="center", rotation="vertical", fontsize=font_size)
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95, hspace=0.6)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes