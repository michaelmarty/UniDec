# chromatogram_plot_utils.py
from __future__ import annotations
from typing import Sequence, Tuple, Optional, List, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Consistent sizing for stacked plots
DEFAULT_WIDTH_IN = 10.0     # used when only height is recomputed
ROW_HEIGHT_IN     = 1.35     # per subplot
MIN_HEIGHT_IN     = 3.8      # minimum total height


def read_chrom_file(
    file_path: str,
    skiprows: int = 3,
    usecols: str = "A:B",
    sheet_name: Union[int, str] = 0
) -> pd.DataFrame:
    """
    Read an Excel chromatogram (A:B, starting at row 4).
    Returns DataFrame with columns ['ml', 'mAU'] and drops NaNs.
    """
    df = pd.read_excel(
        file_path,
        skiprows=skiprows,
        header=None,
        usecols=usecols,
        sheet_name=sheet_name,
        engine="openpyxl",
    )
    df.columns = ["ml", "mAU"]
    df = df.dropna()
    return df


def plot_chrom_stack(
    files: Sequence[str],
    sample_labels: Sequence[str],
    colors: Optional[Sequence[str]] = None,
    standard_volumes: Optional[Sequence[float]] = None,
    normalize: bool = True,
    font_name: str = "Arial",
    font_size: int = 10,
    xlim: Tuple[float, float] = (10.0, 20.0),
    ylim: Tuple[float, float] = (0.0, 1.0),
    figsize: Tuple[float, float] = (8, 10),
    dpi: int = 300,
    tick_length: float = 5.0,
    tick_width: float = 1.0,
    legend_bbox: Tuple[float, float] = (0.02, 0.5),
    tight_rect: Tuple[float, float, float, float] = (0.08, 0.08, 0.98, 0.96),
    h_pad: float = 0.6,
    hspace: float = 0.55,
    save: Optional[str] = None,
    show: bool = False,
):
    n = len(files)
    if len(sample_labels) != n:
        raise ValueError("sample_labels length must match files length.")

    if colors is None or len(colors) < n:
        base = ["#333333", "#00A9A5", "#C76BA0", "#2C7FB8", "#7A7A7A", "#00A9A5", "#F0E442", "#56B4E9"]
        k = (n + len(base) - 1) // len(base)
        colors = (base * k)[:n]
    else:
        colors = list(colors)[:n]

    plt.rcParams.update({
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

    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = max(MIN_HEIGHT_IN, (ROW_HEIGHT_IN * max(1, n)) + 0.8)

    fig, axes = plt.subplots(n, 1, figsize=(width_in, target_h), sharex=True, dpi=dpi)
    if n == 1:
        axes = [axes]

    for i, (file, label) in enumerate(zip(files, sample_labels)):
        df = read_chrom_file(file)

        # Normalize (kept in 0-1 space)
        if normalize:
            max_val = float(df["mAU"].max())
            if max_val > 0:
                df["mAU"] = df["mAU"] / max_val

        ax = axes[i]
        ax.plot(df["ml"], df["mAU"], color=colors[i], linewidth=1.8, label=label)

        # X limits
        ax.set_xlim(*xlim)

        # Y limits: auto-scale for non-normalized if GUI default (0..1) was passed
        if normalize:
            ax.set_ylim(0.0, 1.0)
        else:
            if ylim == (0.0, 1.0):
                ymax = float(df["mAU"].max())
                ax.set_ylim(0.0, max(1.0, ymax))
            else:
                ax.set_ylim(*ylim)

        ax.margins(y=0)

        # Clean spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Y ticks
        if normalize:
            ax.set_yticks([0.0, 0.5, 1.0])
            ax.set_yticklabels(["0", "%", "100"])
        else:
            ymin, ymax = ax.get_ylim()
            ax.set_yticks([ymin, ymax])

        # Tick styling
        ax.tick_params(axis="y", length=tick_length, width=tick_width)
        ax.tick_params(axis="x", length=tick_length, width=tick_width)

        # Legend: frameless, left-center
        ax.legend(loc="center left", bbox_to_anchor=legend_bbox)

        # Standards on top axis (first panel only)
        if i == 0 and standard_volumes:
            top_ax = ax.secondary_xaxis("top")
            top_ax.set_xticks(list(standard_volumes))
            top_ax.set_xticklabels([])
            top_ax.tick_params(axis="x", length=tick_length - 1, width=tick_width)
            top_ax.spines["top"].set_visible(False)

        # Hide x tick labels on intermediate rows
        if i < n - 1:
            ax.tick_params(labelbottom=False)

    # Shared labels
    fig.supxlabel("Volume (mL)")
    fig.supylabel("Normalized Absorbance (%)" if normalize else "Absorbance (mAU)")

    # Layout
    plt.tight_layout(rect=tight_rect, h_pad=h_pad)
    plt.subplots_adjust(hspace=hspace)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()

    return fig, axes