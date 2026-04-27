from __future__ import annotations
from typing import Sequence, Tuple, Optional, List
import pandas as pd
import matplotlib.pyplot as plt

def read_sec_file(file_path: str, skiprows: int = 3, usecols: str = "A:B", sheet_name: int | str = 0) -> pd.DataFrame:
    """
    Read chromatogram from Excel starting at row 4 (skiprows=3), columns A:B.
    Returns DataFrame with columns ['ml', 'mAU'] and drops NaNs.
    """
    df = pd.read_excel(file_path, skiprows=skiprows, header=None, usecols=usecols, sheet_name=sheet_name, engine="openpyxl")
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
    tick_length: float = 3.0,
    tick_width: float = 0.8,
    legend_bbox: Tuple[float, float] = (1.02, 0.5),
    tight_rect: Tuple[float, float, float, float] = (0.08, 0.05, 0.95, 0.97),
    h_pad: float = 0.6,
    hspace: float = 0.55,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Stacked chromatograms with shared x-axis and optional top ticks for standards.
    """
    n = len(files)
    if len(sample_labels) != n:
        raise ValueError("sample_labels length must match files length.")

    if colors is None or len(colors) < n:
        base = ["#333333", "#00A9A5", "#C76BA0", "#2C7FB8", "#7A7A7A", "#00A9A5", "#F0E442", "#56B4E9"]
        k = (n + len(base) - 1) // len(base)
        colors = (base * k)[:n]
    else:
        colors = list(colors)[:n]

    plt.rcParams.update({"font.family": font_name, "font.size": font_size})

    fig, axes = plt.subplots(n, 1, figsize=figsize, sharex=True, dpi=dpi)
    if n == 1:
        axes = [axes]

    for i, (file, label) in enumerate(zip(files, sample_labels)):
        df = read_sec_file(file)
        if normalize:
            max_val = df["mAU"].max()
            if max_val > 0:
                df["mAU"] = df["mAU"] / max_val

        ax = axes[i]
        ax.plot(df["ml"], df["mAU"], color=colors[i], linewidth=1.8, label=label)
        ax.set_ylim(*ylim)
        ax.set_xlim(*xlim)
        ax.margins(y=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_yticks([ylim[0], ylim[1]])
        ax.tick_params(axis="y", length=tick_length, width=tick_width)
        ax.tick_params(axis="x", length=tick_length, width=tick_width)
        ax.legend(loc="center left", bbox_to_anchor=legend_bbox, frameon=False, handlelength=1.5)

        if i == 0 and standard_volumes:
            top_ax = ax.secondary_xaxis("top")
            top_ax.set_xticks(list(standard_volumes))
            top_ax.set_xticklabels([])
            top_ax.tick_params(axis="x", length=4, width=tick_width)
            top_ax.spines["top"].set_visible(False)

        if i < n - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Retention Volume (mL)")

    fig.text(0.04, 0.5, "Normalized Absorbance" if normalize else "Absorbance (mAU)", va="center", rotation="vertical", fontsize=font_size)
    plt.tight_layout(rect=tight_rect, h_pad=h_pad)
    plt.subplots_adjust(hspace=hspace)

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes