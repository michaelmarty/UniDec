from __future__ import annotations
from typing import Sequence, Tuple, Optional, List, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot_helpers import style_context, ensure_colors, compute_target_height, DEFAULT_CHROM_COLORS

DEFAULT_WIDTH_IN = 10.0
ROW_HEIGHT_IN     = 1.35
MIN_HEIGHT_IN     = 3.8


def read_chrom_file(
    file_path: str,
    skiprows: int = 3,
    usecols: str = "A:B",
    sheet_name: Union[int, str] = 0
) -> pd.DataFrame:
    try:
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
    except Exception:
        df = pd.read_csv(file_path, header=None, engine="python")
        num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if len(num_cols) >= 2:
            df = df[num_cols[:2]].copy()
        else:
            df = df.iloc[:, :2].copy()
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

    colors = ensure_colors(colors, n, DEFAULT_CHROM_COLORS)

    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = compute_target_height(n_rows=n)

    with style_context(font_name=font_name, font_size=font_size, tick_length=tick_length, tick_width=tick_width):
        fig, axes = plt.subplots(n, 1, figsize=(width_in, target_h), sharex=True, dpi=dpi)
        if n == 1:
            axes = [axes]

        for i, (fp, label) in enumerate(zip(files, sample_labels)):
            df = read_chrom_file(fp)
            df_local = df.copy()

            if normalize:
                max_val = float(df_local["mAU"].max()) if not df_local["mAU"].empty else 0.0
                if max_val > 0:
                    df_local["mAU"] = df_local["mAU"] / max_val

            ax = axes[i]
            ax.plot(df_local["ml"], df_local["mAU"], color=colors[i], linewidth=1.8, label=label)

            ax.set_xlim(*xlim)
            if normalize:
                ax.set_ylim(0.0, 1.0)
            else:
                if ylim == (0.0, 1.0):
                    ymax = float(df_local["mAU"].max()) if not df_local["mAU"].empty else 1.0
                    ax.set_ylim(0.0, max(1.0, ymax))
                else:
                    ax.set_ylim(*ylim)

            ax.margins(y=0)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if normalize:
                ax.set_yticks([0.0, 0.5, 1.0])
                ax.set_yticklabels(["0", "%", "100"])
            else:
                ymin, ymax = ax.get_ylim()
                ax.set_yticks([ymin, ymax])

            ax.tick_params(axis="y", length=tick_length, width=tick_width)
            ax.tick_params(axis="x", length=tick_length, width=tick_width)

            ax.legend(loc="center left", bbox_to_anchor=legend_bbox)

            if i == 0 and standard_volumes:
                top_ax = ax.secondary_xaxis("top")
                top_ax.set_xticks(list(standard_volumes))
                top_ax.set_xticklabels([])
                top_ax.tick_params(axis="x", length=tick_length - 1, width=tick_width)
                top_ax.spines["top"].set_visible(False)

            if i < n - 1:
                ax.tick_params(labelbottom=False)

        fig.supxlabel("Volume (mL)")
        fig.supylabel("Normalized Absorbance (%)" if normalize else "Absorbance (mAU)")

        plt.tight_layout(rect=tight_rect, h_pad=h_pad)
        plt.subplots_adjust(hspace=hspace)

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, axes