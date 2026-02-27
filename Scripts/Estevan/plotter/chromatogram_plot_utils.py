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
        # Try to read plain text chromatogram with common delimiters.
        # We assume: first column = x (ml), second column = y (mAU).
        seps = [r"\s+", "\t", ","]  # whitespace, tab, comma
        df = None

        for sep in seps:
            try:
                df_try = pd.read_csv(
                    file_path,
                    sep=sep,
                    header=None,
                    engine="python",
                    comment="#",  # ignore commented lines
                )
                # Must have at least two columns of data
                if df_try.shape[1] >= 2:
                    df = df_try
                    break
            except Exception:
                continue

        # If still None, fall back to generic read_csv (whatever delimiter)
        if df is None:
            df = pd.read_csv(file_path, header=None, engine="python")

        if df.shape[1] < 2:
            raise ValueError(
                f"Chromatogram text file '{file_path}' must have at least two columns "
                f"(first = x, second = y). Found only {df.shape[1]} column(s)."
            )

        # Take first two columns as x and y
        df2 = df.iloc[:, :2].copy()
        df2.columns = ["ml", "mAU"]

        # Coerce to numeric and drop NaNs
        df2["ml"] = pd.to_numeric(df2["ml"], errors="coerce")
        df2["mAU"] = pd.to_numeric(df2["mAU"], errors="coerce")
        df2 = df2.dropna()

        return df2


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

def plot_chrom_offset_stack(
    files: Sequence[str],
    sample_labels: Sequence[str],
    colors: Optional[Sequence[str]] = None,
    standard_volumes: Optional[Sequence[float]] = None,
    normalize: bool = True,
    font_name: str = "Arial",
    font_size: int = 10,
    xlim: Tuple[float, float] = (10.0, 20.0),
    y_offset_step: float = 1.0,
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
    """
    Offset-stacked chromatograms: all traces share a single x-axis but are vertically offset
    by y_offset_step. Normalization is applied before offsetting if normalize=True.
    """
    n = len(files)
    if len(sample_labels) != n:
        raise ValueError("sample_labels length must match files length.")

    colors = ensure_colors(colors, n, DEFAULT_CHROM_COLORS)

    # Prepare data
    dfs = []
    for fp in files:
        df = read_chrom_file(fp)
        df_local = df.copy()
        if normalize:
            max_val = float(df_local["mAU"].max()) if not df_local["mAU"].empty else 0.0
            if max_val > 0:
                df_local["mAU"] = df_local["mAU"] / max_val
        dfs.append(df_local)

    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = compute_target_height(n_rows=n)

    with style_context(font_name=font_name, font_size=font_size, tick_length=tick_length, tick_width=tick_width):
        fig, ax = plt.subplots(1, 1, figsize=(width_in, target_h), dpi=dpi)

        # Plot each chromatogram with vertical offset
        for i, (df_local, label, color) in enumerate(zip(dfs, sample_labels, colors)):
            offset = i * y_offset_step
            x = df_local["ml"].values
            y = df_local["mAU"].values + offset

            ax.plot(x, y, color=color, linewidth=1.8, label=label)

        # X limits
        ax.set_xlim(*xlim)
        ax.margins(y=0)

        # Y limits
        if normalize:
            # Each trace is 0–1, shifted by offset * i; last trace top ≈ (n-1)*offset + 1
            ymin = 0.0
            ymax = (n - 1) * y_offset_step + 1.0
            ax.set_ylim(ymin, ymax)
        else:
            # Not normalized: derive bounds from actual data
            all_y = []
            for i, df_local in enumerate(dfs):
                offset = i * y_offset_step
                all_y.append(df_local["mAU"].values + offset)
            if all_y:
                all_y_concat = np.concatenate(all_y)
                ymin = float(np.nanmin(all_y_concat))
                ymax = float(np.nanmax(all_y_concat))
                ax.set_ylim(ymin, ymax)

        # Style
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # For normalized chromatograms, show generic 0 / mid / top ticks relative to first trace
        if normalize:
            # Define ticks for each offset tier, like CD-MS / m/z
            yticks = []
            ylabels = []
            for i in range(n):
                base = i * y_offset_step
                yticks.extend([base, base + 0.5 * y_offset_step, base + y_offset_step])
                ylabels.extend(["0", "%", "100"])
            # To avoid overly dense labels, you might want only the last tier;
            # here we keep it simple and label each tier the same.
            ax.set_yticks(yticks)
            ax.set_yticklabels(ylabels)
        else:
            ymin, ymax = ax.get_ylim()
            ax.set_yticks([ymin, ymax])

        ax.tick_params(axis="y", length=tick_length, width=tick_width)
        ax.tick_params(axis="x", length=tick_length, width=tick_width)

        # Legends: show one legend per line, stacked along left side
        # Simpler: use a single legend listing all labels
        ax.legend(loc="center left", bbox_to_anchor=legend_bbox)

        # Top ticks for standard volumes
        if standard_volumes:
            top_ax = ax.secondary_xaxis("top")
            top_ax.set_xticks(list(standard_volumes))
            top_ax.set_xticklabels([])
            top_ax.tick_params(axis="x", length=tick_length - 1, width=tick_width)
            top_ax.spines["top"].set_visible(False)

        fig.supxlabel("Volume (mL)")
        fig.supylabel("Normalized Absorbance (offset)" if normalize else "Absorbance (mAU, offset)")

        plt.tight_layout(rect=tight_rect, h_pad=h_pad)
        plt.subplots_adjust(hspace=hspace)

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, ax