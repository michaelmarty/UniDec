from __future__ import annotations
from typing import Sequence, Tuple, Optional, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from plot_helpers import style_context, ensure_colors, compute_target_height, DEFAULT_CHROM_COLORS

DEFAULT_WIDTH_IN = 10.0
ROW_HEIGHT_IN = 1.35
MIN_HEIGHT_IN = 3.8


def read_chrom_file(
    file_path: str,
    skiprows: int = 3,
    usecols: str = "A:B",
    sheet_name: Union[int, str] = 0
) -> pd.DataFrame:
    """
    Read a chromatogram from Excel (default) or text.

    Returns a DataFrame with columns:
        - 'ml'  : x-axis (volume or time)
        - 'mAU' : y-axis (signal)
    """
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
        # Fallback: text file with common delimiters; first col = x, second col = y
        seps = [r"\s+", "\t", ","]  # whitespace, tab, comma
        df = None

        for sep in seps:
            try:
                df_try = pd.read_csv(
                    file_path,
                    sep=sep,
                    header=None,
                    engine="python",
                    comment="#",
                )
                if df_try.shape[1] >= 2:
                    df = df_try
                    break
            except Exception:
                continue

        if df is None:
            df = pd.read_csv(file_path, header=None, engine="python")

        if df.shape[1] < 2:
            raise ValueError(
                f"Chromatogram text file '{file_path}' must have at least two columns "
                f"(first = x, second = y). Found only {df.shape[1]} column(s)."
            )

        df2 = df.iloc[:, :2].copy()
        df2.columns = ["ml", "mAU"]
        df2["ml"] = pd.to_numeric(df2["ml"], errors="coerce")
        df2["mAU"] = pd.to_numeric(df2["mAU"], errors="coerce")
        df2 = df2.dropna()

        return df2


def _plain_y_formatter() -> ScalarFormatter:
    """
    Disable scientific notation and use plain decimal for y-axis.
    """
    fmt = ScalarFormatter(useMathText=False)
    fmt.set_scientific(False)
    fmt.set_useOffset(False)
    fmt.set_powerlimits((0, 0))
    return fmt


def plot_chrom_stack(
    files: Sequence[str],
    sample_labels: Sequence[str],
    colors: Optional[Sequence[str]] = None,
    standard_volumes: Optional[Sequence[float]] = None,
    normalize: bool = True,
    norm_mode: str = "percent",  # "percent" or "unit"
    font_name: str = "Arial",
    font_size: int = 10,
    xlim: Tuple[float, float] = (10.0, 20.0),
    ylim: Tuple[float, float] = (0.0, 1.0),  # interpreted only when normalize=True
    figsize: Tuple[float, float] = (8, 10),
    dpi: int = 300,
    tick_length: float = 5.0,
    tick_width: float = 1.0,
    legend_bbox: Tuple[float, float] = (0.02, 0.5),
    tight_rect: Tuple[float, float, float, float] = (0.08, 0.08, 0.98, 0.96),
    h_pad: float = 0.6,
    hspace: float = 0.55,
    line_width: Optional[float] = None,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Multi-panel chromatogram plot (one subplot per file).

    - If normalize=True:
        * Traces are normalized (0–100% or 0–1).
        * ylim=(ymin,ymax) is interpreted in normalized units and applied.
    - If normalize=False:
        * Raw mAU is used.
        * ylim is ignored; each panel gets [0, max(mAU)] from data.
    """
    n = len(files)
    if len(sample_labels) != n:
        raise ValueError("sample_labels length must match files length.")

    colors = ensure_colors(colors, n, DEFAULT_CHROM_COLORS)
    if line_width is None:
        line_width = 1.8

    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = compute_target_height(n_rows=n)

    with style_context(font_name=font_name, font_size=font_size, tick_length=tick_length, tick_width=tick_width):
        fig, axes = plt.subplots(n, 1, figsize=(width_in, target_h), sharex=True, dpi=dpi)
        if n == 1:
            axes = [axes]

        for i, (fp, label) in enumerate(zip(files, sample_labels)):
            df = read_chrom_file(fp)
            df_local = df.copy()

            # Normalization
            if normalize:
                max_val = float(df_local["mAU"].max()) if not df_local["mAU"].empty else 0.0
                if max_val > 0:
                    if norm_mode == "percent":
                        df_local["mAU"] = (df_local["mAU"] / max_val) * 100.0
                    else:
                        df_local["mAU"] = df_local["mAU"] / max_val

            ax = axes[i]
            ax.plot(df_local["ml"], df_local["mAU"], color=colors[i], linewidth=line_width, label=label)
            ax.set_xlim(*xlim)

            # Y-limits & ticks
            if normalize:
                ymin_req, ymax_req = ylim

                if norm_mode == "percent":
                    clamp_lo, clamp_hi = 0.0, 100.0
                else:
                    clamp_lo, clamp_hi = 0.0, 1.0

                ymin = max(clamp_lo, ymin_req)
                ymax = min(clamp_hi, ymax_req)
                if ymax <= ymin:
                    # Fallback if user supplies invalid range
                    ymin, ymax = clamp_lo, clamp_hi

                ax.set_ylim(ymin, ymax)

                mid = ymin + 0.5 * (ymax - ymin)
                if norm_mode == "percent":
                    ax.set_yticks([ymin, mid, ymax])
                    ax.set_yticklabels([
                        f"{ymin:.0f}",
                        "%",
                        f"{ymax:.0f}",
                    ])
                else:
                    ax.set_yticks([ymin, mid, ymax])
                    ax.set_yticklabels([
                        f"{ymin:.2g}",
                        f"{mid:.2g}",
                        f"{ymax:.2g}",
                    ])
            else:
                # Raw mAU
                y_series = df_local["mAU"]
                if y_series.empty:
                    ymin, ymax = 0.0, 1.0
                else:
                    ymin = 0.0
                    ymax = float(y_series.max())
                    if ymax <= 0:
                        ymax = 1.0

                ax.set_ylim(ymin, ymax)
                ax.set_yticks([ymin, ymax])
                ax.yaxis.set_major_formatter(_plain_y_formatter())

            ax.margins(y=0)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

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

        fig.supxlabel(x_label if x_label else "Volume (mL)")
        if y_label:
            fig.supylabel(y_label)
        else:
            fig.supylabel("Normalized Absorbance" if normalize else "Absorbance (mAU)")

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
    line_width: Optional[float] = None,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Offset-stacked chromatograms: all traces share a single x-axis but are vertically offset.
    Not used by the main GUI stack mode, but kept for completeness.
    """
    n = len(files)
    if len(sample_labels) != n:
        raise ValueError("sample_labels length must match files length.")

    colors = ensure_colors(colors, n, DEFAULT_CHROM_COLORS)
    if line_width is None:
        line_width = 1.8

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

        for i, (df_local, label, color) in enumerate(zip(dfs, sample_labels, colors)):
            offset = i * y_offset_step
            x = df_local["ml"].values
            y = df_local["mAU"].values + offset
            ax.plot(x, y, color=color, linewidth=line_width, label=label)

        ax.set_xlim(*xlim)
        ax.margins(y=0)

        if normalize:
            ymin = 0.0
            ymax = (n - 1) * y_offset_step + 1.0
            ax.set_ylim(ymin, ymax)
        else:
            all_y = []
            for i, df_local in enumerate(dfs):
                offset = i * y_offset_step
                all_y.append(df_local["mAU"].values + offset)
            if all_y:
                all_y_concat = np.concatenate(all_y)
                ymin = float(np.nanmin(all_y_concat))
                ymax = float(np.nanmax(all_y_concat))
                ax.set_ylim(ymin, ymax)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if normalize:
            yticks = []
            ylabels = []
            for i in range(n):
                base = i * y_offset_step
                yticks.extend([base, base + 0.5 * y_offset_step, base + y_offset_step])
                ylabels.extend(["0", "%", "100"])
            ax.set_yticks(yticks)
            ax.set_yticklabels(ylabels)
        else:
            ymin, ymax = ax.get_ylim()
            ax.set_yticks([ymin, ymax])

        ax.tick_params(axis="y", length=tick_length, width=tick_width)
        ax.tick_params(axis="x", length=tick_length, width=tick_width)

        ax.legend(loc="center left", bbox_to_anchor=legend_bbox)

        if standard_volumes:
            top_ax = ax.secondary_xaxis("top")
            top_ax.set_xticks(list(standard_volumes))
            top_ax.set_xticklabels([])
            top_ax.tick_params(axis="x", length=tick_length - 1, width=tick_width)
            top_ax.spines["top"].set_visible(False)

        fig.supxlabel(x_label if x_label else "Volume (mL)")
        fig.supylabel(
            y_label if y_label
            else ("Normalized Absorbance (offset)" if normalize else "Absorbance (mAU, offset)")
        )

        plt.tight_layout(rect=tight_rect, h_pad=h_pad)
        plt.subplots_adjust(hspace=hspace)

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, ax