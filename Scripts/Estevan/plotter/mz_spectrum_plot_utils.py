from __future__ import annotations
from typing import Sequence, Tuple, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot_helpers import style_context, place_second_label, ensure_colors, compute_target_height, DEFAULT_MZ_COLORS

DEFAULT_WIDTH_IN = 10.0
ROW_HEIGHT_IN = 1.35
MIN_HEIGHT_IN = 3.8


def read_mz_spectrum_file(file_path: str) -> pd.DataFrame:
    for sep, header in [(None, None), (r"\s+", None), (",", None), (",", 0), (r"\s+", 0)]:
        try:
            df = pd.read_csv(file_path, sep=sep, header=header, engine="python")
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

    arr = np.loadtxt(file_path)
    if arr.ndim == 1:
        arr = np.column_stack([np.arange(len(arr)), arr])
    df = pd.DataFrame(arr[:, :2], columns=["mz", "intensity"]).dropna()
    return df


def plot_single_mz_spectrum(
    dataset: str,
    label: str,
    second_label: Optional[str] = None,
    x_range: Tuple[float, float] = (100.0, 2000.0),
    color: str = "#2C7FB8",
    normalize: bool = True,
    norm_mode: str = "percent",  # "percent" or "unit"
    y_max: Optional[float] = None,
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (10, 4),
    dpi: int = 300,
    tick_length: float = 5.0,
    tick_width: float = 1.0,
    line_width: Optional[float] = None,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
):
    if line_width is None:
        line_width = 2.0

    rc = {
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

    df = read_mz_spectrum_file(dataset)
    if df is None or df.empty:
        raise RuntimeError(f"No data found in spectrum file: {dataset}")

    df = df[(df["mz"] >= x_range[0]) & (df["mz"] <= x_range[1])].copy()
    if len(df) == 0:
        df = read_mz_spectrum_file(dataset)

    x = df["mz"].astype(float).values
    y = df["intensity"].astype(float).values

    if normalize:
        m = float(np.nanmax(y)) if y.size else 0.0
        if m > 0:
            if norm_mode == "percent":
                y = (y / m) * 100.0
            else:
                y = y / m

    with plt.rc_context(rc):
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

        ax.plot(x, y, color=color, linewidth=line_width, label=label)

        ax.set_xlim(*x_range)
        ax.margins(y=0)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if normalize:
            if norm_mode == "percent":
                default_ymax = 100.0
            else:
                default_ymax = 1.0

            if y_max is None:
                ymax = default_ymax
            else:
                if norm_mode == "percent":
                    ymax = max(0.0, min(float(y_max), 100.0))
                else:
                    ymax = max(0.0, min(float(y_max), 1.0))

            ax.set_ylim(0.0, ymax)
            mid = ymax / 2.0
            if norm_mode == "percent":
                ax.set_yticks([0.0, mid, ymax])
                ax.set_yticklabels(["0", "%", f"{ymax:.0f}"])
            else:
                ax.set_yticks([0.0, mid, ymax])
                ax.set_yticklabels(["0", f"{mid:.2g}", f"{ymax:.2g}"])
            ax.set_ylabel(y_label if y_label else "Normalized Intensity")
        else:
            ax.set_ylabel(y_label if y_label else "Intensity (a.u.)")

        ax.set_xlabel(x_label if x_label else "m/z")
        ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.5))

        place_second_label(ax, x, y, second_label, color, x_range, font_size)

        plt.tight_layout()

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, ax


def plot_mz_spectrum_stack(
    datasets: Sequence[str],
    labels: Sequence[str],
    second_labels: Sequence[Optional[str]],
    x_range: Tuple[float, float] = (100.0, 2000.0),
    colors: Optional[Sequence[str]] = None,
    y_max: Optional[float] = None,
    normalize: bool = True,
    norm_mode: str = "percent",  # "percent" or "unit"
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (10, 8),
    dpi: int = 300,
    hspace: float = 0.6,
    tick_length: float = 5.0,
    tick_width: float = 1.0,
    line_width: Optional[float] = None,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Stacked m/z spectra: one panel per dataset, no vertical offsets.
    """
    n = len(datasets)
    if len(labels) != n:
        raise ValueError("labels length must match datasets length.")
    if second_labels is None or len(second_labels) != n:
        second_labels = [None] * n

    colors = ensure_colors(colors, n, DEFAULT_MZ_COLORS)
    if line_width is None:
        line_width = 2.0

    rc = {
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

    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = compute_target_height(n_rows=n)

    if normalize:
        if norm_mode == "percent":
            default_ymax = 100.0
        else:
            default_ymax = 1.0

        if y_max is None:
            norm_ymax = default_ymax
        else:
            if norm_mode == "percent":
                norm_ymax = max(0.0, min(float(y_max), 100.0))
            else:
                norm_ymax = max(0.0, min(float(y_max), 1.0))
    else:
        norm_ymax = None

    with plt.rc_context(rc):
        fig, axes = plt.subplots(n, 1, figsize=(width_in, target_h), sharex=True, dpi=dpi)
        if n == 1:
            axes = [axes]

        for i, (path, lab, sec_lab, color) in enumerate(zip(datasets, labels, second_labels, colors)):
            df = read_mz_spectrum_file(path)
            if df is None or df.empty:
                raise RuntimeError(f"No data found in spectrum file: {path}")

            df = df[(df["mz"] >= x_range[0]) & (df["mz"] <= x_range[1])].copy()
            if len(df) == 0:
                df = read_mz_spectrum_file(path)

            x = df["mz"].astype(float).values
            y_raw = df["intensity"].astype(float).values

            if normalize:
                m = float(np.nanmax(y_raw)) if y_raw.size else 0.0
                if m > 0:
                    if norm_mode == "percent":
                        y = (y_raw / m) * 100.0
                    else:
                        y = y_raw / m
                else:
                    y = np.zeros_like(y_raw)
            else:
                y = y_raw

            ax = axes[i]

            ax.plot(x, y, color=color, linewidth=line_width, label=lab)

            ax.set_xlim(*x_range)
            if normalize:
                ymax_local = norm_ymax
                ax.set_ylim(0.0, ymax_local * 1.015)
                mid = ymax_local / 2.0
                if norm_mode == "percent":
                    ax.set_yticks([0.0, mid, ymax_local])
                    ax.set_yticklabels(["0", "%", f"{ymax_local:.0f}"])
                else:
                    ax.set_yticks([0.0, mid, ymax_local])
                    ax.set_yticklabels([
                        "0",
                        f"{(ymax_local / 2.0):.2g}",
                        f"{ymax_local:.2g}",
                    ])
            else:
                if y.size:
                    ymin = float(np.nanmin(y))
                    ymax = float(np.nanmax(y))
                    ax.set_ylim(max(0.0, ymin), ymax * 1.05)
                    ax.set_yticks([0.0, ymax])

            ax.margins(y=0)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            ax.tick_params(axis="y", length=tick_length, width=tick_width)
            ax.tick_params(axis="x", length=tick_length, width=tick_width)

            ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.5))

            place_second_label(ax, x, y, sec_lab, color, x_range, font_size)

            if i < n - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(x_label if x_label else "m/z")

        fig.supxlabel(x_label if x_label else "m/z")
        if normalize:
            fig.supylabel(y_label if y_label else "Normalized Intensity")
        else:
            fig.supylabel(y_label if y_label else "Intensity (a.u.)")

        plt.tight_layout(rect=(0.08, 0.08, 0.98, 0.96), h_pad=hspace)
        plt.subplots_adjust(hspace=hspace)

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, axes