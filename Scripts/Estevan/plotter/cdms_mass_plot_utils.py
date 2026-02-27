from __future__ import annotations
from typing import Sequence, Tuple, Optional, Union, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot_helpers import (
    style_context, ensure_colors, place_second_label,
    compute_target_height, DEFAULT_CDMS_COLORS
)


def read_cdms_mass_file(file_path: str) -> pd.DataFrame:
    for sep, header in [(None, None), (r"\s+", None), (",", None), (",", 0), (r"\s+", 0)]:
        try:
            df = pd.read_csv(file_path, sep=sep, header=header, engine="python")
            num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            if len(num_cols) >= 2:
                df = df[num_cols[:2]].copy()
            else:
                df = df.iloc[:, :2].copy()
            df.columns = ["mass_kDa", "intensity"]
            df["mass_kDa"] = pd.to_numeric(df["mass_kDa"], errors="coerce")
            df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
            df = df.dropna()
            if len(df) > 2:
                return df
        except Exception:
            continue

    arr = np.loadtxt(file_path)
    if arr.ndim == 1:
        arr = np.column_stack([np.arange(len(arr)), arr])
    return pd.DataFrame(arr[:, :2], columns=["mass_kDa", "intensity"]).dropna()


def _load_mass(source: Union[str, pd.DataFrame]) -> pd.DataFrame:
    if isinstance(source, str):
        try:
            df = read_cdms_mass_file(source)
            return df.sort_values("mass_kDa").reset_index(drop=True)
        except Exception:
            df = pd.read_csv(source, sep=r"\s+", header=None)
            df = df.iloc[:, :2].copy()
            df.columns = ["mass_kDa", "intensity"]
    elif isinstance(source, pd.DataFrame):
        df = source.copy()
        col_map = {}
        if {"mass_kDa", "intensity"}.issubset(df.columns):
            pass
        elif {"Mass_kDa", "Intensity"}.issubset(df.columns):
            col_map = {"Mass_kDa": "mass_kDa", "Intensity": "intensity"}
        elif {"Mass_Da", "Intensity"}.issubset(df.columns):
            df["mass_kDa"] = pd.to_numeric(df["Mass_Da"], errors="coerce") / 1000.0
            df["intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
            df = df[["mass_kDa", "intensity"]]
            return df.dropna().sort_values("mass_kDa").reset_index(drop=True)
        else:
            num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            if len(num_cols) >= 2:
                df = df[num_cols[:2]].copy()
                df.columns = ["mass_kDa", "intensity"]
            else:
                raise ValueError("DataFrame must contain numeric mass and intensity columns.")
        if col_map:
            df = df.rename(columns=col_map)
    else:
        raise TypeError("source must be str or pd.DataFrame.")

    df["mass_kDa"] = pd.to_numeric(df["mass_kDa"], errors="coerce")
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df = df.dropna().sort_values("mass_kDa").reset_index(drop=True)
    return df


def plot_single_cdms_mass(
    dataset: Union[str, pd.DataFrame],
    label: str,
    second_label: Optional[str] = None,
    x_range: Tuple[float, float] = (0, 500),
    color: str = "#2C7FB8",
    normalize: bool = True,
    norm_mode: str = "percent",  # "percent" or "unit"
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False
):
    """
    Plot a single CD-MS mass distribution with optional normalization.
    norm_mode:
        - "percent": 0–100 with ticks at 0, 50, 100
        - "unit":    0–1 with ticks at 0, 0.5, 1
    """
    with style_context(font_name=font_name, font_size=font_size):
        fig, ax = plt.subplots(1, 1, figsize=figsize, dpi=dpi)

        df = _load_mass(dataset)
        x = df["mass_kDa"].values.astype(float)
        y_raw = df["intensity"].values.astype(float)

        if normalize:
            max_i = float(np.nanmax(y_raw)) if y_raw.size else 0.0
            if max_i > 0:
                if norm_mode == "percent":
                    y = (y_raw / max_i) * 100.0
                else:
                    y = y_raw / max_i
            else:
                y = np.zeros_like(y_raw)
        else:
            y = y_raw

        ax.plot(x, y, color=color, linewidth=1.8, label=label)
        ax.set_xlim(*x_range)

        if normalize:
            if norm_mode == "percent":
                ax.set_ylim(0, 100)
                ax.set_yticks([0, 50, 100])
                ax.set_yticklabels(["0", "%", "100"])
            else:
                ax.set_ylim(0, 1.0)
                ax.set_yticks([0.0, 0.5, 1.0])
                ax.set_yticklabels(["0", "0.5", "1.0"])
            ax.set_ylabel("Normalized Intensity")
        else:
            ymax = float(np.nanmax(y)) if y.size else 1.0
            ax.set_ylim(0, ymax * 1.05)
            ax.set_ylabel("Intensity")

        ax.set_xlabel("Mass (kDa)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False, loc="upper left")

        if second_label:
            place_second_label(ax, x, y, second_label, color, x_range, font_size)

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
    y_tick_labels: Optional[Sequence[str]] = None,  # not heavily used now
    normalize: bool = True,
    norm_mode: str = "percent",  # "percent" or "unit"
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False
):
    """
    Offset-stacked CD-MS mass plots.
    norm_mode:
        - "percent": each trace 0–100 over its offset
        - "unit":    each trace 0–1 over its offset
    """
    n = len(datasets)
    if len(labels) != n:
        raise ValueError("labels length must match datasets length.")
    if second_labels is None:
        second_labels = [None] * n
    if len(second_labels) != n:
        raise ValueError("second_labels length must match datasets length (or None).")

    colors = ensure_colors(colors, n, DEFAULT_CDMS_COLORS)

    width_in = figsize[0] if figsize and len(figsize) > 0 else figsize[0]
    target_h = compute_target_height(n_rows=n)

    with style_context(font_name=font_name, font_size=font_size):
        fig, axes = plt.subplots(n, 1, sharex=True, figsize=(width_in, target_h), dpi=dpi)
        if n == 1:
            axes = [axes]

        dfs = [_load_mass(src) for src in datasets]

        for i, (df, label, second, color) in enumerate(zip(dfs, labels, second_labels, colors)):
            ax = axes[i]
            offset = i * y_offset_step

            x = df["mass_kDa"].values.astype(float)
            y_raw = df["intensity"].values.astype(float)

            if normalize:
                max_i = float(np.nanmax(y_raw)) if y_raw.size else 0.0
                if max_i > 0:
                    if norm_mode == "percent":
                        y = (y_raw / max_i) * 100.0
                    else:
                        y = y_raw / max_i
                else:
                    y = np.zeros_like(y_raw)
            else:
                y = y_raw

            ax.plot(x, y + offset, color=color, linewidth=1.5, label=label)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.set_xlim(*x_range)

            if normalize:
                if norm_mode == "percent":
                    ax.set_ylim(offset, offset + 101.5)
                    ax.set_yticks([offset, offset + 50, offset + 100])
                    ax.set_yticklabels(["0", "%", "100"])
                else:
                    ax.set_ylim(offset, offset + 1.05)
                    ax.set_yticks([offset, offset + 0.5, offset + 1.0])
                    ax.set_yticklabels(["0", "0.5", "1.0"])
            else:
                ymax = float(np.nanmax(y)) if y.size else 1.0
                ax.set_ylim(offset, offset + ymax * 1.05)
                ax.set_yticks([offset, offset + ymax])
                ax.set_yticklabels(["", ""])

            leg1 = ax.legend([ax.lines[-1]], [label], loc="upper left", frameon=False, fontsize=font_size)
            ax.add_artist(leg1)

            if second:
                place_second_label(ax, x, y + offset, second, color, x_range, font_size)

            if i < n - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Mass (kDa)")

        fig.text(0.05, 0.5,
                 "Normalized Intensity" if normalize else "Intensity",
                 va="center", rotation="vertical", fontsize=font_size)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95, hspace=0.6)

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, axes