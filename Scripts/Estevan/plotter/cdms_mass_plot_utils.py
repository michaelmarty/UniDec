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
    df = pd.DataFrame(arr[:, :2], columns=["mass_kDa", "intensity"]).dropna()

    df["mass_kDa"] = pd.to_numeric(df["mass_kDa"], errors="coerce")
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df = df.dropna()

    max_mass = float(df["mass_kDa"].max()) if not df.empty else 0.0
    if max_mass > 5000:
        df["mass_kDa"] = df["mass_kDa"] / 1000.0

    return df


def _load_mass(source: Union[str, pd.DataFrame]) -> pd.DataFrame:
    if isinstance(source, str):
        try:
            df = read_cdms_mass_file(source)
        except Exception:
            df = pd.read_csv(source, sep=r"\s+", header=None)
            df = df.iloc[:, :2].copy()
            df.columns = ["mass_kDa", "intensity"]
    elif isinstance(source, pd.DataFrame):
        df = source.copy()
        if {"mass_kDa", "intensity"}.issubset(df.columns):
            df = df[["mass_kDa", "intensity"]]
        elif {"Mass_kDa", "Intensity"}.issubset(df.columns):
            df = df.rename(columns={"Mass_kDa": "mass_kDa", "Intensity": "intensity"})
            df = df[["mass_kDa", "intensity"]]
        elif {"Mass_Da", "Intensity"}.issubset(df.columns):
            df["mass_kDa"] = pd.to_numeric(df["Mass_Da"], errors="coerce") / 1000.0
            df["intensity"] = pd.to_numeric(df["Intensity"], errors="coerce")
            df = df[["mass_kDa", "intensity"]]
        else:
            num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            if len(num_cols) >= 2:
                df = df[num_cols[:2]].copy()
                df.columns = ["mass_kDa", "intensity"]
            else:
                raise ValueError("DataFrame must contain numeric mass and intensity columns.")
    else:
        raise TypeError("source must be str or pd.DataFrame.")

    df["mass_kDa"] = pd.to_numeric(df["mass_kDa"], errors="coerce")
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df = df.dropna()

    max_mass = float(df["mass_kDa"].max()) if not df.empty else 0.0
    if max_mass > 5000:
        df["mass_kDa"] = df["mass_kDa"] / 1000.0

    df = df.sort_values("mass_kDa").reset_index(drop=True)
    return df


def plot_single_cdms_mass(
    dataset: Union[str, pd.DataFrame],
    label: str,
    second_label: Optional[str] = None,
    x_range: Tuple[float, float] = (0, 500),
    color: str = "#2C7FB8",
    normalize: bool = True,
    norm_mode: str = "percent",
    y_max: Optional[float] = None,
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False
):
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


def plot_cdms_mass_stack(
    datasets: Sequence[Union[str, pd.DataFrame]],
    labels: Sequence[str],
    second_labels: Optional[Sequence[Optional[str]]] = None,
    x_range: Tuple[float, float] = (0, 500),
    colors: Optional[Sequence[str]] = None,
    y_max: Optional[float] = None,
    normalize: bool = True,
    norm_mode: str = "percent",
    font_name: str = "Arial",
    font_size: int = 10,
    figsize: Tuple[float, float] = (9, 5),
    dpi: int = 300,
    hspace: float = 0.6,
    save: Optional[str] = None,
    show: bool = False
):
    """
    Stacked CD-MS mass plots: one panel per dataset, each with its own baseline.
    No vertical offset between traces; 'hspace' controls the gaps between panels.
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

    with style_context(font_name=font_name, font_size=font_size):
        fig, axes = plt.subplots(n, 1, sharex=True, figsize=(width_in, target_h), dpi=dpi)
        if n == 1:
            axes = [axes]

        dfs = [_load_mass(src) for src in datasets]

        for i, (df, label, second, color) in enumerate(zip(dfs, labels, second_labels, colors)):
            ax = axes[i]

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

            ax.plot(x, y, color=color, linewidth=1.5, label=label)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
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
                ymax = float(np.nanmax(y)) if y.size else 1.0
                ax.set_ylim(0.0, ymax * 1.05)
                ax.set_yticks([0.0, ymax])

            leg1 = ax.legend([ax.lines[-1]], [label], loc="upper left", frameon=False, fontsize=font_size)
            ax.add_artist(leg1)

            if second:
                place_second_label(ax, x, y, second, color, x_range, font_size)

            if i < n - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Mass (kDa)")

        fig.text(
            0.05,
            0.5,
            "Normalized Intensity" if normalize else "Intensity",
            va="center",
            rotation="vertical",
            fontsize=font_size,
        )
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95, hspace=hspace)

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, axes