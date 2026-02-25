# mwd_chromatogram_utils.py
from __future__ import annotations
from typing import Sequence, Optional, Tuple
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# Consistent sizing
DEFAULT_WIDTH_IN = 10.0
ROW_HEIGHT_IN     = 1.35
MIN_HEIGHT_IN     = 3.8


def load_chrom_file(path: str, signal_keywords: Sequence[str]) -> Optional[pd.DataFrame]:
    try:
        df = pd.read_excel(path, header=None, engine="openpyxl")
    except Exception as e:
        print(f"Could not read {path}: {e}")
        return None

    header_row = None
    for idx in range(min(10, len(df))):
        row = df.iloc[idx].astype(str).str.lower().str.strip().values
        if any("ml" in s for s in row) and any(any(key in s for key in signal_keywords) for s in row):
            header_row = idx
            break

    if header_row is not None:
        df = pd.read_excel(path, header=header_row, engine="openpyxl")
    else:
        df = pd.read_excel(path, header=2, engine="openpyxl")

    cols_lower = [str(c).lower().strip() for c in df.columns]
    ml_col, sig_col = None, None
    for j, cname in enumerate(cols_lower):
        if "ml" in cname:
            ml_col = df.columns[j]
        for key in signal_keywords:
            if key in cname:
                sig_col = df.columns[j]
    if ml_col is None or sig_col is None:
        ml_col, sig_col = df.columns[0], df.columns[1]

    df[ml_col] = pd.to_numeric(df[ml_col], errors="coerce")
    df[sig_col] = pd.to_numeric(df[sig_col], errors="coerce")
    clean = df[[ml_col, sig_col]].dropna()
    clean.columns = ["ml", "signal"]
    clean = clean.sort_values("ml").reset_index(drop=True)
    return clean


def plot_mwd_chromatograms(
    mau_file: str,
    mv_file: str,
    mau_label: str = "mAU",
    mv_label: str = "mV",
    colors: Optional[Sequence[str]] = None,
    xlim: Tuple[float, float] = (0.0, 26.0),
    normalize: bool = True,
    alpha_mv: float = 0.9,
    font_name: str = "Arial",
    font_size: int = 12,
    figsize: Tuple[float, float] = (10, 4.2),
    dpi: int = 300,
    tick_length: float = 5.0,
    tick_width: float = 1.0,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Overlay MWD chromatograms (mAU and mV) using twin y-axes.

    Styling and consistency:
      - Unified fonts and tick styling.
      - Percent y-ticks (0, %, 100) when normalize=True, independently per axis.
      - Frameless, combined legend at left center.
      - Clean spines; width fixed, height comfortable.

    Returns:
      fig, (ax1, ax2) where ax1 is mAU (left y-axis) and ax2 is mV (right y-axis).
    """
    if colors is None or len(colors) < 2:
        colors = ["#0F0C08", "#3CAEA3"]

    mpl.rcParams.update({
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

    # Respect provided width; ensure comfortable height
    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = max(MIN_HEIGHT_IN, ROW_HEIGHT_IN + 0.8)

    fig, ax1 = plt.subplots(1, 1, figsize=(width_in, target_h), dpi=dpi)

    # mAU trace (left y-axis)
    df_mau = load_chrom_file(mau_file, ["mau"])
    if df_mau is not None:
        if normalize:
            max_val = float(df_mau["signal"].max())
            if max_val > 0:
                df_mau["signal"] = (df_mau["signal"] / max_val) * 100.0
        ax1.plot(df_mau["ml"], df_mau["signal"], color=colors[0], label=mau_label, linewidth=1.8)

    # mV trace (right y-axis)
    ax2 = ax1.twinx()
    df_mv = load_chrom_file(mv_file, ["mv"])
    if df_mv is not None:
        if normalize:
            max_val = float(df_mv["signal"].max())
            if max_val > 0:
                df_mv["signal"] = (df_mv["signal"] / max_val) * 100.0
        ax2.plot(df_mv["ml"], df_mv["signal"], color=colors[1], label=mv_label, linewidth=1.8, alpha=alpha_mv)

    # Limits
    ax1.set_xlim(*xlim)
    ax1.margins(y=0)
    ax2.margins(y=0)

    # Spines cleanup
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)  # avoid double left spine

    # Y ticks
    if normalize:
        ax1.set_ylim(0, 100)
        ax2.set_ylim(0, 100)
        ax1.set_yticks([0, 50, 100])
        ax2.set_yticks([0, 50, 100])
        ax1.set_yticklabels(["0", "%", "100"])
        ax2.set_yticklabels(["0", "%", "100"])
        ax1.set_ylabel("Relative Absorbance (%)")
        ax2.set_ylabel("Relative Voltage (%)")
    else:
        # endpoints based on data
        if df_mau is not None and len(df_mau):
            ax1.set_ylim(0, max(1.0, float(df_mau["signal"].max())))
        if df_mv is not None and len(df_mv):
            ax2.set_ylim(0, max(1.0, float(df_mv["signal"].max())))
        ax1.set_ylabel("Absorbance (mAU)")
        ax2.set_ylabel("Voltage (mV)")

    # X label
    ax1.set_xlabel("Volume (mL)")

    # Tick styling
    ax1.tick_params(axis="x", length=tick_length, width=tick_width)
    ax1.tick_params(axis="y", length=tick_length, width=tick_width)
    ax2.tick_params(axis="y", length=tick_length, width=tick_width)

    # Combined legend (frameless), left-center
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="center left", bbox_to_anchor=(0.02, 0.5))

    # Layout
    plt.tight_layout(rect=(0.08, 0.08, 0.98, 0.96))

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()

    return fig, (ax1, ax2)