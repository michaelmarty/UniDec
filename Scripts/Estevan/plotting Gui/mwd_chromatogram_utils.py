# mwd_chromatogram_utils.py
from __future__ import annotations
from typing import Sequence, Optional, Tuple
import pandas as pd
import matplotlib.pyplot as plt

def load_chrom_file(path: str, signal_keywords: Sequence[str]) -> Optional[pd.DataFrame]:
    """
    Load chromatogram from Excel and detect 'ml' + signal column using header sniffing.
    Returns DataFrame ['ml', 'signal'] sorted by ml, or None on failure.
    """
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
    alpha_mv: float = 0.8,
    font_name: str = "Arial",
    font_size: int = 12,
    figsize: Tuple[float, float] = (6, 4),
    dpi: int = 300,
    save: Optional[str] = None,
    show: bool = False,
):
    """
    Dual-axis chromatograms: mAU (left) + mV (right), normalized and combined legend.
    Returns (fig, (ax1, ax2)).
    """
    if colors is None:
        colors = ["#0F0C08", "#3CAEA3"]

    plt.rcParams.update({
        "font.family": font_name,
        "font.size": font_size,
        "axes.labelsize": font_size,
        "xtick.labelsize": max(font_size - 2, 6),
        "ytick.labelsize": max(font_size - 2, 6),
    })

    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(111)

    df_mau = load_chrom_file(mau_file, ["mau"])
    if df_mau is not None:
        if normalize:
            max_val = df_mau["signal"].max()
            if max_val > 0:
                df_mau["signal"] = df_mau["signal"] / max_val
        ax1.plot(df_mau["ml"], df_mau["signal"], color=colors[0], label=mau_label, linewidth=1.3)

    ax2 = ax1.twinx()
    df_mv = load_chrom_file(mv_file, ["mv"])
    if df_mv is not None:
        if normalize:
            max_val = df_mv["signal"].max()
            if max_val > 0:
                df_mv["signal"] = df_mv["signal"] / max_val
        ax2.plot(df_mv["ml"], df_mv["signal"], color=colors[1], label=mv_label, linewidth=1.3, alpha=alpha_mv)

    ax1.set_xlim(*xlim)
    ax1.set_ylim(0, 1 if normalize else max(1, df_mau["signal"].max() if df_mau is not None else 1))
    ax2.set_ylim(0, 1 if normalize else max(1, df_mv["signal"].max() if df_mv is not None else 1))
    ax1.set_xlabel("Retention Volume (mL)")
    ax1.set_ylabel("Normalized Absorbance (mAU)" if normalize else "Absorbance (mAU)")
    ax2.set_ylabel("Normalized Voltage (mV)" if normalize else "Voltage (mV)")

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_position(("data", 0))
    ax1.spines["bottom"].set_position(("data", 0))

    l1, lab1 = ax1.get_legend_handles_labels()
    l2, lab2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lab1 + lab2, frameon=False, loc="upper right")

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return fig, (ax1, ax2)