from __future__ import annotations
from typing import Sequence, Optional, Tuple
import pandas as pd
import matplotlib.pyplot as plt

from plot_helpers import style_context, ensure_colors, compute_target_height, DEFAULT_MWD_COLORS

DEFAULT_WIDTH_IN = 10.0
ROW_HEIGHT_IN     = 1.35
MIN_HEIGHT_IN     = 3.8


def load_chrom_file(path: str, signal_keywords: Sequence[str]) -> Optional[pd.DataFrame]:
    try:
        raw = pd.read_excel(path, header=None, engine="openpyxl")
    except Exception:
        raw = None

    header_row = None
    df = None

    if raw is not None:
        for idx in range(min(10, len(raw))):
            row = raw.iloc[idx].astype(str).str.lower().str.strip().values
            if any("ml" in s for s in row) and any(any(key in s for key in signal_keywords) for s in row):
                header_row = idx
                break
        try:
            if header_row is not None:
                df = pd.read_excel(path, header=header_row, engine="openpyxl")
            else:
                df = pd.read_excel(path, header=2, engine="openpyxl")
        except Exception:
            df = None

    if df is None:
        try:
            df_try = pd.read_csv(path, header=None, engine="python")
            num_cols = [c for c in df_try.columns if pd.api.types.is_numeric_dtype(df_try[c])]
            if len(num_cols) >= 2:
                df = df_try[num_cols[:2]].copy()
            else:
                df = df_try.iloc[:, :2].copy()
        except Exception:
            return None

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
    colors = ensure_colors(colors, 2, DEFAULT_MWD_COLORS)

    width_in = figsize[0] if figsize and len(figsize) > 0 else DEFAULT_WIDTH_IN
    target_h = compute_target_height(n_rows=1)

    with style_context(font_name=font_name, font_size=font_size, tick_length=tick_length, tick_width=tick_width):
        fig, ax1 = plt.subplots(1, 1, figsize=(width_in, target_h), dpi=dpi)

        df_mau = load_chrom_file(mau_file, ["mau"])
        if df_mau is not None:
            df_mau_local = df_mau.copy()
            if normalize:
                max_val = float(df_mau_local["signal"].max()) if not df_mau_local["signal"].empty else 0.0
                if max_val > 0:
                    df_mau_local["signal"] = (df_mau_local["signal"] / max_val) * 100.0
            ax1.plot(df_mau_local["ml"], df_mau_local["signal"], color=colors[0], label=mau_label, linewidth=1.8)

        ax2 = ax1.twinx()
        df_mv = load_chrom_file(mv_file, ["mv"])
        if df_mv is not None:
            df_mv_local = df_mv.copy()
            if normalize:
                max_val = float(df_mv_local["signal"].max()) if not df_mv_local["signal"].empty else 0.0
                if max_val > 0:
                    df_mv_local["signal"] = (df_mv_local["signal"] / max_val) * 100.0
            ax2.plot(df_mv_local["ml"], df_mv_local["signal"], color=colors[1], label=mv_label, linewidth=1.8, alpha=alpha_mv)

        ax1.set_xlim(*xlim)
        ax1.margins(y=0)
        ax2.margins(y=0)

        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax2.spines["top"].set_visible(False)
        ax2.spines["left"].set_visible(False)

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
            if df_mau is not None and len(df_mau):
                ax1.set_ylim(0, max(1.0, float(df_mau["signal"].max())))
            if df_mv is not None and len(df_mv):
                ax2.set_ylim(0, max(1.0, float(df_mv["signal"].max())))
            ax1.set_ylabel("Absorbance (mAU)")
            ax2.set_ylabel("Voltage (mV)")

        ax1.set_xlabel("Volume (mL)")

        ax1.tick_params(axis="x", length=tick_length, width=tick_width)
        ax1.tick_params(axis="y", length=tick_length, width=tick_width)
        ax2.tick_params(axis="y", length=tick_length, width=tick_width)

        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        if h1 or h2:
            ax1.legend(h1 + h2, l1 + l2, loc="center left", bbox_to_anchor=(0.02, 0.5))

        plt.tight_layout(rect=(0.08, 0.08, 0.98, 0.96))

        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()

    return fig, (ax1, ax2)