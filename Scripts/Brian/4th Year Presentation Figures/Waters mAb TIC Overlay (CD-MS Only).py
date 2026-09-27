from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATA_ROOT = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\Waters mAb 21 Injections Set 1"
)
MERGED_ROOT = DATA_ROOT / "merged"
TWENTY_ONE_INJECTION_DIR = MERGED_ROOT / (
    "20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_0o1mgmLs_merged_21_"
    "unidecfiles"
)
ONE_INJECTION_DIR = MERGED_ROOT / (
    "20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_0o1mgmLs_merged_1_"
    "unidecfiles"
)
TIC_21 = TWENTY_ONE_INJECTION_DIR / "TIC_decon.txt"
TIC_1 = ONE_INJECTION_DIR / "TIC_chrom.txt"

TIME_LIMITS = (0.0, 5.5)
FIGURE_SIZE_INCHES = (8, 6)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Presentations\Brian\4th Year Presentation\Figures"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Waters mAb TIC Overlay (CD-MS Only).png"


def load_chromatogram(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load and validate a two-column time/intensity chromatogram."""
    if not path.is_file():
        raise FileNotFoundError(f"Chromatogram not found: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in: {path}")

    time = data[:, 0]
    intensity = data[:, 1]
    finite = np.isfinite(time) & np.isfinite(intensity)
    if not np.any(finite):
        raise ValueError(f"No finite time/intensity pairs in: {path}")

    # Deconvolution can leave numerical noise around -1e-14. Ion counts are
    # physically non-negative, so remove only those negative values here.
    return time[finite], np.clip(intensity[finite], 0.0, None)


def style_axis(axis: plt.Axes) -> None:
    """Apply presentation styling without a panel letter."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(2.0)
    axis.spines["bottom"].set_linewidth(2.0)
    axis.tick_params(direction="out", width=2.0, length=6)
    axis.set_xlim(*TIME_LIMITS)
    axis.set_xlabel("Time (min)")


def make_figure() -> plt.Figure:
    """Compare raw one- and twenty-one-injection CD-MS ion counts."""
    comparison = [
        (TIC_1, "1 INJ CD-MS", "#0072B2"),
        (TIC_21, "21 INJ CD-MS", "#D55E00"),
    ]
    figure, axis = plt.subplots(
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )
    maximum_ion_count = 0.0
    for path, label, color in comparison:
        time, intensity = load_chromatogram(path)
        maximum_ion_count = max(maximum_ion_count, float(np.max(intensity)))
        axis.plot(
            time,
            intensity,
            color=color,
            linewidth=2.5,
            label=label,
        )

    style_axis(axis)
    axis.set_ylabel("Ion count")
    axis.set_ylim(0.0, maximum_ion_count * 1.40)
    axis.legend(
        loc="upper right",
        frameon=False,
        handlelength=1.5,
    )

    return figure


if __name__ == "__main__":
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 24,
            "axes.labelsize": 28,
            "legend.fontsize": 20,
            "savefig.transparent": True,
        }
    )
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    figure = make_figure()
    # Preserve the exact 8 by 6 inch canvas without bbox_inches="tight".
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
