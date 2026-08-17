"""Create the three-panel chromatogram comparison used in Figure 2."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATA_ROOT = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\Protein Mix 21 Injections Set 1"
)
MERGED_ROOT = DATA_ROOT / "merged"
TWENTY_ONE_INJECTION_DIR = MERGED_ROOT / (
    "20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21_"
    "unidecfiles"
)
ONE_INJECTION_DIR = MERGED_ROOT / (
    "20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_1_"
    "unidecfiles"
)
ENSEMBLE_CHROMATOGRAM = Path(
    r"Z:\Group Share\BHT\Temp\Paper\Ensemble Protein Mix Triplicate"
) / "UniDec_Figures_and_Files" / "Protein mix Rep 1 chrom.txt"

TIC_21 = TWENTY_ONE_INJECTION_DIR / "TIC_decon.txt"
TIC_1 = ONE_INJECTION_DIR / "TIC_chrom.txt"

# The one-injection TIC is zero-padded after the chromatographic run. All of the
# actual signals in the supplied files lie within this common time window.
TIME_LIMITS = (0.0, 5.5)

# Panel B traces are independently normalized to 0--100% first, then shifted
# slightly so that coincident features remain visible.
STAGGER_STEP = 10.0

OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 2.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 2.pdf"
FIGURE_SIZE_INCHES = (3.25, 6.5)
OUTPUT_DPI = 300


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


def normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Min-max normalize a signal to exactly 0--100%."""
    minimum = float(np.min(intensity))
    maximum = float(np.max(intensity))
    span = maximum - minimum
    if span <= 0:
        return np.zeros_like(intensity)
    return 100.0 * (intensity - minimum) / span


def style_axis(
    axis: plt.Axes, panel_letter: str, panel_letter_y: float = 0.97
) -> None:
    """Apply shared publication styling and add a panel letter."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.5)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.tick_params(direction="out", width=1.5, length=6)
    axis.set_xlim(*TIME_LIMITS)
    axis.text(
        0.015,
        panel_letter_y,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def hide_y_axis(axis: plt.Axes) -> None:
    """Hide the y-axis while retaining the data scaling."""
    axis.set_ylabel("")
    axis.set_yticks([])
    axis.spines["left"].set_visible(False)
    axis.tick_params(axis="y", left=False, labelleft=False)


def style_percent_axis(axis: plt.Axes, upper_limit: float = 104.0) -> None:
    """Set the normalized range while keeping the y-axis hidden."""
    axis.set_ylim(0.0, upper_limit)
    hide_y_axis(axis)


def make_figure() -> plt.Figure:
    """Build and return the complete three-row, one-column figure."""
    panel_a = [
        (TIC_21, "Total ion chromatogram", "#000000"),
        (
            TWENTY_ONE_INJECTION_DIR / "Swoop_mz__1130_z__11_decon.txt",
            "Cytochrome C EIC",
            "#0072B2",
        ),
        (
            TWENTY_ONE_INJECTION_DIR / "Swoop_mz__1210_z__14_decon.txt",
            "Myoglobin EIC",
            "#009E73",
        ),
        (
            TWENTY_ONE_INJECTION_DIR / "Swoop_mz__1337_z__22_decon.txt",
            "Carbonic Anhydrase EIC",
            "#CC79A7",
        ),
    ]
    comparison = [
        (ENSEMBLE_CHROMATOGRAM, "Ensemble MS", "#333333"),
        (TIC_1, "1-injection CD-MS", "#0072B2"),
        (TIC_21, "21-injection CD-MS", "#D55E00"),
    ]

    figure, axes = plt.subplots(
        nrows=3,
        ncols=1,
        sharex=True,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    # A: raw TIC and protein-specific EIC overlays (no normalization).
    for path, label, color in panel_a:
        time, intensity = load_chromatogram(path)
        if path != TIC_21:
            retention_time = float(time[np.argmax(intensity)])
            axes[0].axvline(
                retention_time,
                color=color,
                linestyle=(0, (4, 3)),
                linewidth=1.5,
                alpha=1.0,
                zorder=3,
            )
            protein_name = label.removesuffix(" EIC")
            if protein_name in {"Cytochrome C", "Myoglobin"}:
                label_x = retention_time - (
                    0.01 if protein_name == "Myoglobin" else 0.08
                )
                label_y = 0.79 if protein_name == "Myoglobin" else 0.96
                horizontal_alignment = "right"
            else:
                label_x = retention_time + (
                    0.12 if protein_name == "Carbonic Anhydrase" else 0.08
                )
                label_y = 0.96
                horizontal_alignment = "left"
            axes[0].text(
                label_x,
                label_y,
                f"{protein_name}\n{retention_time:.2f} min",
                transform=axes[0].get_xaxis_transform(),
                color=color,
                fontsize=7.0,
                fontweight="normal",
                ha=horizontal_alignment,
                va="top",
                zorder=4,
            )
        axes[0].plot(
            time,
            intensity,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=2,
        )
    style_axis(axes[0], "A", panel_letter_y=0.78)
    axes[0].set_ylim(bottom=0.0)
    hide_y_axis(axes[0])

    # B: normalized comparison with a small vertical offset between traces.
    offsets = STAGGER_STEP * np.arange(len(comparison) - 1, -1, -1)
    for (path, label, color), offset in zip(comparison, offsets):
        time, intensity = load_chromatogram(path)
        axes[1].plot(
            time,
            normalize_to_percent(intensity) + offset,
            color=color,
            linewidth=2.0,
            label=label,
        )
    style_axis(axes[1], "B")
    # Leave a clear band above the staggered traces for the three-row legend.
    style_percent_axis(axes[1], upper_limit=180.0)
    axes[1].legend(
        loc="upper left",
        bbox_to_anchor=(0.12, 1.0),
        ncols=1,
        frameon=False,
        fontsize=9.5,
        handlelength=2.0,
        borderaxespad=0.0,
    )

    # C: raw (not normalized) one- and twenty-one-injection CD-MS TICs.
    maximum_raw_ion_count = 0.0
    for path, label, color in comparison[1:]:
        time, intensity = load_chromatogram(path)
        maximum_raw_ion_count = max(maximum_raw_ion_count, float(np.max(intensity)))
        axes[2].plot(
            time,
            intensity,
            color=color,
            linewidth=2.0,
            label=label,
        )
    style_axis(axes[2], "C")
    axes[2].set_ylabel("Ion count")
    axes[2].set_xlabel("Time (min)")
    # Reserve headroom for the legend instead of placing it over a raw TIC peak.
    axes[2].set_ylim(0.0, maximum_raw_ion_count * 1.40)
    axes[2].legend(
        loc="upper left",
        bbox_to_anchor=(0.12, 1.0),
        frameon=False,
        fontsize=9.5,
        borderaxespad=0.0,
    )

    return figure


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.labelsize": 13,
            "legend.fontsize": 9.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )

    figure = make_figure()
    # Do not use bbox_inches="tight" here: it changes the requested physical
    # dimensions. The PNG is exactly 3.25 inches wide at 300 DPI (975 pixels).
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
