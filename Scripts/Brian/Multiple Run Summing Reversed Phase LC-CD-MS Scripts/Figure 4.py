"""Create the three-panel staggered TIC comparison used in Figure 4."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TIME_LIMITS = (0.0, 5.5)
STAGGER_STEP = 10.0
FIGURE_SIZE_INCHES = (3.25, 6.5)
OUTPUT_DPI = 300
MINIMUM_FONT_SIZE = 10
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin\MartyLab - General\Papers"
    r"\Muliple Run Summing Reversed Phase LC-CD-MS\Figures for Manuscript\Figure 4"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 4.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 4.pdf"

ENSEMBLE_COLOR = "#333333"
ONE_INJECTION_COLOR = "#0072B2"
TWENTY_ONE_INJECTION_COLOR = "#D55E00"

PANELS = [
    (
        "A",
        "Waters Intact mAb\nMass Check Standard",
        [
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\Ensemble Waters mAb Triplicate"
                    r"\Waters mAb 1 TIC.txt"
                ),
                "Ensemble MS",
                ENSEMBLE_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT"
                    r"\Waters mAb 21 Injections Set 1\merged"
                    r"\20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_"
                    r"0o1mgmLs_merged_1_unidecfiles\TIC_chrom.txt"
                ),
                "1 INJ CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT"
                    r"\Waters mAb 21 Injections Set 1\merged"
                    r"\20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_"
                    r"0o1mgmLs_merged_21_unidecfiles\TIC_decon.txt"
                ),
                "21 INJ CD-MS",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
    (
        "B",
        "Human IgG",
        [
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper"
                    r"\Ensemble Amgen Human IgG Triplicate\Human IgG 1 TIC.txt"
                ),
                "Ensemble MS",
                ENSEMBLE_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT"
                    r"\Amgen Human IgG 21 Injections Set 1\merged"
                    r"\20260331_BHT_agilentc8column_AmgenhumanmAb_PRO35304_CMDS_"
                    r"20msINJ_set1_0o1mgmLs_merged_1_unidecfiles\TIC_chrom.txt"
                ),
                "1 INJ CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT"
                    r"\Amgen Human IgG 21 Injections Set 1\merged"
                    r"\20260331_BHT_agilentc8column_AmgenhumanmAb_PRO35304_CMDS_"
                    r"20msINJ_set1_0o1mgmLs_merged_21_unidecfiles\TIC_decon.txt"
                ),
                "21 INJ CD-MS",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
    (
        "C",
        "Rabbit IgG",
        [
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper"
                    r"\Ensemble Amgen Rabbit IgG Triplicate\Rabbit IgG 1 TIC.txt"
                ),
                "Ensemble MS",
                ENSEMBLE_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT"
                    r"\Amgen Rabbit IgG 21 Injections Set 1\merged"
                    r"\20260401_BHT_agilentc8column_AmgenrabbitmAb_EPA13724_CMDS_"
                    r"20msINJ_set1_0o1mgmLs_merged_1_unidecfiles\TIC_chrom.txt"
                ),
                "1 INJ CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT"
                    r"\Amgen Rabbit IgG 21 Injections Set 1\merged"
                    r"\20260401_BHT_agilentc8column_AmgenrabbitmAb_EPA13724_CMDS_"
                    r"20msINJ_set1_0o1mgmLs_merged_21_unidecfiles\TIC_decon.txt"
                ),
                "21 INJ CD-MS",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
]


def load_chromatogram(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a two-column chromatogram and retain the common time window."""
    if not path.is_file():
        raise FileNotFoundError(f"Chromatogram not found: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in: {path}")

    time = data[:, 0]
    intensity = data[:, 1]
    valid = (
        np.isfinite(time)
        & np.isfinite(intensity)
        & (time >= TIME_LIMITS[0])
        & (time <= TIME_LIMITS[1])
    )
    if not np.any(valid):
        raise ValueError(f"No finite values in the common time window: {path}")

    return time[valid], np.clip(intensity[valid], 0.0, None)


def normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Min-max normalize a chromatogram to exactly 0--100%."""
    minimum = float(np.min(intensity))
    maximum = float(np.max(intensity))
    span = maximum - minimum
    if span <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * (intensity - minimum) / span


def style_axis(axis: plt.Axes, panel_letter: str, sample_label: str) -> None:
    """Apply shared publication styling to one Figure 4 panel."""
    axis.set_xlim(*TIME_LIMITS)
    axis.set_ylim(0.0, 180.0)
    axis.set_ylabel(
        sample_label,
        fontsize=13,
        rotation=90,
        labelpad=22,
        ha="center",
        va="center",
    )
    # Center the sample name on the visible 0--100 portion of the y-axis.
    axis.yaxis.set_label_coords(-0.27, 50.0 / 180.0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(True)
    axis.spines["left"].set_linewidth(1.5)
    axis.spines["left"].set_bounds(0.0, 100.0)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.set_yticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
    axis.tick_params(
        direction="out",
        width=1.5,
        length=6,
        labelsize=MINIMUM_FONT_SIZE,
    )
    axis.text(
        0.015,
        0.80,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def make_figure() -> plt.Figure:
    """Build and return the complete three-panel Figure 4."""
    figure, axes = plt.subplots(
        nrows=3,
        ncols=1,
        sharex=True,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    offsets = STAGGER_STEP * np.arange(2, -1, -1)
    for axis, (panel_letter, sample_label, datasets) in zip(axes, PANELS):
        for (path, label, color), offset in zip(datasets, offsets):
            time, intensity = load_chromatogram(path)
            axis.plot(
                time,
                normalize_to_percent(intensity) + offset,
                color=color,
                linewidth=2.0,
                label=label,
                zorder=2,
            )

        style_axis(axis, panel_letter, sample_label)
        if panel_letter == "A":
            axis.legend(
                loc="upper right",
                bbox_to_anchor=(0.98, 0.80),
                frameon=False,
                fontsize=MINIMUM_FONT_SIZE,
                handlelength=2.0,
                borderaxespad=0.0,
            )

    axes[-1].set_xlabel("Time (min)")
    axes[-1].set_xticks(np.arange(0.0, 5.1, 1.0))

    return figure


def main() -> None:
    """Render Figure 4 as both PNG and PDF."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.labelsize": 13,
            "legend.fontsize": MINIMUM_FONT_SIZE,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )

    figure = make_figure()
    # Preserve the requested physical dimensions in both output formats.
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
