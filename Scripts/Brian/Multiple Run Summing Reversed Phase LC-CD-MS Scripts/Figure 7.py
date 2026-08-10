"""Create the high-throughput three-panel staggered TIC comparison."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


TIME_LIMITS = (0.0, 3.0)
STAGGER_STEP = 10.0
FIGURE_SIZE_INCHES = (3.25, 6.5)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 7.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 7.pdf"

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
                    r"Z:\Group Share\BHT\Temp\Paper\High Throughput"
                    r"\Ensemble Waters mAb Triplicate\Waters mAb 1 TIC.txt"
                ),
                "Ensemble MS",
                ENSEMBLE_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
                    r"\Waters mAb 21 Injections Set 1\merged"
                    r"\20260501_BHT_haloc4column_watersmAb_set1_CDMS_0o1mgmL_"
                    r"merged_1_unidecfiles\TIC_chrom.txt"
                ),
                "1-injection CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
                    r"\Waters mAb 21 Injections Set 1\merged"
                    r"\20260501_BHT_haloc4column_watersmAb_set1_CDMS_0o1mgmL_"
                    r"merged_21_unidecfiles\TIC_decon.txt"
                ),
                "21-injection CD-MS",
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
                    r"Z:\Group Share\BHT\Temp\Paper\High Throughput"
                    r"\Ensemble Amgen Human IgG Triplicate\Human IgG 1 TIC.txt"
                ),
                "Ensemble MS",
                ENSEMBLE_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
                    r"\Amgen Human IgG 21 Injections Set 1\merged"
                    r"\20260501_BHT_haloc4column_amgenHuman_PRO35304_set1_CDMS_"
                    r"0o1mgmL_merged_1_unidecfiles\TIC_chrom.txt"
                ),
                "1-injection CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
                    r"\Amgen Human IgG 21 Injections Set 1\merged"
                    r"\20260501_BHT_haloc4column_amgenHuman_PRO35304_set1_CDMS_"
                    r"0o1mgmL_merged_21_unidecfiles\TIC_decon.txt"
                ),
                "21-injection CD-MS",
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
                    r"Z:\Group Share\BHT\Temp\Paper\High Throughput"
                    r"\Ensemble Amgen Rabbit IgG Triplicate\Rabbit IgG 1 TIC.txt"
                ),
                "Ensemble MS",
                ENSEMBLE_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
                    r"\Amgen Rabbit IgG 21 Injections Set 1\merged"
                    r"\20260605_BHT_haloc4column_amgenRabbitmAb_EPA13724_set1_"
                    r"CDMS_0o01mgmL_merged_1_unidecfiles\TIC_chrom.txt"
                ),
                "1-injection CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                Path(
                    r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
                    r"\Amgen Rabbit IgG 21 Injections Set 1\merged"
                    r"\20260605_BHT_haloc4column_amgenRabbitmAb_EPA13724_set1_"
                    r"CDMS_0o01mgmL_merged_21_unidecfiles\TIC_decon.txt"
                ),
                "21-injection CD-MS",
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
    """Apply the same publication styling used in Figure 4."""
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
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.set_yticks([])
    axis.tick_params(axis="x", direction="out", width=1.5, length=6)
    axis.tick_params(axis="y", left=False, labelleft=False)
    axis.text(
        0.015,
        0.97,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def make_figure() -> plt.Figure:
    """Build and return the complete three-panel Figure 7."""
    figure, axes = plt.subplots(
        nrows=3,
        ncols=1,
        sharex=True,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    offsets = STAGGER_STEP * np.arange(2, -1, -1)
    for axis, (panel_letter, sample_label, datasets) in zip(axes, PANELS):
        retention_time = None
        for (path, label, color), offset in zip(datasets, offsets):
            time, intensity = load_chromatogram(path)
            if (
                panel_letter == "C" and label == "Ensemble MS"
            ) or (
                panel_letter != "C" and label == "21-injection CD-MS"
            ):
                retention_time = float(time[np.argmax(intensity)])
            axis.plot(
                time,
                normalize_to_percent(intensity) + offset,
                color=color,
                linewidth=2.0,
                label=label,
                zorder=2,
            )

        style_axis(axis, panel_letter, sample_label)
        if retention_time is not None:
            axis.axvline(
                retention_time,
                ymin=0.0,
                ymax=0.72,
                color="#000000",
                linestyle=(0, (4, 3)),
                linewidth=1.5,
                zorder=3,
            )
            axis.text(
                retention_time - 0.05,
                0.73,
                f"{retention_time:.2f} min",
                transform=axis.get_xaxis_transform(),
                color="#000000",
                fontsize=7.0,
                ha="right",
                va="bottom",
                zorder=4,
            )
        if panel_letter == "A":
            axis.legend(
                loc="upper right",
                frameon=False,
                fontsize=9.5,
                handlelength=2.0,
                borderaxespad=0.2,
            )

    axes[-1].set_xlabel("Time (min)")
    axes[-1].set_xticks(np.arange(0.0, 3.1, 0.5))

    return figure


def main() -> None:
    """Render Figure 7 as a transparent PNG and PDF."""
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
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
