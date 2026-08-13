"""Compare Waters mAb TICs from the 8- and 3-minute separations."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "New Figure 8 vs 3 Min TIC Waters mAb.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "New Figure 8 vs 3 Min TIC Waters mAb.pdf"
FIGURE_SIZE_INCHES = (6.5, 3.25)

Y_LIMITS = (0.0, 115.0)

STANDARD_METHOD_COLOR = "#0072B2"
HIGH_THROUGHPUT_COLOR = "#D55E00"

TIC_DATASETS = [
    (
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Waters mAb 21 Injections Set 1\merged"
            r"\20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_"
            r"0o1mgmLs_merged_21_unidecfiles\TIC_decon.txt"
        ),
        (0.0, 5.5),
        "8-min method (Agilent C8)",
        STANDARD_METHOD_COLOR,
    ),
    (
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
            r"\Waters mAb 21 Injections Set 1\merged"
            r"\20260501_BHT_haloc4column_watersmAb_set1_CDMS_0o1mgmL_"
            r"merged_21_unidecfiles\TIC_decon.txt"
        ),
        (0.0, 3.0),
        "3-min method (HALO C4)",
        HIGH_THROUGHPUT_COLOR,
    ),
]


def load_chromatogram(
    path: Path,
    time_limits: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Load and validate one two-column TIC file."""
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
        & (time >= time_limits[0])
        & (time <= time_limits[1])
    )
    if not np.any(valid):
        raise ValueError(f"No finite values in the comparison window: {path}")

    return time[valid], np.clip(intensity[valid], 0.0, None)


def normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Min-max normalize a TIC to exactly 0--100%."""
    minimum = float(np.min(intensity))
    maximum = float(np.max(intensity))
    span = maximum - minimum
    if span <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * (intensity - minimum) / span


def style_axis(
    axis: plt.Axes,
    panel_letter: str,
    time_limits: tuple[float, float],
    show_y_axis: bool,
) -> None:
    """Apply the publication styling used by the related TIC figures."""
    display_limits = (
        (time_limits[0], time_limits[1] + 0.05)
        if time_limits[1] <= 3.0
        else time_limits
    )
    axis.set_xlim(*display_limits)
    axis.set_ylim(*Y_LIMITS)
    if time_limits[1] <= 3.0:
        axis.set_xticks(np.arange(0.0, 3.1, 0.5))
    else:
        axis.set_xticks(np.arange(0.0, 5.1, 1.0))
    axis.set_xlabel("Time (min)")

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.tick_params(direction="out", width=1.5, length=6)

    if show_y_axis:
        axis.spines["left"].set_linewidth(1.5)
        axis.spines["left"].set_bounds(0.0, 100.0)
        axis.set_yticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
    else:
        axis.spines["left"].set_visible(False)
        axis.tick_params(axis="y", left=False, labelleft=False)

    axis.text(
        0.04,
        0.90,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def make_figure() -> plt.Figure:
    """Build and return the two-panel separation-method comparison."""
    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
        sharey=True,
    )

    panel_datasets = (TIC_DATASETS[1], TIC_DATASETS[0])
    for index, (axis, dataset) in enumerate(zip(axes, panel_datasets)):
        path, time_limits, _, color = dataset
        time, intensity = load_chromatogram(path, time_limits)
        normalized = normalize_to_percent(intensity)
        retention_time = float(time[np.argmax(normalized)])

        axis.plot(
            time,
            normalized,
            color=color,
            linewidth=2.0,
            zorder=2,
        )
        axis.axvline(
            retention_time,
            ymin=0.0,
            ymax=0.86,
            color=color,
            linestyle=(0, (4, 3)),
            linewidth=1.5,
            zorder=1,
        )
        axis.text(
            retention_time,
            103.0,
            f"{retention_time:.2f} min",
            color=color,
            fontsize=8.0,
            ha="center",
            va="bottom",
        )

        style_axis(
            axis,
            chr(ord("A") + index),
            time_limits,
            show_y_axis=index == 0,
        )

    return figure


def main() -> None:
    """Render the comparison as a transparent PNG and PDF."""
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
