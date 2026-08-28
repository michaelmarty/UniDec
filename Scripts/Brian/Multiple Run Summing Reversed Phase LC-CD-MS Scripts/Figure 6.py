"""Compare Waters mAb TICs and mass distributions for Figure 6."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin\MartyLab - General"
    r"\Papers\Muliple Run Summing Reversed Phase LC-CD-MS\Figures for Manuscript"
    r"\Figure 6"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 6.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 6.pdf"
FIGURE_SIZE_INCHES = (7, 6.5)
MINIMUM_FONT_SIZE = 10

TIME_LIMITS = (0.0, 5.5)
MASS_LIMITS_KDA = (147.0, 150.0)
Y_LIMITS = (0.0, 115.0)
MAX_MASS_POINTS = 12_000

STANDARD_METHOD_COLOR = "#5B6770"
HIGH_THROUGHPUT_COLOR = "#EE6677"

TIC_DATASETS = [
    (
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Waters mAb 21 Injections Set 1\merged"
            r"\20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_"
            r"0o1mgmLs_merged_21_unidecfiles\TIC_decon.txt"
        ),
        (0.0, 5.5),
        "8.1-min method (Agilent C8)",
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
        "3.0-min method (HALO C4)",
        HIGH_THROUGHPUT_COLOR,
    ),
]

MASS_DATASETS = [
    (
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Waters mAb 21 Injections Set 1\merged"
            r"\20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_"
            r"0o1mgmLs_merged_21_unidecfiles"
            r"\20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_"
            r"0o1mgmLs_merged_21_mass.txt"
        ),
        "8-min method (Agilent C8)",
        STANDARD_METHOD_COLOR,
    ),
    (
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
            r"\Waters mAb 21 Injections Set 1\merged"
            r"\20260501_BHT_haloc4column_watersmAb_set1_CDMS_0o1mgmL_"
            r"merged_21_unidecfiles"
            r"\20260501_BHT_haloc4column_watersmAb_set1_CDMS_0o1mgmL_"
            r"merged_21_mass.txt"
        ),
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


def load_mass_distribution(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load and validate one two-column mass-distribution file."""
    if not path.is_file():
        raise FileNotFoundError(f"Mass distribution not found: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in: {path}")

    mass_da = data[:, 0]
    intensity = data[:, 1]
    valid = np.isfinite(mass_da) & np.isfinite(intensity)
    if not np.any(valid):
        raise ValueError(f"No finite mass/intensity pairs in: {path}")

    return mass_da[valid], np.clip(intensity[valid], 0.0, None)


def normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Min-max normalize a TIC to exactly 0--100%."""
    minimum = float(np.min(intensity))
    maximum = float(np.max(intensity))
    span = maximum - minimum
    if span <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * (intensity - minimum) / span


def max_normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Normalize a mass distribution so its maximum is exactly 100%."""
    maximum = float(np.max(intensity))
    if maximum <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * intensity / maximum


def minmax_downsample(
    x_values: np.ndarray,
    y_values: np.ndarray,
    max_points: int = MAX_MASS_POINTS,
) -> tuple[np.ndarray, np.ndarray]:
    """Reduce dense lines while retaining each display bin's extrema."""
    if len(x_values) <= max_points:
        return x_values, y_values

    block_size = int(np.ceil(len(x_values) / (max_points // 2)))
    indices: list[int] = [0]
    for start in range(0, len(x_values), block_size):
        stop = min(start + block_size, len(x_values))
        block = y_values[start:stop]
        indices.extend((start + int(np.argmin(block)), start + int(np.argmax(block))))
    indices.append(len(x_values) - 1)

    ordered_indices = np.unique(np.asarray(indices, dtype=int))
    return x_values[ordered_indices], y_values[ordered_indices]


def style_axis(axis: plt.Axes, panel_letter: str) -> None:
    """Apply shared publication styling and a panel identifier."""
    axis.set_ylim(*Y_LIMITS)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.spines["left"].set_linewidth(1.5)
    axis.spines["left"].set_bounds(0.0, 100.0)
    axis.set_yticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
    axis.tick_params(direction="out", width=1.5, length=6)
    axis.text(
        0.015,
        0.92,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def style_tic_axis(axis: plt.Axes) -> None:
    """Apply TIC-specific limits and labels."""
    style_axis(axis, "A")
    axis.set_xlim(*TIME_LIMITS)
    axis.set_xticks(np.arange(0.0, 5.1, 1.0))
    axis.set_xlabel("Time (min)")


def style_mass_axis(axis: plt.Axes) -> None:
    """Apply mass-distribution-specific limits and labels."""
    style_axis(axis, "B")
    axis.set_xlim(*MASS_LIMITS_KDA)
    axis.set_xticks(np.arange(147.0, 150.1, 1.0))
    axis.set_xlabel("Mass (kDa)")


def make_figure() -> plt.Figure:
    """Build and return the TIC and mass-distribution comparison."""
    figure, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )
    tic_axis, mass_axis = axes

    for path, time_limits, label, color in TIC_DATASETS:
        time, intensity = load_chromatogram(path, time_limits)
        normalized = normalize_to_percent(intensity)
        retention_time = float(time[np.argmax(normalized)])

        tic_axis.plot(
            time,
            normalized,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=2,
        )
        tic_axis.axvline(
            retention_time,
            ymin=0.0,
            ymax=0.86,
            color=color,
            linestyle=(0, (4, 3)),
            linewidth=1.5,
            zorder=1,
        )
        tic_axis.text(
            retention_time,
            103.0,
            f"{retention_time:.2f} min",
            color=color,
            fontsize=MINIMUM_FONT_SIZE,
            ha="center",
            va="bottom",
        )

    style_tic_axis(tic_axis)
    tic_axis.legend(
        loc="upper right",
        frameon=False,
        fontsize=MINIMUM_FONT_SIZE,
        handlelength=2.0,
    )

    for path, _, color in MASS_DATASETS:
        mass_da, intensity = load_mass_distribution(path)
        mass_kda = mass_da / 1000.0
        normalized = max_normalize_to_percent(intensity)
        display_mask = (
            (mass_kda >= MASS_LIMITS_KDA[0])
            & (mass_kda <= MASS_LIMITS_KDA[1])
        )
        display_mass, display_intensity = minmax_downsample(
            mass_kda[display_mask], normalized[display_mask]
        )
        mass_axis.plot(
            display_mass,
            display_intensity,
            color=color,
            linewidth=1.5,
        )

    style_mass_axis(mass_axis)

    return figure


def main() -> None:
    """Render Figure 6 as a transparent PNG and PDF."""
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
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
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
