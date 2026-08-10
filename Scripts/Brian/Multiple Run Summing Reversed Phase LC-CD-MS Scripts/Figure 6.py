"""Create the TIC and mass-distribution comparisons used in Figure 6."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATA_ROOT = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT"
    r"\Amgen IgM-like Fusion Protein 21 Injections Set 1\merged"
)
FILE_PREFIX = "20260424_BHT_agilentc8column_AmgenIgM_CDMS_set1_0o1mgmLs_merged"

ENSEMBLE_TIC = Path(
    r"Z:\Group Share\BHT\Temp\Paper\Ensemble Amgen IgM-like Fusion Triplicate"
    r"\IgM-like Fusion Protein TIC.txt"
)

TIC_DATASETS = [
    (ENSEMBLE_TIC, "Ensemble MS", "#333333"),
    (
        DATA_ROOT / f"{FILE_PREFIX}_1_unidecfiles" / "TIC_chrom.txt",
        "1-injection CD-MS",
        "#0072B2",
    ),
    (
        DATA_ROOT / f"{FILE_PREFIX}_21_unidecfiles" / "TIC_decon.txt",
        "21-injection CD-MS",
        "#D55E00",
    ),
]

MASS_DATASETS = [
    (
        DATA_ROOT
        / f"{FILE_PREFIX}_1_unidecfiles"
        / f"{FILE_PREFIX}_1_mass.txt",
        "1 Injection",
        "#0072B2",
    ),
    (
        DATA_ROOT
        / f"{FILE_PREFIX}_10_unidecfiles"
        / f"{FILE_PREFIX}_10_mass.txt",
        "10 Injections",
        "#785EF0",
    ),
    (
        DATA_ROOT
        / f"{FILE_PREFIX}_21_unidecfiles"
        / f"{FILE_PREFIX}_21_mass.txt",
        "21 Injections",
        "#D55E00",
    ),
]

TIME_LIMITS = (0.0, 5.5)
MASS_LIMITS_KDA = (0.0, 2000.0)
STAGGER_STEP = 10.0
Y_LIMITS = (0.0, 124.0)
FIGURE_SIZE_INCHES = (7.0, 3.25)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 6.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 6.pdf"


def load_distribution(
    path: Path,
    x_limits: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Load and validate a two-column x/intensity distribution."""
    if not path.is_file():
        raise FileNotFoundError(f"Distribution not found: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in: {path}")

    x_values = data[:, 0]
    intensity = data[:, 1]
    valid = np.isfinite(x_values) & np.isfinite(intensity)
    if x_limits is not None:
        valid &= (x_values >= x_limits[0]) & (x_values <= x_limits[1])
    if not np.any(valid):
        raise ValueError(f"No finite values in the requested range: {path}")

    return x_values[valid], np.clip(intensity[valid], 0.0, None)


def minmax_normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Min-max normalize a chromatogram to exactly 0--100%."""
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


def style_axis(
    axis: plt.Axes,
    panel_letter: str,
    show_y_axis: bool,
    letter_x: float,
    letter_y: float,
    letter_alignment: str,
) -> None:
    """Apply shared publication styling to one Figure 6 panel."""
    axis.set_ylim(*Y_LIMITS)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.tick_params(direction="out", width=1.5, length=6)

    if show_y_axis:
        axis.spines["left"].set_linewidth(1.5)
        axis.spines["left"].set_bounds(0.0, 100.0)
        axis.set_yticks([0.0, 50.0, 100.0], labels=["0%", "%", "100%"])
    else:
        axis.spines["left"].set_visible(False)
        axis.set_yticks([])
        axis.tick_params(axis="y", left=False, labelleft=False)

    axis.text(
        letter_x,
        letter_y,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha=letter_alignment,
        va="top",
        zorder=10,
    )


def make_figure() -> plt.Figure:
    """Build and return the complete two-panel Figure 6."""
    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )
    offsets = STAGGER_STEP * np.arange(2, -1, -1)

    retention_time = None
    for (path, label, color), offset in zip(TIC_DATASETS, offsets):
        time, intensity = load_distribution(path, TIME_LIMITS)
        if label == "21-injection CD-MS":
            retention_time = float(time[np.argmax(intensity)])
        axes[0].plot(
            time,
            minmax_normalize_to_percent(intensity) + offset,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=2,
        )

    axes[0].set_xlim(*TIME_LIMITS)
    axes[0].set_xticks(np.arange(0.0, 5.1, 1.0))
    axes[0].set_xlabel("Time (min)")
    style_axis(
        axes[0],
        "A",
        show_y_axis=True,
        letter_x=0.985,
        letter_y=0.97,
        letter_alignment="right",
    )
    if retention_time is not None:
        axes[0].axvline(
            retention_time,
            ymin=0.0,
            ymax=0.86,
            color="#000000",
            linestyle=(0, (4, 3)),
            linewidth=1.5,
            zorder=3,
        )
        axes[0].text(
            retention_time - 0.08,
            0.87,
            f"{retention_time:.2f} min",
            transform=axes[0].get_xaxis_transform(),
            color="#000000",
            fontsize=7.0,
            ha="right",
            va="bottom",
            zorder=4,
        )
    axes[0].legend(
        loc="upper right",
        bbox_to_anchor=(1.0, 0.80),
        frameon=False,
        fontsize=9.5,
        handlelength=2.0,
        borderaxespad=0.0,
    )

    for (path, label, color), offset in zip(MASS_DATASETS, offsets):
        mass_da, intensity = load_distribution(path)
        axes[1].plot(
            mass_da / 1000.0,
            max_normalize_to_percent(intensity) + offset,
            color=color,
            linewidth=1.5,
            label=label,
        )

    axes[1].set_xlim(*MASS_LIMITS_KDA)
    axes[1].set_xticks(np.arange(0.0, 2000.1, 500.0))
    axes[1].set_xlabel("Mass (kDa)")
    style_axis(
        axes[1],
        "B",
        show_y_axis=False,
        letter_x=0.985,
        letter_y=0.97,
        letter_alignment="right",
    )
    axes[1].legend(
        loc="upper left",
        bbox_to_anchor=(0.02, 0.80),
        frameon=False,
        fontsize=9.5,
        handlelength=2.0,
        borderaxespad=0.0,
    )

    return figure


def main() -> None:
    """Render Figure 6 as a transparent PNG and PDF."""
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
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
