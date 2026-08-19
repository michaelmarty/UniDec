"""Create a standalone combined TIC and zoomed mass-distribution figure.

The three rows contain the Waters, human, and rabbit samples, respectively.
TIC panels are in the left column and zoomed mass-distribution panels are in
the right column, giving the row-major panel order A--F.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


FIGURE_SIZE_INCHES = (7.0, 9.0)
OUTPUT_DPI = 300
MINIMUM_FONT_SIZE = 10
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 4 and 5 Combined.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 4 and 5 Combined.pdf"

TIME_LIMITS = (0.0, 5.5)
TIC_Y_LIMITS = (0.0, 180.0)
MASS_Y_LIMITS = (0.0, 150.0)
STAGGER_STEP = 10.0
PANEL_LETTER_X = 0.015
PANEL_LETTER_Y = 0.85
PANEL_LETTERS = (("A", "B"), ("C", "D"), ("E", "F"))

ENSEMBLE_COLOR = "#333333"
ONE_INJECTION_COLOR = "#0072B2"
TEN_INJECTION_COLOR = "#785EF0"
TWENTY_ONE_INJECTION_COLOR = "#D55E00"

WATERS_ROOT = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\Waters mAb 21 Injections Set 1\merged"
)
WATERS_PREFIX = (
    "20260330_BHT_agilentc8column_watersmAb_CMDS_20msINJ_set1_"
    "0o1mgmLs_merged"
)

HUMAN_ROOT = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT"
    r"\Amgen Human IgG 21 Injections Set 1\merged"
)
HUMAN_PREFIX = (
    "20260331_BHT_agilentc8column_AmgenhumanmAb_PRO35304_CMDS_20msINJ_"
    "set1_0o1mgmLs_merged"
)

RABBIT_ROOT = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT"
    r"\Amgen Rabbit IgG 21 Injections Set 1\merged"
)
RABBIT_PREFIX = (
    "20260401_BHT_agilentc8column_AmgenrabbitmAb_EPA13724_CMDS_20msINJ_"
    "set1_0o1mgmLs_merged"
)


def tic_path(
    root: Path,
    prefix: str,
    injections: int,
    filename: str,
) -> Path:
    """Construct a UniDec TIC path for one injection count."""
    return root / f"{prefix}_{injections}_unidecfiles" / filename


def mass_path(root: Path, prefix: str, injections: int) -> Path:
    """Construct a UniDec mass-distribution path for one injection count."""
    return (
        root
        / f"{prefix}_{injections}_unidecfiles"
        / f"{prefix}_{injections}_mass.txt"
    )


TIC_PANELS = [
    (
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
                tic_path(WATERS_ROOT, WATERS_PREFIX, 1, "TIC_chrom.txt"),
                "1 INJ CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                tic_path(WATERS_ROOT, WATERS_PREFIX, 21, "TIC_decon.txt"),
                "21 INJ CD-MS",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
    (
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
                tic_path(HUMAN_ROOT, HUMAN_PREFIX, 1, "TIC_chrom.txt"),
                "1 INJ CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                tic_path(HUMAN_ROOT, HUMAN_PREFIX, 21, "TIC_decon.txt"),
                "21 INJ CD-MS",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
    (
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
                tic_path(RABBIT_ROOT, RABBIT_PREFIX, 1, "TIC_chrom.txt"),
                "1 INJ CD-MS",
                ONE_INJECTION_COLOR,
            ),
            (
                tic_path(RABBIT_ROOT, RABBIT_PREFIX, 21, "TIC_decon.txt"),
                "21 INJ CD-MS",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
]

MASS_SAMPLES = [
    (
        "Waters Intact mAb\nMass Check Standard",
        (147.0, 150.0),
        [
            (
                mass_path(WATERS_ROOT, WATERS_PREFIX, 1),
                "1 INJ",
                ONE_INJECTION_COLOR,
            ),
            (
                mass_path(WATERS_ROOT, WATERS_PREFIX, 10),
                "10 INJ",
                TEN_INJECTION_COLOR,
            ),
            (
                mass_path(WATERS_ROOT, WATERS_PREFIX, 21),
                "21 INJ",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
    (
        "Human IgG",
        (147.0, 150.0),
        [
            (
                mass_path(HUMAN_ROOT, HUMAN_PREFIX, 1),
                "1 INJ",
                ONE_INJECTION_COLOR,
            ),
            (
                mass_path(HUMAN_ROOT, HUMAN_PREFIX, 10),
                "10 INJ",
                TEN_INJECTION_COLOR,
            ),
            (
                mass_path(HUMAN_ROOT, HUMAN_PREFIX, 21),
                "21 INJ",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
    (
        "Rabbit IgG",
        (142.0, 147.0),
        [
            (
                mass_path(RABBIT_ROOT, RABBIT_PREFIX, 1),
                "1 INJ",
                ONE_INJECTION_COLOR,
            ),
            (
                mass_path(RABBIT_ROOT, RABBIT_PREFIX, 10),
                "10 INJ",
                TEN_INJECTION_COLOR,
            ),
            (
                mass_path(RABBIT_ROOT, RABBIT_PREFIX, 21),
                "21 INJ",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
]


def load_chromatogram(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a two-column chromatogram within the shared time window."""
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


def load_mass_distribution(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load and validate a two-column mass/intensity distribution."""
    if not path.is_file():
        raise FileNotFoundError(f"Mass distribution not found: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in: {path}")

    mass = data[:, 0]
    intensity = data[:, 1]
    finite = np.isfinite(mass) & np.isfinite(intensity)
    if not np.any(finite):
        raise ValueError(f"No finite mass/intensity pairs in: {path}")

    return mass[finite], np.clip(intensity[finite], 0.0, None)


def normalize_chromatogram(intensity: np.ndarray) -> np.ndarray:
    """Min-max normalize a chromatogram to exactly 0--100%."""
    minimum = float(np.min(intensity))
    maximum = float(np.max(intensity))
    span = maximum - minimum
    if span <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * (intensity - minimum) / span


def normalize_mass_distribution(intensity: np.ndarray) -> np.ndarray:
    """Normalize a mass distribution to a maximum of exactly 100%."""
    maximum = float(np.max(intensity))
    if maximum <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * intensity / maximum


def add_panel_letter(axis: plt.Axes, panel_letter: str) -> None:
    """Place every panel letter at one shared axes-relative position."""
    axis.text(
        PANEL_LETTER_X,
        PANEL_LETTER_Y,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def style_tic_axis(
    axis: plt.Axes,
    panel_letter: str,
    sample_label: str,
) -> None:
    """Apply the combined figure's standalone TIC styling."""
    axis.set_xlim(*TIME_LIMITS)
    axis.set_ylim(*TIC_Y_LIMITS)
    axis.set_ylabel(
        sample_label,
        fontsize=13,
        rotation=90,
        labelpad=22,
        ha="center",
        va="center",
    )
    axis.yaxis.set_label_coords(-0.27, 50.0 / TIC_Y_LIMITS[1])
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
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
    add_panel_letter(axis, panel_letter)


def style_mass_axis(axis: plt.Axes, panel_letter: str) -> None:
    """Apply the combined figure's standalone zoomed-mass styling."""
    axis.set_ylim(*MASS_Y_LIMITS)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.set_yticks([])
    axis.tick_params(
        direction="out",
        width=1.5,
        length=6,
        labelsize=MINIMUM_FONT_SIZE,
    )
    axis.tick_params(axis="y", left=False, labelleft=False)
    add_panel_letter(axis, panel_letter)


def plot_tic_panel(
    axis: plt.Axes,
    panel_letter: str,
    sample_label: str,
    datasets: list[tuple[Path, str, str]],
) -> None:
    """Plot one staggered TIC panel."""
    offsets = STAGGER_STEP * np.arange(len(datasets) - 1, -1, -1)
    for (path, label, color), offset in zip(datasets, offsets):
        time, intensity = load_chromatogram(path)
        axis.plot(
            time,
            normalize_chromatogram(intensity) + offset,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=2,
        )

    style_tic_axis(axis, panel_letter, sample_label)


def plot_zoomed_mass_panel(
    axis: plt.Axes,
    panel_letter: str,
    zoom_limits: tuple[float, float],
    datasets: list[tuple[Path, str, str]],
) -> None:
    """Plot one staggered zoomed mass-distribution panel."""
    offsets = STAGGER_STEP * np.arange(len(datasets) - 1, -1, -1)
    for (path, label, color), offset in zip(datasets, offsets):
        mass_da, intensity = load_mass_distribution(path)
        mass_kda = mass_da / 1000.0
        normalized = normalize_mass_distribution(intensity) + offset
        zoom_mask = (mass_kda >= zoom_limits[0]) & (mass_kda <= zoom_limits[1])
        axis.plot(
            mass_kda[zoom_mask],
            normalized[zoom_mask],
            color=color,
            linewidth=1.5,
            label=label,
        )

    axis.set_xlim(*zoom_limits)
    axis.set_xticks(np.arange(zoom_limits[0], zoom_limits[1] + 0.1, 1.0))
    style_mass_axis(axis, panel_letter)


def make_figure() -> plt.Figure:
    """Build and return the six-panel combined figure."""
    figure, axes = plt.subplots(
        nrows=3,
        ncols=2,
        figsize=FIGURE_SIZE_INCHES,
        gridspec_kw={"height_ratios": [1, 1, 1]},
        constrained_layout=True,
    )

    for row, (tic_panel, mass_sample, letters) in enumerate(
        zip(TIC_PANELS, MASS_SAMPLES, PANEL_LETTERS)
    ):
        tic_sample_label, tic_datasets = tic_panel
        mass_sample_label, zoom_limits, mass_datasets = mass_sample
        if tic_sample_label != mass_sample_label:
            raise ValueError(
                "TIC and mass-distribution sample rows do not match: "
                f"{tic_sample_label!r} != {mass_sample_label!r}"
            )

        plot_tic_panel(axes[row, 0], letters[0], tic_sample_label, tic_datasets)
        plot_zoomed_mass_panel(
            axes[row, 1],
            letters[1],
            zoom_limits,
            mass_datasets,
        )

    axes[0, 0].legend(
        loc="upper right",
        bbox_to_anchor=(0.98, 0.70),
        frameon=False,
        fontsize=MINIMUM_FONT_SIZE,
        handlelength=1.5,
        borderaxespad=0.0,
    )
    axes[0, 1].legend(
        loc="upper right",
        bbox_to_anchor=(0.98, 0.70),
        frameon=False,
        fontsize=MINIMUM_FONT_SIZE,
        handlelength=1.5,
        borderaxespad=0.0,
    )

    axes[-1, 0].set_xlabel("Time (min)", fontsize=MINIMUM_FONT_SIZE)
    axes[-1, 0].set_xticks(np.arange(0.0, 5.1, 1.0))
    axes[-1, 0].tick_params(labelsize=MINIMUM_FONT_SIZE)
    axes[-1, 1].set_xlabel("Mass (kDa)", fontsize=MINIMUM_FONT_SIZE)
    axes[-1, 1].tick_params(labelsize=MINIMUM_FONT_SIZE)

    return figure


def main() -> None:
    """Render the combined figure as a transparent PNG and PDF."""
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
