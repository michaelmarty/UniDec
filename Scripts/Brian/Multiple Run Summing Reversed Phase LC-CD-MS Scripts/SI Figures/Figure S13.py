"""Create the three full-range mass distributions for Figure S13.

These are the former Figure 5 panels A, C, and E, relabeled A, B, and C for
their new single-column position.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


FIGURE_SIZE_INCHES = (3.25, 7.5)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S13"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S13.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S13.pdf"

FULL_MASS_LIMITS_KDA = (10.0, 200.0)
STAGGER_STEP = 10.0
Y_LIMITS = (0.0, 150.0)
MAX_FULL_RANGE_POINTS = 12_000

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

PANEL_LETTERS = ("A", "B", "C")


def mass_path(root: Path, prefix: str, injections: int) -> Path:
    """Construct a UniDec mass-distribution path for an injection count."""
    return (
        root
        / f"{prefix}_{injections}_unidecfiles"
        / f"{prefix}_{injections}_mass.txt"
    )


SAMPLES = [
    (
        "Waters Intact mAb\nMass Check Standard",
        [
            (
                mass_path(WATERS_ROOT, WATERS_PREFIX, 1),
                "1 Injection",
                ONE_INJECTION_COLOR,
            ),
            (
                mass_path(WATERS_ROOT, WATERS_PREFIX, 10),
                "10 Injections",
                TEN_INJECTION_COLOR,
            ),
            (
                mass_path(WATERS_ROOT, WATERS_PREFIX, 21),
                "21 Injections",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
    (
        "Human IgG",
        [
            (
                mass_path(HUMAN_ROOT, HUMAN_PREFIX, 1),
                "1 Injection",
                ONE_INJECTION_COLOR,
            ),
            (
                mass_path(HUMAN_ROOT, HUMAN_PREFIX, 10),
                "10 Injections",
                TEN_INJECTION_COLOR,
            ),
            (
                mass_path(HUMAN_ROOT, HUMAN_PREFIX, 21),
                "21 Injections",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
    (
        "Rabbit IgG",
        [
            (
                mass_path(RABBIT_ROOT, RABBIT_PREFIX, 1),
                "1 Injection",
                ONE_INJECTION_COLOR,
            ),
            (
                mass_path(RABBIT_ROOT, RABBIT_PREFIX, 10),
                "10 Injections",
                TEN_INJECTION_COLOR,
            ),
            (
                mass_path(RABBIT_ROOT, RABBIT_PREFIX, 21),
                "21 Injections",
                TWENTY_ONE_INJECTION_COLOR,
            ),
        ],
    ),
]


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


def normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Normalize a distribution so its maximum intensity is exactly 100%."""
    maximum = float(np.max(intensity))
    if maximum <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * intensity / maximum


def minmax_downsample(
    x_values: np.ndarray,
    y_values: np.ndarray,
    max_points: int = MAX_FULL_RANGE_POINTS,
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


def add_panel_letter(axis: plt.Axes, panel_letter: str) -> None:
    """Add a bold panel identifier at the upper-left of an axis."""
    axis.text(
        0.015,
        0.64,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def style_axis(axis: plt.Axes) -> None:
    """Apply the shared Figure S13 axis styling."""
    axis.set_ylim(*Y_LIMITS)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.spines["left"].set_linewidth(1.5)
    axis.spines["left"].set_bounds(0.0, 100.0)
    axis.set_yticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
    axis.tick_params(direction="out", width=1.5, length=6)


def plot_full_mass_panel(
    axis: plt.Axes,
    panel_letter: str,
    sample_label: str,
    datasets: list[tuple[Path, str, str]],
) -> None:
    """Plot one full-range staggered mass-distribution panel."""
    offsets = STAGGER_STEP * np.arange(len(datasets) - 1, -1, -1)
    for (path, label, color), offset in zip(datasets, offsets):
        mass_da, intensity = load_mass_distribution(path)
        mass_kda = mass_da / 1000.0
        normalized = normalize_to_percent(intensity) + offset
        full_mask = (
            (mass_kda >= FULL_MASS_LIMITS_KDA[0])
            & (mass_kda <= FULL_MASS_LIMITS_KDA[1])
        )
        full_x, full_y = minmax_downsample(
            mass_kda[full_mask], normalized[full_mask]
        )
        axis.plot(
            full_x,
            full_y,
            color=color,
            linewidth=1.5,
            label=label,
        )

    axis.set_xlim(*FULL_MASS_LIMITS_KDA)
    axis.set_xticks(np.arange(50.0, 200.1, 50.0))
    style_axis(axis)
    add_panel_letter(axis, panel_letter)
    y_minimum, y_maximum = axis.get_ylim()
    spine_minimum, spine_maximum = axis.spines["left"].get_bounds()
    spine_center = (spine_minimum + spine_maximum) / 2.0
    label_y = (spine_center - y_minimum) / (y_maximum - y_minimum)
    axis.set_ylabel(
        sample_label,
        rotation=90,
        labelpad=20,
        ha="center",
        va="center",
    )
    axis.yaxis.label.set_y(label_y)


def make_figure() -> plt.Figure:
    """Build and return the three-panel SI figure."""
    figure, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
        sharex=True,
    )
    figure.get_layout_engine().set(h_pad=0.01, hspace=0.0)

    for axis, panel_letter, sample in zip(axes, PANEL_LETTERS, SAMPLES):
        sample_label, datasets = sample
        plot_full_mass_panel(axis, panel_letter, sample_label, datasets)

    axes[0].legend(
        loc="upper left",
        bbox_to_anchor=(0.12, 0.66),
        frameon=False,
        fontsize=9.5,
        handlelength=2.0,
        borderaxespad=0.0,
    )
    axes[-1].set_xlabel("Mass (kDa)")

    return figure


def main() -> None:
    """Render Figure S13 as a transparent PNG and PDF."""
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
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
