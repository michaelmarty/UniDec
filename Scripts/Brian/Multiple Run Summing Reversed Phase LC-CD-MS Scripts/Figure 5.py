"""Create the six-panel full-range and zoomed mass plots for Figure 5."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


FIGURE_SIZE_INCHES = (7.0, 9.0)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 5.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 5.pdf"

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
        (147.0, 150.0),
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
        (147.0, 150.0),
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
        (142.0, 147.0),
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

PANEL_LETTERS = (("A", "B"), ("C", "D"), ("E", "F"))


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


def style_axis(axis: plt.Axes, show_y_axis: bool) -> None:
    """Apply the shared Figure 5 axis styling."""
    axis.set_ylim(*Y_LIMITS)
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
        axis.set_yticks([])
        axis.tick_params(axis="y", left=False, labelleft=False)


def plot_sample_row(
    full_axis: plt.Axes,
    zoom_axis: plt.Axes,
    sample_label: str,
    zoom_limits: tuple[float, float],
    datasets: list[tuple[Path, str, str]],
    panel_letters: tuple[str, str],
) -> None:
    """Plot one sample's full and zoomed staggered mass distributions."""
    offsets = STAGGER_STEP * np.arange(len(datasets) - 1, -1, -1)
    full_axis.axvspan(
        zoom_limits[0],
        zoom_limits[1],
        facecolor="#2CA02C",
        alpha=0.25,
        linewidth=0.0,
        zorder=0,
    )

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
        full_axis.plot(
            full_x,
            full_y,
            color=color,
            linewidth=1.5,
            label=label,
        )

        zoom_mask = (mass_kda >= zoom_limits[0]) & (mass_kda <= zoom_limits[1])
        zoom_axis.plot(
            mass_kda[zoom_mask],
            normalized[zoom_mask],
            color=color,
            linewidth=1.5,
            label=label,
        )

    full_axis.set_xlim(*FULL_MASS_LIMITS_KDA)
    full_axis.set_xticks(np.arange(25.0, 200.1, 25.0))
    zoom_axis.set_xlim(*zoom_limits)
    zoom_axis.set_xticks(np.arange(zoom_limits[0], zoom_limits[1] + 0.1, 1.0))

    style_axis(full_axis, show_y_axis=True)
    style_axis(zoom_axis, show_y_axis=False)
    add_panel_letter(full_axis, panel_letters[0])
    add_panel_letter(zoom_axis, panel_letters[1])

    full_axis.set_ylabel(
        sample_label,
        rotation=90,
        labelpad=20,
        ha="center",
        va="center",
    )


def make_figure() -> plt.Figure:
    """Build and return the complete six-panel Figure 5."""
    figure, axes = plt.subplots(
        nrows=3,
        ncols=2,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    for row, (sample_label, zoom_limits, datasets) in enumerate(SAMPLES):
        plot_sample_row(
            axes[row, 0],
            axes[row, 1],
            sample_label,
            zoom_limits,
            datasets,
            PANEL_LETTERS[row],
        )

    axes[0, 0].legend(
        loc="upper left",
        bbox_to_anchor=(0.12, 0.66),
        frameon=False,
        fontsize=9.5,
        handlelength=2.0,
        borderaxespad=0.0,
    )
    axes[-1, 0].set_xlabel("Mass (kDa)")
    axes[-1, 1].set_xlabel("Mass (kDa)")

    return figure


def main() -> None:
    """Render Figure 5 as a transparent PNG and PDF."""
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
