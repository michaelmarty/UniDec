"""Create a standalone figure of the Waters mAb mass distribution."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


FIGURE_SIZE_INCHES = (8.0, 6.0)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Presentations\Brian\4th Year Presentation\Figures"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Waters mAb Mass Distribution.png"

MASS_LIMITS_KDA = (147.0, 150.0)
MASS_Y_LIMITS = (0.0, 125.0)
STAGGER_STEP = 10.0
GLYCOFORM_MASS_DA = 148868.0
GLYCOFORM_MARKER_WINDOW_DA = 50.0
GLYCOFORM_MARKER_OFFSET = 4.0

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


def mass_path(injections: int) -> Path:
    """Construct a UniDec mass-distribution path for one injection count."""
    return (
        WATERS_ROOT
        / f"{WATERS_PREFIX}_{injections}_unidecfiles"
        / f"{WATERS_PREFIX}_{injections}_mass.txt"
    )


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


def normalize_mass_distribution(intensity: np.ndarray) -> np.ndarray:
    """Normalize a mass distribution to a maximum of exactly 100%."""
    maximum = float(np.max(intensity))
    if maximum <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * intensity / maximum


def make_figure() -> plt.Figure:
    """Build and return the Waters mAb mass-distribution figure."""
    figure, axis = plt.subplots(figsize=FIGURE_SIZE_INCHES, constrained_layout=True)
    datasets = [
        (mass_path(1), "1 INJ", ONE_INJECTION_COLOR),
        (mass_path(10), "10 INJ", TEN_INJECTION_COLOR),
        (mass_path(21), "21 INJ", TWENTY_ONE_INJECTION_COLOR),
    ]
    offsets = STAGGER_STEP * np.arange(len(datasets) - 1, -1, -1)
    marker_peak_height = -np.inf

    for (path, label, color), offset in zip(datasets, offsets):
        mass_da, intensity = load_mass_distribution(path)
        mass_kda = mass_da / 1000.0
        normalized = normalize_mass_distribution(intensity) + offset
        zoom_mask = (
            (mass_kda >= MASS_LIMITS_KDA[0])
            & (mass_kda <= MASS_LIMITS_KDA[1])
        )
        axis.plot(
            mass_kda[zoom_mask],
            normalized[zoom_mask],
            color=color,
            linewidth=2.0,
            label=label,
        )

        marker_mask = zoom_mask & (
            np.abs(mass_da - GLYCOFORM_MASS_DA) <= GLYCOFORM_MARKER_WINDOW_DA
        )
        if np.any(marker_mask):
            marker_peak_height = max(
                marker_peak_height,
                float(np.max(normalized[marker_mask])),
            )

    if not np.isfinite(marker_peak_height):
        raise ValueError(
            f"No data found near the selected {GLYCOFORM_MASS_DA:,.0f} Da glycoform"
        )
    axis.text(
        GLYCOFORM_MASS_DA / 1000.0,
        marker_peak_height + GLYCOFORM_MARKER_OFFSET,
        "*",
        color="black",
        fontsize=24,
        fontweight="bold",
        ha="center",
        va="bottom",
        zorder=5,
    )

    axis.set_xlim(*MASS_LIMITS_KDA)
    axis.set_ylim(*MASS_Y_LIMITS)
    axis.set_xticks(np.arange(MASS_LIMITS_KDA[0], MASS_LIMITS_KDA[1] + 0.1, 1.0))
    axis.set_xlabel("Mass (kDa)")
    axis.set_yticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(2.0)
    axis.spines["left"].set_bounds(0.0, 100.0)
    axis.spines["bottom"].set_linewidth(2.0)
    axis.tick_params(direction="out", width=2.0, length=6)
    axis.legend(
        loc="upper right",
        bbox_to_anchor=(0.98, 0.70),
        frameon=False,
        handlelength=1.5,
        borderaxespad=0.0,
    )
    return figure


def main() -> None:
    """Render the mass distribution as a transparent 300 DPI PNG."""
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
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
