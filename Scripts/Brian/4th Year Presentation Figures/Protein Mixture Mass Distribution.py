"""Create a standalone figure of the protein mixture mass distributions."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATA_ROOT = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\Protein Mix 21 Injections Set 1\merged"
)
FILE_PREFIX = (
    "20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged"
)
DATASETS = [
    (
        DATA_ROOT / f"{FILE_PREFIX}_1_unidecfiles" / f"{FILE_PREFIX}_1_mass.txt",
        "1 INJ (2,263 ions)",
        "#0072B2",
    ),
    (
        DATA_ROOT / f"{FILE_PREFIX}_10_unidecfiles" / f"{FILE_PREFIX}_10_mass.txt",
        "10 INJ (22,321 ions)",
        "#785EF0",
    ),
    (
        DATA_ROOT / f"{FILE_PREFIX}_21_unidecfiles" / f"{FILE_PREFIX}_21_mass.txt",
        "21 INJ (47,334 ions)",
        "#D55E00",
    ),
]

MASS_LIMITS_KDA = (5.0, 40.0)
MASS_Y_LIMITS = (0.0, 125.0)
STAGGER_STEP = 10.0
FIGURE_SIZE_INCHES = (8.0, 6.0)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Presentations\Brian\4th Year Presentation\Figures"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Protein Mixture Mass Distribution.png"


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
    """Build the Figure 3 overlay with presentation styling."""
    figure, axis = plt.subplots(figsize=FIGURE_SIZE_INCHES, constrained_layout=True)
    offsets = STAGGER_STEP * np.arange(len(DATASETS) - 1, -1, -1)

    for (path, label, color), offset in zip(DATASETS, offsets):
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

    axis.set_xlim(*MASS_LIMITS_KDA)
    axis.set_ylim(*MASS_Y_LIMITS)
    axis.set_xticks(np.arange(5.0, 40.1, 5.0))
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
        bbox_to_anchor=(0.98, 0.98),
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
    try:
        figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    finally:
        plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
