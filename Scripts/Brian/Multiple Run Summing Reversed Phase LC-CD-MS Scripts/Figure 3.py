"""Create the staggered, normalized mass-distribution comparison for Figure 3."""

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
        "1 Injection (2,263 ions)",
        "#0072B2",  # Figure 2 blue
    ),
    (
        DATA_ROOT / f"{FILE_PREFIX}_10_unidecfiles" / f"{FILE_PREFIX}_10_mass.txt",
        "10 Injections (22,321 ions)",
        "#785EF0",  # Colorblind-friendly violet; distinct from Figure 2 green
    ),
    (
        DATA_ROOT / f"{FILE_PREFIX}_21_unidecfiles" / f"{FILE_PREFIX}_21_mass.txt",
        "21 Injections (47,334 ions)",
        "#D55E00",  # Figure 2 vermillion
    ),
]

MASS_LIMITS_KDA = (5.0, 40.0)
STAGGER_STEP = 10.0
FIGURE_SIZE_INCHES = (3.25, 2.5)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 3.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 3.pdf"


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

    # Intensities are physically non-negative; clipping also removes tiny
    # negative values that can arise from numerical processing.
    return mass[finite], np.clip(intensity[finite], 0.0, None)


def normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Normalize a distribution so its maximum intensity is exactly 100%."""
    maximum = float(np.max(intensity))
    if maximum <= 0.0:
        return np.zeros_like(intensity)
    return 100.0 * intensity / maximum


def make_figure() -> plt.Figure:
    """Build and return the publication-sized mass-distribution overlay."""
    figure, axis = plt.subplots(
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    # Put the lowest-injection trace on top and the 21-injection trace at the
    # base. Each trace spans a full 0--100% before its offset is applied, giving
    # displayed ranges of 20--120%, 10--110%, and 0--100%, respectively.
    offsets = STAGGER_STEP * np.arange(len(DATASETS) - 1, -1, -1)
    for (path, label, color), offset in zip(DATASETS, offsets):
        mass, intensity = load_mass_distribution(path)
        axis.plot(
            mass / 1000.0,
            normalize_to_percent(intensity) + offset,
            color=color,
            linewidth=1.5,
            label=label,
        )

    axis.set_xlim(*MASS_LIMITS_KDA)
    axis.set_ylim(0.0, 124.0)
    axis.set_xticks(np.arange(5.0, 40.1, 5.0))
    axis.set_yticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
    axis.set_xlabel("Mass (kDa)")

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.5)
    axis.spines["left"].set_bounds(0.0, 100.0)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.tick_params(direction="out", width=1.5, length=6)
    axis.legend(
        loc="upper right",
        frameon=False,
        fontsize=7.0,
        handlelength=1.5,
        borderaxespad=0.3,
    )

    return figure


def main() -> None:
    """Render Figure 3 as both PNG and PDF."""
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
    # Avoid bbox_inches="tight" so the saved files retain the requested size.
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
