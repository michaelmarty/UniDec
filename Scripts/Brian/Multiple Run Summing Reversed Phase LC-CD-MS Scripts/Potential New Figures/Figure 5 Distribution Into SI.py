"""Create the three full-range Figure 5 mass distributions for the SI.

These are the former Figure 5 panels A, C, and E, relabeled A, B, and C for
their new single-column position.
"""

from importlib import util
from pathlib import Path
from types import ModuleType

import matplotlib.pyplot as plt
import numpy as np


FIGURE_SIZE_INCHES = (3.5, 9.0)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 5 Distribution Into SI.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 5 Distribution Into SI.pdf"

SOURCE_DIRECTORY = OUTPUT_DIRECTORY.parent
PANEL_LETTERS = ("A", "B", "C")


def load_figure_module(filename: str, module_name: str) -> ModuleType:
    """Load a source figure script from the parent directory."""
    source_path = SOURCE_DIRECTORY / filename
    specification = util.spec_from_file_location(module_name, source_path)
    if specification is None or specification.loader is None:
        raise ImportError(f"Could not load figure source: {source_path}")

    module = util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


FIGURE_5 = load_figure_module("Figure 5.py", "figure_5_source_for_si")


def plot_full_mass_panel(
    axis: plt.Axes,
    panel_letter: str,
    sample_label: str,
    datasets: list[tuple[Path, str, str]],
) -> None:
    """Plot one full-range staggered mass-distribution panel."""
    offsets = FIGURE_5.STAGGER_STEP * np.arange(len(datasets) - 1, -1, -1)
    for (path, label, color), offset in zip(datasets, offsets):
        mass_da, intensity = FIGURE_5.load_mass_distribution(path)
        mass_kda = mass_da / 1000.0
        normalized = FIGURE_5.normalize_to_percent(intensity) + offset
        full_mask = (
            (mass_kda >= FIGURE_5.FULL_MASS_LIMITS_KDA[0])
            & (mass_kda <= FIGURE_5.FULL_MASS_LIMITS_KDA[1])
        )
        full_x, full_y = FIGURE_5.minmax_downsample(
            mass_kda[full_mask], normalized[full_mask]
        )
        axis.plot(
            full_x,
            full_y,
            color=color,
            linewidth=1.5,
            label=label,
        )

    axis.set_xlim(*FIGURE_5.FULL_MASS_LIMITS_KDA)
    axis.set_xticks(np.arange(25.0, 200.1, 25.0))
    FIGURE_5.style_axis(axis, show_y_axis=True)
    FIGURE_5.add_panel_letter(axis, panel_letter)
    axis.set_ylabel(
        sample_label,
        rotation=90,
        labelpad=20,
        ha="center",
        va="center",
    )


def make_figure() -> plt.Figure:
    """Build and return the three-panel SI figure."""
    figure, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
        sharex=True,
    )

    for axis, panel_letter, sample in zip(axes, PANEL_LETTERS, FIGURE_5.SAMPLES):
        sample_label, _, datasets = sample
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
    """Render the SI figure as a transparent PNG and PDF."""
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
