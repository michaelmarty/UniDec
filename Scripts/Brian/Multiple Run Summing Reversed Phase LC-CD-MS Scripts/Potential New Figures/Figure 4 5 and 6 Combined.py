"""Combine Figures 4, 5, and 6 into a consistent eight-panel figure.

The first three rows reproduce the combined Figure 4/5 panels A--F. Figure 6
is appended as the fourth row, with its TIC and mass distribution relabeled G
and H, respectively.
"""

from importlib import util
from pathlib import Path
from types import ModuleType

import matplotlib.pyplot as plt
import numpy as np


FIGURE_SIZE_INCHES = (7.0, 12.0)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 4 5 and 6 Combined.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 4 5 and 6 Combined.pdf"

SOURCE_DIRECTORY = OUTPUT_DIRECTORY.parent
PANEL_LETTERS = (("A", "B"), ("C", "D"), ("E", "F"))
FIGURE_6_SAMPLE_LABEL = "IgM-like Fusion\nProtein"


def load_figure_module(path: Path, module_name: str) -> ModuleType:
    """Load a figure script without requiring import-safe filenames."""
    specification = util.spec_from_file_location(module_name, path)
    if specification is None or specification.loader is None:
        raise ImportError(f"Could not load figure source: {path}")

    module = util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


BASE_COMBINED = load_figure_module(
    OUTPUT_DIRECTORY / "Figure 4 and 5 Combined.py",
    "figure_4_and_5_combined_source",
)
FIGURE_6 = load_figure_module(
    SOURCE_DIRECTORY / "Figure 6.py",
    "figure_6_source",
)


def add_panel_letter(axis: plt.Axes, panel_letter: str) -> None:
    """Place a panel letter consistently with the combined A--F panels."""
    axis.text(
        0.015,
        0.97,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def plot_figure_6_mass_panel(axis: plt.Axes) -> None:
    """Plot Figure 6B as panel H using the combined mass-panel styling."""
    offsets = FIGURE_6.STAGGER_STEP * np.arange(
        len(FIGURE_6.MASS_DATASETS) - 1, -1, -1
    )

    for (path, label, color), offset in zip(FIGURE_6.MASS_DATASETS, offsets):
        mass_da, intensity = FIGURE_6.load_distribution(path)
        axis.plot(
            mass_da / 1000.0,
            FIGURE_6.max_normalize_to_percent(intensity) + offset,
            color=color,
            linewidth=1.5,
            label=label,
        )

    axis.set_xlim(*FIGURE_6.MASS_LIMITS_KDA)
    axis.set_xticks(np.arange(0.0, 2000.1, 500.0))
    BASE_COMBINED.FIGURE_5.style_axis(axis, show_y_axis=False)
    add_panel_letter(axis, "H")


def make_figure() -> plt.Figure:
    """Build and return the complete eight-panel combined figure."""
    figure, axes = plt.subplots(
        nrows=4,
        ncols=2,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    for row, (tic_panel, mass_sample, letters) in enumerate(
        zip(
            BASE_COMBINED.FIGURE_4.PANELS,
            BASE_COMBINED.FIGURE_5.SAMPLES,
            PANEL_LETTERS,
        )
    ):
        _, tic_sample_label, tic_datasets = tic_panel
        mass_sample_label, zoom_limits, mass_datasets = mass_sample
        if tic_sample_label != mass_sample_label:
            raise ValueError(
                "Figure 4 and Figure 5 sample rows do not match: "
                f"{tic_sample_label!r} != {mass_sample_label!r}"
            )

        BASE_COMBINED.plot_tic_panel(
            axes[row, 0], letters[0], tic_sample_label, tic_datasets
        )
        BASE_COMBINED.plot_zoomed_mass_panel(
            axes[row, 1], letters[1], zoom_limits, mass_datasets
        )

    BASE_COMBINED.plot_tic_panel(
        axes[3, 0],
        "G",
        FIGURE_6_SAMPLE_LABEL,
        FIGURE_6.TIC_DATASETS,
    )
    plot_figure_6_mass_panel(axes[3, 1])

    axes[0, 0].legend(
        loc="upper right",
        frameon=False,
        fontsize=9.5,
        handlelength=2.0,
        borderaxespad=0.2,
    )
    axes[0, 1].legend(
        loc="upper right",
        frameon=False,
        fontsize=9.5,
        handlelength=2.0,
        borderaxespad=0.2,
    )

    axes[-1, 0].set_xlabel("Time (min)")
    axes[-1, 0].set_xticks(np.arange(0.0, 5.1, 1.0))
    axes[-1, 1].set_xlabel("Mass (kDa)")

    return figure


def main() -> None:
    """Render the eight-panel combined figure as a transparent PNG and PDF."""
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
