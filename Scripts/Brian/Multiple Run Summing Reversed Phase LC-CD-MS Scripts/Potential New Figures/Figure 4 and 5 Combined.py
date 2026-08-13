"""Combine the Figure 4 TICs with the zoomed Figure 5 mass distributions.

The three rows contain the Waters, human, and rabbit samples, respectively.
TIC panels are in the left column and zoomed mass-distribution panels are in
the right column, giving the row-major panel order A--F.
"""

from importlib import util
from pathlib import Path
from types import ModuleType

import matplotlib.pyplot as plt
import numpy as np


FIGURE_SIZE_INCHES = (7.0, 9.0)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure 4 and 5 Combined.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure 4 and 5 Combined.pdf"

SOURCE_DIRECTORY = OUTPUT_DIRECTORY.parent
PANEL_LETTERS = (("A", "B"), ("C", "D"), ("E", "F"))


def load_figure_module(filename: str, module_name: str) -> ModuleType:
    """Load one of the source figure scripts from the parent directory."""
    source_path = SOURCE_DIRECTORY / filename
    specification = util.spec_from_file_location(module_name, source_path)
    if specification is None or specification.loader is None:
        raise ImportError(f"Could not load figure source: {source_path}")

    module = util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


FIGURE_4 = load_figure_module("Figure 4.py", "figure_4_source")
FIGURE_5 = load_figure_module("Figure 5.py", "figure_5_source")


def plot_tic_panel(
    axis: plt.Axes,
    panel_letter: str,
    sample_label: str,
    datasets: list[tuple[Path, str, str]],
) -> None:
    """Plot one staggered TIC panel using the Figure 4 styling."""
    offsets = FIGURE_4.STAGGER_STEP * np.arange(len(datasets) - 1, -1, -1)
    retention_time = None

    for (path, label, color), offset in zip(datasets, offsets):
        time, intensity = FIGURE_4.load_chromatogram(path)
        if label == "21-injection CD-MS":
            retention_time = float(time[np.argmax(intensity)])
        axis.plot(
            time,
            FIGURE_4.normalize_to_percent(intensity) + offset,
            color=color,
            linewidth=2.0,
            label=label,
            zorder=2,
        )

    FIGURE_4.style_axis(axis, panel_letter, sample_label)
    if retention_time is not None:
        axis.axvline(
            retention_time,
            ymin=0.0,
            ymax=0.72,
            color="#000000",
            linestyle=(0, (4, 3)),
            linewidth=1.5,
            zorder=3,
        )
        axis.text(
            retention_time - 0.08,
            0.73,
            f"{retention_time:.2f} min",
            transform=axis.get_xaxis_transform(),
            color="#000000",
            fontsize=7.0,
            ha="right",
            va="bottom",
            zorder=4,
        )


def plot_zoomed_mass_panel(
    axis: plt.Axes,
    panel_letter: str,
    zoom_limits: tuple[float, float],
    datasets: list[tuple[Path, str, str]],
) -> None:
    """Plot one staggered zoom panel using the Figure 5 styling."""
    offsets = FIGURE_5.STAGGER_STEP * np.arange(len(datasets) - 1, -1, -1)

    for (path, label, color), offset in zip(datasets, offsets):
        mass_da, intensity = FIGURE_5.load_mass_distribution(path)
        mass_kda = mass_da / 1000.0
        normalized = FIGURE_5.normalize_to_percent(intensity) + offset
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
    FIGURE_5.style_axis(axis, show_y_axis=False)
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


def make_figure() -> plt.Figure:
    """Build and return the six-panel combined figure."""
    figure, axes = plt.subplots(
        nrows=3,
        ncols=2,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    for row, (tic_panel, mass_sample, letters) in enumerate(
        zip(FIGURE_4.PANELS, FIGURE_5.SAMPLES, PANEL_LETTERS)
    ):
        _, tic_sample_label, tic_datasets = tic_panel
        mass_sample_label, zoom_limits, mass_datasets = mass_sample
        if tic_sample_label != mass_sample_label:
            raise ValueError(
                "Figure 4 and Figure 5 sample rows do not match: "
                f"{tic_sample_label!r} != {mass_sample_label!r}"
            )

        plot_tic_panel(axes[row, 0], letters[0], tic_sample_label, tic_datasets)
        plot_zoomed_mass_panel(
            axes[row, 1], letters[1], zoom_limits, mass_datasets
        )

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
    """Render the combined figure as a transparent PNG and PDF."""
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
