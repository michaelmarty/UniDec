"""Compare three UniDec mass-distribution regions for Figure S4."""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import StrMethodFormatter
import numpy as np


@dataclass(frozen=True)
class SpectrumSource:
    """Input path and plotting style for one mass distribution."""

    label: str
    path: Path
    color: str
    linestyle: str


SPECTRUM_SOURCES = (
    SpectrumSource(
        label=r"1 $m/z$ bin",
        path=Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 1\merged"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21_unidecfiles"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21_mass.txt"
        ),
        color="#0072B2",
        linestyle="-",
    ),
    SpectrumSource(
        label=r"0.1 $m/z$ bin",
        path=Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 1\merged"
            r"\Protien Mix 21 Injection 0.1 mz bin test_unidecfiles"
            r"\Protien Mix 21 Injection 0.1 mz bin test_mass.txt"
        ),
        color="#D55E00",
        linestyle="-",
    ),
    SpectrumSource(
        label=r"1 $m/z$ bin +FC",
        path=Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 1 Frequency Correction Test\merged"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21_unidecfiles"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21_mass.txt"
        ),
        color="#009E73",
        linestyle=(0, (5, 2)),
    ),
    SpectrumSource(
        label=r"0.1 $m/z$ bin +FC",
        path=Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 1 Frequency Correction Test\merged"
            r"\Protien Mix 21 Injection 0.1 mz bin test_unidecfiles"
            r"\Protien Mix 21 Injection 0.1 mz bin test_mass.txt"
        ),
        color="#CC79A7",
        linestyle=(0, (5, 2)),
    ),
)

# Theoretical protein masses and zoom limits, in Da.
PANELS = (
    ("A", "Cytochrome C", (12_320.0, 12_400.0), 12_360.08),
    ("B", "Myoglobin", (16_890.0, 17_010.0), 16_951.48),
    ("C", "Carbonic Anhydrase", (28_960.0, 29_080.0), 29_024.63),
)

FIGURE_SIZE_INCHES = (6.5, 7.0)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S4"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S4.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S4.pdf"
AXIS_STROKE_WIDTH = 0.8
TRACE_LINE_WIDTH = 1.15


def load_mass_distribution(path: Path) -> np.ndarray:
    """Load a two-column mass distribution and normalize it to 100%."""
    if not path.is_file():
        raise FileNotFoundError(f"Mass distribution not found: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two numeric columns in: {path}")

    data = data[:, :2]
    if not np.all(np.isfinite(data)):
        raise ValueError(f"Non-finite values found in: {path}")
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError(f"Mass values must be strictly increasing in: {path}")

    maximum = float(np.max(data[:, 1]))
    if maximum <= 0:
        raise ValueError(f"The intensity maximum must be positive in: {path}")
    data[:, 1] = data[:, 1] / maximum * 100.0
    return data


def style_axis(axis: plt.Axes, panel_letter: str) -> None:
    """Apply the manuscript styling shared by the SI figures."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(AXIS_STROKE_WIDTH)
    axis.spines["left"].set_bounds(0.0, 100.0)
    axis.spines["bottom"].set_linewidth(AXIS_STROKE_WIDTH)
    axis.tick_params(direction="out", width=AXIS_STROKE_WIDTH, length=4)
    axis.text(
        0.012,
        0.97,
        panel_letter,
        transform=axis.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def plot_trace(
    axis: plt.Axes,
    data: np.ndarray,
    source: SpectrumSource,
    limits: tuple[float, float],
    *,
    rasterized: bool,
) -> None:
    """Plot only the portion of a spectrum visible within an axis."""
    first = int(np.searchsorted(data[:, 0], limits[0], side="left"))
    last = int(np.searchsorted(data[:, 0], limits[1], side="right"))
    visible = data[first:last]
    if visible.size == 0:
        raise ValueError(f"No data from {source.path} falls within {limits}")

    axis.plot(
        visible[:, 0],
        visible[:, 1],
        color=source.color,
        linestyle=source.linestyle,
        linewidth=TRACE_LINE_WIDTH,
        label=source.label,
        rasterized=rasterized,
        zorder=2,
    )


def make_figure(distributions: tuple[np.ndarray, ...]) -> plt.Figure:
    """Build and return the three protein-region plots."""
    figure, axes = plt.subplots(
        nrows=len(PANELS),
        ncols=1,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    for axis, (panel_letter, title, limits, theoretical_mass) in zip(axes, PANELS):
        for data, source in zip(distributions, SPECTRUM_SOURCES):
            plot_trace(
                axis,
                data,
                source,
                limits,
                rasterized=False,
            )

        if theoretical_mass is not None:
            axis.axvline(
                theoretical_mass,
                color="#000000",
                linestyle=(0, (3, 2)),
                linewidth=1.0,
                zorder=3,
            )

        axis.set_xlim(*limits)
        axis.set_ylim(0.0, 105.0)
        display_title = (
            title
            if theoretical_mass is None
            else f"{title} ({theoretical_mass:,.2f} Da)"
        )
        axis.set_title(display_title, fontsize=12.5, fontweight="bold", pad=6)
        axis.set_xlabel("Mass (Da)")
        axis.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
        axis.set_yticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
        style_axis(axis, panel_letter)

    trace_handles = [
        Line2D(
            [0],
            [0],
            color=source.color,
            linestyle=source.linestyle,
            linewidth=TRACE_LINE_WIDTH,
            label=source.label,
        )
        for source in SPECTRUM_SOURCES
    ]
    theoretical_legend = axes[0].legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="#000000",
                linestyle=(0, (3, 2)),
                linewidth=1.0,
                label="Theo. mass",
            )
        ],
        loc="upper right",
        frameon=False,
        fontsize=10,
        handlelength=2.5,
    )
    axes[0].add_artist(theoretical_legend)
    axes[0].legend(
        handles=trace_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.20),
        ncols=2,
        frameon=False,
        fontsize=10,
        columnspacing=1.3,
        handlelength=3.0,
        borderaxespad=0.0,
    )
    return figure


def main() -> None:
    """Load all four spectra and save Figure S4 as PNG and PDF."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10.5,
            "axes.labelsize": 11.5,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    distributions = tuple(
        load_mass_distribution(source.path) for source in SPECTRUM_SOURCES
    )
    figure = make_figure(distributions)
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
