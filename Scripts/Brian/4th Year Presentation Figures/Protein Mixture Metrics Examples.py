"""Illustrate the S3 metrics on the 1 m/z-bin carbonic anhydrase trace from S5.

The spectrum is replicate 1, 21 injections, without frequency correction.
The displayed peak height is normalized to 100%; no smoothing is used.
Exported apex, centroid, and FWHM come from peakparam columns 1, 7, and 8.
The S5 spectrum and S3 peak parameters are different exports: the spectrum has
a 1 Da mass grid, whereas the saved apex is specified to 0.1 Da. Reference
markers show the saved masses. The FWHM arrow uses interpolated half-height
crossings (25.24 Da); its label reports the saved width (25.4 Da). The small
difference reflects the source exports' different mass sampling.
The shaded area illustrates the saved apex +/-30 Da integration window;
integrating the display-normalized trace would not reproduce the raw Area
column. S3 divides each raw area by the sum for the three matched proteins.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np


SOURCE_STEM = "20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21"
SOURCE_DIRECTORY = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT"
    r"\Protein Mix 21 Injections Set 1\merged"
) / f"{SOURCE_STEM}_unidecfiles"
MASS_DISTRIBUTION = SOURCE_DIRECTORY / f"{SOURCE_STEM}_mass.txt"
PEAK_PARAMETERS = (
    Path(r"Z:\Group Share\BHT\Temp\Paper\Stats\Protein Mix Rep 1")
    / f"{SOURCE_STEM}_peakparam.dat"
)

TARGET_MASS = 29024.63
MASS_LIMITS = (28960.0, 29080.0)
INTEGRATION_OFFSETS = (-30.0, 30.0)
FIGURE_SIZE_INCHES = (8, 6)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Presentations\Brian\4th Year Presentation\Figures"
)

TRACE_COLOR = "#0072B2"
CENTROID_COLOR = "#D55E00"
FWHM_COLOR = "#8E44AD"
APEX_COLOR = "#000000"
AREA_COLOR = "#009E73"
TRACE_WIDTH = 2.5
ANNOTATION_FONT_SIZE = 18


def load_example() -> tuple[np.ndarray, np.ndarray]:
    """Load the S5 trace and select the S3 carbonic anhydrase peak by apex mass."""
    data = np.loadtxt(MASS_DISTRIBUTION, dtype=float)
    if data.ndim != 2 or data.shape[1] != 2 or len(data) < 2:
        raise ValueError("Expected a two-column mass distribution")
    if not np.all(np.isfinite(data)) or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError("Mass distribution must be finite and ordered by mass")

    # S3's statistics copy has an uncommented header; also accept a # header.
    header = PEAK_PARAMETERS.read_text(encoding="utf-8").splitlines()[0]
    columns = header.lstrip("# ").split()
    if len(columns) < 8 or [columns[index] for index in (0, 5, 6, 7)] != [
        "Mass", "Area", "MassCentroid", "MassFWHM"
    ]:
        raise ValueError("Unexpected peak-parameter column order")
    parameters = np.loadtxt(PEAK_PARAMETERS, skiprows=1, ndmin=2)
    if parameters.shape[1] < 8 or not np.all(np.isfinite(parameters)):
        raise ValueError("Expected at least eight finite peak-parameter columns")
    peak = parameters[np.argmin(np.abs(parameters[:, 0] - TARGET_MASS))]
    if abs(peak[0] - TARGET_MASS) > 100.0 or peak[7] <= 0:
        raise ValueError("No valid carbonic anhydrase peak within 100 Da of the target")

    visible = data[
        (data[:, 0] >= MASS_LIMITS[0]) & (data[:, 0] <= MASS_LIMITS[1])
    ].copy()
    if len(visible) < 3 or np.max(visible[:, 1]) <= 0:
        raise ValueError("No positive carbonic anhydrase distribution in the plot window")
    visible[:, 1] *= 100.0 / np.max(visible[:, 1])
    return visible, peak


def half_height_bounds(data: np.ndarray) -> tuple[float, float]:
    """Interpolate the first half-height crossings on either side of the apex."""
    mass, intensity = data.T
    apex_index = int(np.argmax(intensity))
    left = apex_index
    right = apex_index
    while left > 0 and intensity[left] > 50.0:
        left -= 1
    while right < len(data) - 1 and intensity[right] > 50.0:
        right += 1
    if intensity[left] > 50.0 or intensity[right] > 50.0:
        raise ValueError("The plot window does not include both half-height crossings")
    lower = np.interp(50.0, intensity[left:left + 2], mass[left:left + 2])
    upper = np.interp(
        50.0, intensity[right - 1:right + 1][::-1], mass[right - 1:right + 1][::-1]
    )
    return float(lower), float(upper)


def make_base_figure(data: np.ndarray, title: str) -> tuple[plt.Figure, plt.Axes]:
    """Use the same mass window and intensity scale for all three examples."""
    figure, axis = plt.subplots(figsize=FIGURE_SIZE_INCHES, constrained_layout=True)
    axis.plot(data[:, 0], data[:, 1], color=TRACE_COLOR, linewidth=TRACE_WIDTH)
    axis.set_xlim(*MASS_LIMITS)
    axis.set_ylim(0.0, 150.0)
    axis.set_xticks([28970, 29020, 29070])
    axis.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    axis.set_yticks([0, 50, 100], labels=["0", "%", "100"])
    axis.set_xlabel("Mass (Da)")
    axis.set_title(title, fontsize=26, pad=12)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_bounds(0, 100)
    for side in ("left", "bottom"):
        axis.spines[side].set_linewidth(2.0)
    axis.tick_params(direction="out", width=2.0, length=6)
    return figure, axis


def make_centroid_figure(data: np.ndarray, peak: np.ndarray) -> plt.Figure:
    """Contrast the saved centroid within the FWHM window with the apex mass."""
    apex, centroid = peak[0], peak[6]
    figure, axis = make_base_figure(data, "Mass Centroid vs Peak Apex")
    for mass, color, linestyle in (
        (apex, APEX_COLOR, "--"), (centroid, CENTROID_COLOR, ":")
    ):
        axis.vlines(mass, 0, 104, color=color, linestyle=linestyle, linewidth=2)
    axis.text(
        0.02, 0.97, f"Apex: {apex:,.3f} Da",
        transform=axis.transAxes, fontsize=ANNOTATION_FONT_SIZE, va="top",
    )
    axis.text(
        0.02, 0.89, f"Centroid: {centroid:,.3f} Da",
        transform=axis.transAxes, color=CENTROID_COLOR,
        fontsize=ANNOTATION_FONT_SIZE, va="top",
    )

    # A mass ruler magnifies the saved values without inventing sub-bin data.
    inset = axis.inset_axes([0.64, 0.57, 0.33, 0.27])
    inset.set_title("Mass (zoom)", fontsize=12, pad=8)
    inset_padding = abs(centroid - apex) * 0.55
    inset.set_xlim(min(centroid, apex) - inset_padding, max(centroid, apex) + inset_padding)
    inset.set_ylim(0, 1)
    for mass, label, color, linestyle in (
        (centroid, "Centroid", CENTROID_COLOR, ":"),
        (apex, "Apex", APEX_COLOR, "--"),
    ):
        inset.vlines(mass, 0, 0.78, color=color, linestyle=linestyle, linewidth=2)
        inset.text(mass, 0.94, label, color=color, fontsize=11, ha="center", va="top")
    inset.annotate(
        "", xy=(apex, 0.45), xytext=(centroid, 0.45),
        arrowprops={"arrowstyle": "<->", "color": CENTROID_COLOR, "lw": 1.8},
    )
    inset.text(
        (apex + centroid) / 2, 0.25, f"{abs(apex - centroid):.3f} Da",
        ha="center", fontsize=12,
    )
    inset.set_xticks([centroid, apex], labels=[f"{centroid:.2f}", f"{apex:.2f}"])
    inset.tick_params(axis="x", labelsize=10)
    inset.set_yticks([])
    for side in ("top", "right", "left"):
        inset.spines[side].set_visible(False)
    return figure


def make_fwhm_figure(data: np.ndarray, peak: np.ndarray) -> plt.Figure:
    """Mark the half-height crossings and report peakparam column 8."""
    lower, upper = half_height_bounds(data)
    figure, axis = make_base_figure(data, "Full Width at Half Maximum")
    axis.axhline(50, color="#777777", linestyle=":", linewidth=1.5, zorder=0)
    axis.vlines([lower, upper], 0, 50, color=FWHM_COLOR, linestyle="--", linewidth=2)
    axis.annotate(
        "", xy=(upper, 50), xytext=(lower, 50),
        arrowprops={"arrowstyle": "<->", "color": FWHM_COLOR, "lw": 2.5},
    )
    axis.text(
        0.5, 0.94, f"FWHM = {peak[7]:.1f} Da", transform=axis.transAxes,
        fontsize=ANNOTATION_FONT_SIZE, color=FWHM_COLOR, ha="center", va="top",
    )
    axis.text(
        MASS_LIMITS[0] + 3, 54, "Half maximum",
        fontsize=ANNOTATION_FONT_SIZE, color="#555555", va="bottom",
    )
    return figure


def make_area_figure(data: np.ndarray, peak: np.ndarray) -> plt.Figure:
    """Shade the area inside the apex-relative integration bounds from S3."""
    lower, upper = peak[0] + np.asarray(INTEGRATION_OFFSETS)
    inside = (data[:, 0] > lower) & (data[:, 0] < upper)
    mass = np.concatenate(([lower], data[inside, 0], [upper]))
    intensity = np.interp(mass, data[:, 0], data[:, 1])
    figure, axis = make_base_figure(data, "Peak Area")
    axis.fill_between(mass, intensity, 0, color=AREA_COLOR, alpha=0.35)
    axis.vlines([lower, upper], 0, 112, color=AREA_COLOR, linestyle="--", linewidth=2)
    axis.annotate(
        "", xy=(upper, 112), xytext=(lower, 112),
        arrowprops={"arrowstyle": "<->", "color": AREA_COLOR, "lw": 2.5},
    )
    axis.text(
        0.5, 0.97, "Integration Window: Apex ± 30 Da",
        transform=axis.transAxes, ha="center", va="top", fontsize=ANNOTATION_FONT_SIZE,
    )
    axis.text(
        0.5, 0.88, f"{lower:,.1f}–{upper:,.1f} Da",
        transform=axis.transAxes, ha="center", va="top",
        fontsize=ANNOTATION_FONT_SIZE, color=AREA_COLOR,
    )
    return figure


if __name__ == "__main__":
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 22,
            "axes.labelsize": 28,
            "savefig.transparent": True,
        }
    )
    distribution, carbonic_anhydrase = load_example()
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    for name, plotter in (
        ("Mass Centroid", make_centroid_figure),
        ("FWHM", make_fwhm_figure),
        ("Peak Area", make_area_figure),
    ):
        figure = plotter(distribution, carbonic_anhydrase)
        output = OUTPUT_DIRECTORY / f"Protein Mixture Metrics - {name}.png"
        # Keep the requested canvas dimensions; do not use bbox_inches="tight".
        figure.savefig(output, dpi=OUTPUT_DPI, transparent=True)
        plt.close(figure)
        print(f"Saved {output}")
