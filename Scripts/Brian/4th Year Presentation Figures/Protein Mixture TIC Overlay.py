"""Create the normalized protein mixture TIC overlay copied from Figure 2B."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATA_ROOT = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\Protein Mix 21 Injections Set 1"
)
MERGED_ROOT = DATA_ROOT / "merged"
TWENTY_ONE_INJECTION_DIR = MERGED_ROOT / (
    "20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_21_"
    "unidecfiles"
)
ONE_INJECTION_DIR = MERGED_ROOT / (
    "20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_set1_merged_1_"
    "unidecfiles"
)
ENSEMBLE_CHROMATOGRAM = Path(
    r"Z:\Group Share\BHT\Temp\Paper\Ensemble Protein Mix Triplicate"
) / "UniDec_Figures_and_Files" / "Protein mix Rep 1 chrom.txt"

TIC_21 = TWENTY_ONE_INJECTION_DIR / "TIC_decon.txt"
TIC_1 = ONE_INJECTION_DIR / "TIC_chrom.txt"

TIME_LIMITS = (0.0, 5.5)
FIGURE_SIZE_INCHES = (8, 6)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Presentations\Brian\4th Year Presentation\Figures"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Protein Mixture TIC Overlay.png"

# Normalize each trace independently, then offset by 20, 10, and 0 percent.
STAGGER_STEP = 10.0


def load_chromatogram(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load and validate a two-column time/intensity chromatogram."""
    if not path.is_file():
        raise FileNotFoundError(f"Chromatogram not found: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in: {path}")

    time = data[:, 0]
    intensity = data[:, 1]
    finite = np.isfinite(time) & np.isfinite(intensity)
    if not np.any(finite):
        raise ValueError(f"No finite time/intensity pairs in: {path}")

    # Deconvolution can leave numerical noise around -1e-14. Ion counts are
    # physically non-negative, so remove only those negative values here.
    return time[finite], np.clip(intensity[finite], 0.0, None)


def normalize_to_percent(intensity: np.ndarray) -> np.ndarray:
    """Min-max normalize a signal to exactly 0--100%."""
    minimum = float(np.min(intensity))
    maximum = float(np.max(intensity))
    span = maximum - minimum
    if span <= 0:
        return np.zeros_like(intensity)
    return 100.0 * (intensity - minimum) / span


def style_axis(axis: plt.Axes) -> None:
    """Apply presentation styling without a panel letter."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(2.0)
    axis.spines["bottom"].set_linewidth(2.0)
    axis.tick_params(direction="out", width=2.0, length=6)
    axis.set_xlim(*TIME_LIMITS)
    axis.set_xlabel("Time (min)")


def label_protein_peaks(axis: plt.Axes) -> None:
    """Identify TIC peaks using the protein-specific EICs from Figure 2A."""
    proteins = (
        ("Cytochrome C", "Swoop_mz__1150_z__11_decon.txt", (-30, 30), "center", "bottom"),
        ("Myoglobin", "Swoop_mz__1260_z__14_decon.txt", (-40, 35), "center", "bottom"),
        ("Carbonic\nAnhydrase", "Swoop_mz__1387_z__21_decon.txt", (32, 0), "left", "top"),
    )
    for name, filename, offset, horizontal, vertical in proteins:
        time, intensity = load_chromatogram(TWENTY_ONE_INJECTION_DIR / filename)
        retention_time = time[np.argmax(intensity)]
        candidates = []
        for trace in axis.lines:
            trace_time, trace_intensity = trace.get_data()
            nearby = np.flatnonzero(np.abs(trace_time - retention_time) <= 0.05)
            if nearby.size:
                index = nearby[np.argmax(trace_intensity[nearby])]
                candidates.append((trace_time[index], trace_intensity[index]))
        if not candidates:
            raise ValueError(f"No TIC peak near the {name} EIC retention time")
        peak_time, peak_intensity = max(candidates, key=lambda point: point[1])
        axis.annotate(
            name,
            xy=(peak_time, peak_intensity),
            xytext=offset,
            textcoords="offset points",
            ha=horizontal,
            va=vertical,
            fontsize=18,
            arrowprops={"arrowstyle": "->", "color": "#555555", "lw": 1.5},
        )


def make_figure() -> plt.Figure:
    """Compare Ensemble MS with one- and twenty-one-injection CD-MS."""
    comparison = [
        (ENSEMBLE_CHROMATOGRAM, "Ensemble MS", "#333333"),
        (TIC_1, "1 INJ CD-MS", "#0072B2"),
        (TIC_21, "21 INJ CD-MS", "#D55E00"),
    ]
    figure, axis = plt.subplots(
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )
    offsets = STAGGER_STEP * np.arange(len(comparison) - 1, -1, -1)
    for (path, label, color), offset in zip(comparison, offsets):
        time, intensity = load_chromatogram(path)
        axis.plot(
            time,
            normalize_to_percent(intensity) + offset,
            color=color,
            linewidth=2.5,
            label=label,
        )

    style_axis(axis)
    axis.set_ylim(0.0, 124.0)
    axis.set_ylabel("")
    axis.set_yticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
    axis.spines["left"].set_bounds(0.0, 100.0)
    axis.legend(
        loc="upper left",
        frameon=False,
        handlelength=1.5,
    )
    label_protein_peaks(axis)
    return figure


if __name__ == "__main__":
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
    # Preserve the exact 8 by 6 inch canvas without bbox_inches="tight".
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
