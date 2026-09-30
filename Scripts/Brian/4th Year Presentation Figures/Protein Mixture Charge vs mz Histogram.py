"""Deconvolve protein mixture datasets with UniChromCD and plot charge vs. m/z."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from unidec.modules.HTEng import UniChromCDEng


DATASETS = (
    (
        "1 INJ",
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 1\merged"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_"
            r"0o1mgmLs_set1_merged_1.npz"
        ),
    ),
    (
        "21 INJ",
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 1\merged"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_"
            r"0o1mgmLs_set1_merged_21.npz"
        ),
    ),
)
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Presentations\Brian\4th Year Presentation\Figures"
)
FIGURE_SIZE_INCHES = (8.0, 6.0)
OUTPUT_DPI = 300
CONTOUR_LEVELS = np.linspace(0.0, 100.0, 101)


def deconvolve_dataset(
    label: str,
    path: Path,
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """Run standard UniChromCD deconvolution on one summed dataset."""
    if not path.is_file():
        raise FileNotFoundError(f"CD-MS dataset not found: {path}")

    print(f"Running regular UniChromCD deconvolution for {label}: {path}")
    engine = UniChromCDEng()
    engine.open_file(str(path))
    engine.run_deconvolution(process_data=True)

    histogram = np.asarray(engine.harray, dtype=float)
    if histogram.ndim != 2 or not np.any(histogram > 0):
        raise ValueError(f"Deconvolution produced no charge/m/z histogram for {path}")
    if histogram.shape != (len(engine.ztab), len(engine.mz)):
        raise ValueError(
            f"Unexpected charge/m/z histogram shape for {path}: "
            f"{histogram.shape}"
        )

    return len(engine.farray), engine.mz.copy(), engine.ztab.copy(), histogram


def make_figure(
    result: tuple[str, int, np.ndarray, np.ndarray, np.ndarray],
) -> plt.Figure:
    """Build one deconvolved charge/m/z histogram plot."""
    label, ion_count, mz, charge, histogram = result
    figure, axis = plt.subplots(figsize=FIGURE_SIZE_INCHES, constrained_layout=True)
    colormap = plt.get_cmap("nipy_spectral")
    norm = Normalize(vmin=0.0, vmax=100.0)

    normalized = histogram / np.max(histogram) * 100.0
    contour = axis.contourf(
        mz,
        charge,
        normalized,
        levels=CONTOUR_LEVELS,
        cmap=colormap,
        norm=norm,
        antialiased=False,
    )
    axis.text(
        0.98,
        0.98,
        f"{label} ({ion_count:,} ions)",
        transform=axis.transAxes,
        color="white",
        fontsize=20,
        ha="right",
        va="top",
    )
    axis.set_xlabel("m/z (Th)")
    axis.set_ylabel("Charge")
    axis.tick_params(direction="out", width=2.0, length=6, labelsize=20)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    colorbar = figure.colorbar(contour, ax=axis, pad=0.03)
    colorbar.set_ticks([0.0, 50.0, 100.0], labels=["0", "%", "100"])
    colorbar.set_label("Relative intensity")
    return figure


def main() -> None:
    """Deconvolve both datasets and save separate transparent 300 DPI PNGs."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 20,
            "axes.labelsize": 24,
            "savefig.transparent": True,
        }
    )

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    for label, path in DATASETS:
        ion_count, mz, charge, histogram = deconvolve_dataset(label, path)
        result = (label, ion_count, mz, charge, histogram)
        output_png = OUTPUT_DIRECTORY / f"Protein Mixture Charge vs mz Histogram {label}.png"
        figure = make_figure(result)
        try:
            figure.savefig(output_png, dpi=OUTPUT_DPI, transparent=True)
        finally:
            plt.close(figure)
        print(f"Saved {output_png}")


if __name__ == "__main__":
    main()
