"""Create the Human IgG pre/post-deconvolution plot for Figure S7."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, Normalize
import numpy as np


DATA_FILE = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\Amgen Human IgG 21 Injections Set 1\merged\20260331_BHT_agilentc8column_AmgenhumanmAb_PRO35304_CMDS_20msINJ_set1_0o1mgmLs_merged_21.npz"
)
UNIDEC_DIRECTORY = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\Amgen Human IgG 21 Injections Set 1\merged\20260331_BHT_agilentc8column_AmgenhumanmAb_PRO35304_CMDS_20msINJ_set1_0o1mgmLs_merged_21_unidecfiles"
)
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S7"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S7.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S7.pdf"

UCCD_INPUT = UNIDEC_DIRECTORY / "uccd_input.bin"
UCCD_DECONVOLUTION = UNIDEC_DIRECTORY / "uccd_decon.bin"

FIGURE_SIZE_INCHES = (6.5, 3.0)
OUTPUT_DPI = 300
CONTOUR_LEVELS = np.linspace(0.0, 100.0, 101)

# Make the repository copy of UniDec importable when this script is launched
# directly from the supplementary-figure directory.
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
UNIDEC_SOURCE = REPOSITORY_ROOT / "public" / "UniDec"
if str(UNIDEC_SOURCE) not in sys.path:
    sys.path.insert(0, str(UNIDEC_SOURCE))

from unidec.modules.CDEng import (  # noqa: E402
    UCCD_BINARY_HEADER,
    UCCD_BINARY_MAGIC,
    UCCD_BINARY_RECORD,
)
from unidec.modules.HTEng import UniChromCDEng  # noqa: E402


def require_file(path: Path) -> None:
    """Raise a clear error when a required input is unavailable."""
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")


def load_unichrom_engine() -> UniChromCDEng:
    """Open the supplied NPZ with UniChromCD and load its configuration."""
    require_file(DATA_FILE)
    engine = UniChromCDEng()
    # The plots use the saved configuration and UCCD cubes directly. Avoid the
    # optional filename-token parser because "rabbit" contains "bit".
    engine.open_cdms_file(str(DATA_FILE))

    loaded_directory = Path(engine.config.udir).resolve()
    expected_directory = UNIDEC_DIRECTORY.resolve()
    if loaded_directory != expected_directory:
        raise ValueError(
            "UniChromCD selected an unexpected result directory: "
            f"{loaded_directory} (expected {expected_directory})"
        )
    return engine


def load_summed_histogram(
    path: Path,
) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Stream a sparse UCCD cube and sum its chromatographic dimension."""
    require_file(path)
    with path.open("rb") as infile:
        header = infile.read(UCCD_BINARY_HEADER.size)
        if len(header) != UCCD_BINARY_HEADER.size:
            raise ValueError(f"Incomplete UCCD header in {path}")

        magic, nchrom, nmz, nz, nonzero_count = UCCD_BINARY_HEADER.unpack(header)
        if magic != UCCD_BINARY_MAGIC:
            raise ValueError(f"Invalid UCCD binary file: {path}")

        axis_count = nchrom + nmz + nz
        axis_byte_count = axis_count * np.dtype("<f4").itemsize
        axis_bytes = infile.read(axis_byte_count)
        if len(axis_bytes) != axis_byte_count:
            raise ValueError(f"Incomplete UCCD axes in {path}")

        axes = np.frombuffer(axis_bytes, dtype="<f4")
        chrom_axis = axes[:nchrom].copy()
        mz_axis = axes[nchrom : nchrom + nmz].copy()
        charge_axis = axes[nchrom + nmz :].copy()

        histogram = np.zeros((nz, nmz), dtype=np.float64)
        records_left = nonzero_count
        total_cube_size = nchrom * nmz * nz
        chunk_size = 1_000_000

        while records_left:
            count = min(chunk_size, records_left)
            record_byte_count = count * UCCD_BINARY_RECORD.itemsize
            record_bytes = infile.read(record_byte_count)
            if len(record_bytes) != record_byte_count:
                raise ValueError(f"Incomplete UCCD sparse data in {path}")

            records = np.frombuffer(record_bytes, dtype=UCCD_BINARY_RECORD)
            indexes = records["index"].astype(np.int64)
            if np.any(indexes >= total_cube_size):
                raise ValueError(f"Out-of-range sparse index in {path}")

            remainder = indexes % (nmz * nz)
            mz_indexes, charge_indexes = divmod(remainder, nz)
            flattened_indexes = mz_indexes * nz + charge_indexes
            chunk_sum = np.bincount(
                flattened_indexes,
                weights=records["intensity"],
                minlength=nmz * nz,
            )
            histogram += chunk_sum.reshape(nmz, nz).T
            records_left -= count

    return np.clip(histogram, 0.0, None), (
        chrom_axis,
        mz_axis,
        charge_axis,
    )


def validate_axes(
    first: tuple[np.ndarray, np.ndarray, np.ndarray],
    second: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
    """Ensure the pre- and post-deconvolution cubes use identical axes."""
    for name, first_axis, second_axis in zip(
        ("chromatography", "m/z", "charge"), first, second
    ):
        if first_axis.shape != second_axis.shape or not np.allclose(
            first_axis, second_axis
        ):
            raise ValueError(f"Pre/post {name} axes do not match")


def normalize_histogram(histogram: np.ndarray) -> np.ndarray:
    """Normalize one histogram to UniChromCD's 0--100% intensity scale."""
    maximum = float(np.max(histogram))
    if not np.isfinite(maximum) or maximum <= 0.0:
        raise ValueError("Cannot normalize an empty histogram")
    return histogram * (100.0 / maximum)


def draw_histogram(
    axis: plt.Axes,
    mz_axis: np.ndarray,
    charge_axis: np.ndarray,
    histogram: np.ndarray,
    title: str,
    cmap: Colormap,
    norm: Normalize,
) -> None:
    """Draw one normalized charge-versus-m/z contour plot."""
    axis.contourf(
        mz_axis,
        charge_axis,
        normalize_histogram(histogram),
        levels=CONTOUR_LEVELS,
        cmap=cmap,
        norm=norm,
        antialiased=False,
    )
    axis.set_xlim(float(mz_axis[0]), float(mz_axis[-1]))
    axis.set_ylim(float(charge_axis[0]), float(charge_axis[-1]))
    axis.set_xlabel(r"$m/z$ (Th)")
    axis.set_title(title, pad=4.0)
    axis.tick_params(direction="out", length=2.5, width=0.8)
    for spine in axis.spines.values():
        spine.set_linewidth(0.8)


def make_figure() -> plt.Figure:
    """Build the side-by-side pre/post-deconvolution Figure S7."""
    engine = load_unichrom_engine()
    pre_histogram, pre_axes = load_summed_histogram(UCCD_INPUT)
    post_histogram, post_axes = load_summed_histogram(UCCD_DECONVOLUTION)
    validate_axes(pre_axes, post_axes)
    _, mz_axis, charge_axis = pre_axes

    cmap = plt.get_cmap(engine.config.cmap)
    norm = Normalize(vmin=0.0, vmax=100.0)

    figure = plt.figure(figsize=FIGURE_SIZE_INCHES)
    grid = figure.add_gridspec(
        1,
        5,
        width_ratios=(1.0, 0.12, 1.0, 0.04, 0.04),
        left=0.075,
        right=0.945,
        bottom=0.17,
        top=0.90,
        wspace=0.0,
    )
    axes = (figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 2]))
    colorbar_axis = figure.add_subplot(grid[0, 4])

    draw_histogram(
        axes[0],
        mz_axis,
        charge_axis,
        pre_histogram,
        "Pre-deconvolution",
        cmap,
        norm,
    )
    draw_histogram(
        axes[1],
        mz_axis,
        charge_axis,
        post_histogram,
        "Post-deconvolution",
        cmap,
        norm,
    )
    axes[0].set_ylabel("Charge")
    axes[1].tick_params(labelleft=False)

    axes[0].text(
        -0.17,
        1.08,
        "A",
        transform=axes[0].transAxes,
        fontsize=14.0,
        fontweight="bold",
        ha="left",
        va="top",
    )
    axes[1].text(
        -0.08,
        1.08,
        "B",
        transform=axes[1].transAxes,
        fontsize=14.0,
        fontweight="bold",
        ha="left",
        va="top",
    )

    colorbar = figure.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        cax=colorbar_axis,
        ticks=(0.0, 50.0, 100.0),
    )
    colorbar.ax.set_yticklabels(("0", "%", "100"))
    colorbar.ax.tick_params(
        direction="out",
        length=2.5,
        width=0.8,
        labelsize=10.0,
        pad=2.0,
    )
    colorbar.outline.set_linewidth(0.8)
    return figure


def main() -> None:
    """Render Figure S7 as transparent PNG and PDF files."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10.0,
            "axes.labelsize": 10.0,
            "axes.titlesize": 10.0,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    figure = make_figure()
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()

