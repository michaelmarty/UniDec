"""Create the pre/post-deconvolution charge-versus-m/z plots for Figure S1."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np


DATA_FILE = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT"
    r"\Protein Mix 21 Injections Set 1\merged"
    r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_"
    r"set1_merged_21.npz"
)
UNIDEC_DIRECTORY = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT"
    r"\Protein Mix 21 Injections Set 1\merged"
    r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_"
    r"set1_merged_21_unidecfiles"
)
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S1"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S1.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S1.pdf"

UCCD_INPUT = UNIDEC_DIRECTORY / "uccd_input.bin"
UCCD_DECONVOLUTION = UNIDEC_DIRECTORY / "uccd_decon.bin"
CHROMATOGRAM_CONFIG = UNIDEC_DIRECTORY / "chroms.txt"

FIGURE_SIZE_INCHES = (6.5, 3.0)
OUTPUT_DPI = 300
CONTOUR_LEVELS = np.linspace(0.0, 100.0, 101)

# Make the repository copy of UniDec importable when this script is launched
# directly from the supplementary-figure directory.
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
UNIDEC_SOURCE = REPOSITORY_ROOT / "public" / "UniDec"
if str(UNIDEC_SOURCE) not in sys.path:
    sys.path.insert(0, str(UNIDEC_SOURCE))

from unidec import tools as ud  # noqa: E402
from unidec.modules.CDEng import read_uccd_binary  # noqa: E402
from unidec.modules.HTEng import UniChromCDEng  # noqa: E402


@dataclass(frozen=True)
class Swoop:
    """One saved UniChromCD swoop selection."""

    parameters: np.ndarray
    color: tuple[float, float, float, float]


def require_file(path: Path) -> None:
    """Raise a clear error when a required input is unavailable."""
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")


def load_unichrom_engine() -> UniChromCDEng:
    """Open the supplied NPZ with UniChromCD and load its saved configuration."""
    require_file(DATA_FILE)
    engine = UniChromCDEng()
    engine.open_file(str(DATA_FILE))

    loaded_directory = Path(engine.config.udir)
    if loaded_directory != UNIDEC_DIRECTORY:
        raise ValueError(
            "UniChromCD selected an unexpected result directory: "
            f"{loaded_directory} (expected {UNIDEC_DIRECTORY})"
        )
    return engine


def load_swoops(path: Path) -> list[Swoop]:
    """Read active swoop parameters and RGBA colors from UniChromCD chroms.txt."""
    require_file(path)
    swoops: list[Swoop] = []

    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split("\t")
        if len(fields) < 13:
            raise ValueError(f"Malformed chromatogram row {line_number}: {path}")

        parameters = np.asarray(fields[9:13], dtype=float)
        if parameters[0] < 0.0:
            continue

        color_values = np.fromstring(fields[1].strip("[]"), sep=" ", dtype=float)
        if color_values.size not in (3, 4):
            raise ValueError(f"Invalid swoop color on row {line_number}: {fields[1]}")
        if color_values.size == 3:
            color_values = np.append(color_values, 1.0)

        swoops.append(Swoop(parameters, tuple(color_values)))

    if not swoops:
        raise ValueError(f"No swoop selections found in: {path}")
    return swoops


def load_summed_histogram(
    path: Path,
) -> tuple[np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Load a sparse UniChromCD cube and sum its chromatographic dimension."""
    require_file(path)
    cube, axes = read_uccd_binary(path)
    histogram = np.sum(cube, axis=0, dtype=np.float64)
    del cube
    return np.clip(histogram, 0.0, None), axes


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
    """Normalize one histogram to the 0--100% intensity scale used by UniChromCD."""
    maximum = float(np.max(histogram))
    if not np.isfinite(maximum) or maximum <= 0.0:
        raise ValueError("Cannot normalize an empty histogram")
    return histogram * (100.0 / maximum)


def draw_swoops(
    axis: plt.Axes,
    swoops: list[Swoop],
    adduct_mass: float,
) -> None:
    """Overlay the saved UniChromCD swoop regions on one histogram."""
    for swoop in swoops:
        swoop_mz, _, upper_charge, lower_charge = ud.calc_swoop(
            swoop.parameters,
            adduct_mass=adduct_mass,
        )
        axis.fill_between(
            swoop_mz,
            lower_charge,
            upper_charge,
            facecolor=swoop.color,
            edgecolor=swoop.color,
            linewidth=0.65,
            alpha=0.42,
            zorder=4,
        )


def draw_histogram(
    axis: plt.Axes,
    mz_axis: np.ndarray,
    charge_axis: np.ndarray,
    histogram: np.ndarray,
    title: str,
    cmap: Colormap,
    norm: Normalize,
    swoops: list[Swoop],
    adduct_mass: float,
) -> None:
    """Draw one normalized charge-versus-m/z contour plot and its swoops."""
    axis.contourf(
        mz_axis,
        charge_axis,
        normalize_histogram(histogram),
        levels=CONTOUR_LEVELS,
        cmap=cmap,
        norm=norm,
        antialiased=False,
    )
    draw_swoops(axis, swoops, adduct_mass)

    axis.set_xlim(float(mz_axis[0]), float(mz_axis[-1]))
    axis.set_ylim(float(charge_axis[0]), float(charge_axis[-1]))
    axis.set_xlabel(r"$m/z$ (Th)")
    axis.set_title(title, pad=4.0)
    axis.tick_params(direction="out", length=2.5, width=0.8)
    for spine in axis.spines.values():
        spine.set_linewidth(0.8)


def make_figure() -> plt.Figure:
    """Build the side-by-side pre/post-deconvolution Figure S1."""
    engine = load_unichrom_engine()
    swoops = load_swoops(CHROMATOGRAM_CONFIG)

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
        swoops,
        float(engine.config.adductmass),
    )
    draw_histogram(
        axes[1],
        mz_axis,
        charge_axis,
        post_histogram,
        "Post-deconvolution",
        cmap,
        norm,
        swoops,
        float(engine.config.adductmass),
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
    """Render Figure S1 as transparent PNG and PDF files."""
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
