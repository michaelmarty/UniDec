"""Create Figure S25 from the Rabbit IgG demultiplexed UniDec stack.

The main figure contains charge-versus-m/z histograms and mass distributions
for five manually defined TIC windows.  The TIC selection plot, individual
peak plots, numeric results, and a compact archive of the plotted arrays are
also exported.  Pass ``--rerun-stack`` to repeat the GUI's "Run UniDec
Demult. Stack" workflow before rendering the figure.
"""

from __future__ import annotations

import argparse
import csv
import gc
import struct
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as scipy_filter

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
UNIDEC_PACKAGE_ROOT = REPOSITORY_ROOT / "public" / "UniDec"
sys.path.insert(0, str(UNIDEC_PACKAGE_ROOT))

from unidec.modules.HTEng import UniChromCDEng


DATA_FILE = Path(
    r"Z:\Group Share\BHT\Temp\Paper\DMT\High Throughput"
    r"\Amgen Rabbit IgG 21 Injections Set 1\merged"
    r"\20260605_BHT_haloc4column_amgenRabbitmAb_EPA13724_set1_"
    r"CDMS_0o01mgmL_merged_21.npz"
)
UNIDEC_DIRECTORY = DATA_FILE.with_name(f"{DATA_FILE.stem}_unidecfiles")
CONFIG_FILE = UNIDEC_DIRECTORY / "conf.dat"
UCCD_INPUT_FILE = UNIDEC_DIRECTORY / "uccd_input.bin"
UCCD_DECON_FILE = UNIDEC_DIRECTORY / "uccd_decon.bin"

OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S25"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S25.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S25.pdf"
TIC_PNG = OUTPUT_DIRECTORY / "Figure S25 TIC.png"
TIC_PDF = OUTPUT_DIRECTORY / "Figure S25 TIC.pdf"
RESULTS_CSV = OUTPUT_DIRECTORY / "Figure S25 results.csv"
RESULTS_NPZ = OUTPUT_DIRECTORY / "Figure S25 data.npz"

PEAK_WINDOWS = (
    ("Peak 1", 0.25, 0.50),
    ("Peak 2", 0.60, 1.36),
    ("Peak 3", 1.42, 1.84),
    ("Peak 4", 1.86, 2.24),
    ("Peak 5", 2.26, 2.44),
)
PEAK_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")

OUTPUT_DPI = 300
FIGURE_SIZE_INCHES = (7.0, 9.0)
INDIVIDUAL_FIGURE_SIZE_INCHES = (7.0, 3.0)
TIC_FIGURE_SIZE_INCHES = (7.0, 3.0)
MAX_MASS_PLOT_POINTS = 20_000
UCCD_MAGIC = b"UCCDBIN1"
UCCD_HEADER = struct.Struct("<8sIIIQ")
UCCD_RECORD = np.dtype([("index", "<u4"), ("intensity", "<f4")])
UCCD_CHUNK_RECORDS = 1_000_000


def validate_paths(require_stack: bool = True) -> None:
    """Validate the source, configuration, and optional stack artifacts."""
    required = [DATA_FILE, CONFIG_FILE]
    if require_stack:
        required.extend((UCCD_INPUT_FILE, UCCD_DECON_FILE))
    missing = [path for path in required if not path.is_file()]
    if missing:
        formatted = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(f"Required Figure S25 input files are missing:\n{formatted}")


def rerun_demultiplexed_stack() -> None:
    """Repeat the UniChromCD Run UniDec Demult. Stack button workflow."""
    validate_paths(require_stack=False)
    engine = UniChromCDEng()
    engine.open_file(str(DATA_FILE))

    print("Running full MRS demultiplexing stack")
    engine.run_all_ht()
    if engine.fullhstack_ht is None:
        raise RuntimeError("Full-stack demultiplexing did not produce an output stack")

    print("Running UCCD on the demultiplexed stack")
    deconvolved_stack = engine.decon_full_stack(
        hstack=engine.fullhstack_ht,
        chromaxis=engine.decontime,
    )
    if deconvolved_stack is None:
        raise RuntimeError("UCCD did not produce a deconvolved stack")

    del deconvolved_stack
    del engine
    gc.collect()


def read_sparse_stack(
    path: Path,
    collect_peak_histograms: bool,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], np.ndarray, list[np.ndarray]]:
    """Read a sparse UCCD file while accumulating its TIC and peak windows."""
    with path.open("rb") as infile:
        header = infile.read(UCCD_HEADER.size)
        if len(header) != UCCD_HEADER.size:
            raise ValueError(f"Incomplete UCCD header: {path}")
        magic, time_count, mz_count, charge_count, record_count = UCCD_HEADER.unpack(header)
        if magic != UCCD_MAGIC:
            raise ValueError(f"Invalid UCCD file signature: {path}")

        axis_count = time_count + mz_count + charge_count
        axes = np.fromfile(infile, dtype="<f4", count=axis_count)
        if len(axes) != axis_count:
            raise ValueError(f"Incomplete UCCD axes: {path}")
        time_axis = axes[:time_count]
        mz_axis = axes[time_count : time_count + mz_count]
        charge_axis = axes[time_count + mz_count :]

        tic = np.zeros(time_count, dtype=np.float64)
        peak_histograms: list[np.ndarray] = []
        peak_for_time = np.full(time_count, -1, dtype=np.int8)
        if collect_peak_histograms:
            peak_histograms = [
                np.zeros((charge_count, mz_count), dtype=np.float64)
                for _ in PEAK_WINDOWS
            ]
            for peak_index, (_, start, end) in enumerate(PEAK_WINDOWS):
                in_window = (time_axis >= start) & (time_axis <= end)
                if np.any(peak_for_time[in_window] >= 0):
                    raise ValueError("Figure S25 TIC windows must not overlap")
                peak_for_time[in_window] = peak_index

        records_left = record_count
        while records_left:
            requested = min(UCCD_CHUNK_RECORDS, records_left)
            records = np.fromfile(infile, dtype=UCCD_RECORD, count=requested)
            if len(records) != requested:
                raise ValueError(f"Incomplete sparse UCCD records: {path}")

            indexes = records["index"].astype(np.int64)
            maximum_index = time_count * mz_count * charge_count
            if np.any(indexes >= maximum_index):
                raise ValueError(f"Out-of-range sparse UCCD index: {path}")
            time_index, remainder = divmod(indexes, mz_count * charge_count)
            mz_index, charge_index = divmod(remainder, charge_count)
            intensity = records["intensity"].astype(np.float64)
            np.add.at(tic, time_index, intensity)

            if collect_peak_histograms:
                record_peak = peak_for_time[time_index]
                for peak_index, histogram in enumerate(peak_histograms):
                    selected = record_peak == peak_index
                    np.add.at(
                        histogram,
                        (charge_index[selected], mz_index[selected]),
                        intensity[selected],
                    )

            records_left -= requested

        if infile.read(1):
            raise ValueError(f"Unexpected trailing data in UCCD file: {path}")

    return (time_axis, mz_axis, charge_axis), tic, peak_histograms


def assert_matching_axes(
    first: tuple[np.ndarray, np.ndarray, np.ndarray],
    second: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
    """Ensure the MRS input and UCCD output describe the same cube."""
    for label, first_axis, second_axis in zip(
        ("time", "m/z", "charge"), first, second
    ):
        if first_axis.shape != second_axis.shape or not np.allclose(
            first_axis, second_axis
        ):
            raise ValueError(f"UCCD input and output {label} axes do not match")


def build_mass_transform_engine(
    mz_axis: np.ndarray,
    charge_axis: np.ndarray,
) -> UniChromCDEng:
    """Create a lightweight engine configured for the saved mass transform."""
    engine = UniChromCDEng()
    engine.config.config_import(str(CONFIG_FILE))
    engine.mz = np.asarray(mz_axis, dtype=float)
    engine.ztab = np.asarray(charge_axis, dtype=float)
    engine.X, engine.Y = np.meshgrid(engine.mz, engine.ztab, indexing="xy")
    engine.mass = (engine.X - engine.config.adductmass) * engine.Y
    engine.data.ztab = engine.ztab
    return engine


def calculate_mass_distributions(
    mz_axis: np.ndarray,
    charge_axis: np.ndarray,
    peak_histograms: list[np.ndarray],
) -> list[np.ndarray]:
    """Apply UniDec's configured m/z-to-mass transform to each TIC window."""
    engine = build_mass_transform_engine(mz_axis, charge_axis)
    distributions: list[np.ndarray] = []
    for (label, start, end), histogram in zip(PEAK_WINDOWS, peak_histograms):
        if not np.any(histogram > 0):
            raise ValueError(f"{label} ({start:g}-{end:g} min) contains no UCCD signal")
        engine.transform(harray=histogram)
        mass_data = np.asarray(engine.data.massdat, dtype=float).copy()
        mass_data[:, 1] = np.clip(mass_data[:, 1], 0.0, None)
        distributions.append(mass_data)
    return distributions


def minmax_downsample(
    x_values: np.ndarray,
    y_values: np.ndarray,
    max_points: int = MAX_MASS_PLOT_POINTS,
) -> tuple[np.ndarray, np.ndarray]:
    """Reduce dense lines while retaining the extrema in every display bin."""
    if len(x_values) <= max_points:
        return x_values, y_values

    block_size = int(np.ceil(len(x_values) / (max_points // 2)))
    indices: list[int] = [0]
    for start in range(0, len(x_values), block_size):
        stop = min(start + block_size, len(x_values))
        block = y_values[start:stop]
        indices.extend((start + int(np.argmin(block)), start + int(np.argmax(block))))
    indices.append(len(x_values) - 1)
    ordered_indices = np.unique(np.asarray(indices, dtype=int))
    return x_values[ordered_indices], y_values[ordered_indices]


def style_axis(axis: plt.Axes) -> None:
    """Apply the manuscript's common open-axis styling."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.0)
    axis.spines["bottom"].set_linewidth(1.0)
    axis.tick_params(direction="out", width=1.0, length=4)


def plot_charge_mz(
    axis: plt.Axes,
    figure: plt.Figure,
    mz_axis: np.ndarray,
    charge_axis: np.ndarray,
    histogram: np.ndarray,
) -> None:
    """Plot one charge-versus-m/z histogram with UniDec's contour interpolation."""
    maximum = float(np.max(histogram))
    if maximum <= 0.0:
        normalized = np.zeros_like(histogram)
    else:
        normalized = 100.0 * histogram / maximum

    # Match Plot2dBase.contourplot, which the UniChromCD GUI uses for the
    # supplied example. Sparse grids are expanded by one neighboring bin and
    # interpolated on a triangular mesh; denser grids use regular contours.
    plot_grid = normalized.transpose()
    active = plot_grid > 1.0
    if np.count_nonzero(active) / active.size < 0.1:
        active = scipy_filter.uniform_filter(active.astype(float), size=3) > 0
        active[0, 0] = True
        active[0, -1] = True
        active[-1, 0] = True
        active[-1, -1] = True
        mz_grid, charge_grid = np.meshgrid(
            mz_axis, charge_axis, indexing="ij"
        )
        image = axis.tricontourf(
            mz_grid[active],
            charge_grid[active],
            plot_grid[active],
            100,
            cmap="nipy_spectral",
            vmin=0.0,
            vmax=100.0,
        )
    else:
        image = axis.contourf(
            mz_axis,
            charge_axis,
            normalized,
            100,
            cmap="nipy_spectral",
            vmin=0.0,
            vmax=100.0,
        )
    axis.set_xlim(float(mz_axis[0]), float(mz_axis[-1]))
    axis.set_ylim(float(charge_axis[0]), float(charge_axis[-1]))
    colorbar = figure.colorbar(image, ax=axis, ticks=(0.0, 50.0, 100.0), pad=0.02)
    colorbar.ax.set_yticklabels(("0", "%", "100"))
    axis.set_ylabel("Charge")
    style_axis(axis)


def plot_mass_distribution(axis: plt.Axes, mass_data: np.ndarray) -> None:
    """Plot one downsampled zero-charge mass distribution in kDa."""
    mass_kda = mass_data[:, 0] / 1000.0
    mass_kda, intensity = minmax_downsample(mass_kda, mass_data[:, 1])
    axis.plot(mass_kda, intensity, color="black", linewidth=0.8)
    axis.set_xlim(float(mass_kda[0]), float(mass_kda[-1]))
    axis.set_ylabel("Intensity")
    style_axis(axis)


def make_main_figure(
    mz_axis: np.ndarray,
    charge_axis: np.ndarray,
    peak_histograms: list[np.ndarray],
    mass_distributions: list[np.ndarray],
) -> plt.Figure:
    """Build the five-row charge/m/z and mass-distribution figure."""
    figure, axes = plt.subplots(
        len(PEAK_WINDOWS),
        2,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )
    panel_index = 0
    for row, ((label, start, end), histogram, mass_data) in enumerate(
        zip(PEAK_WINDOWS, peak_histograms, mass_distributions)
    ):
        plot_charge_mz(axes[row, 0], figure, mz_axis, charge_axis, histogram)
        plot_mass_distribution(axes[row, 1], mass_data)
        axes[row, 0].set_title(
            f"{label}: {start:.2f}-{end:.2f} min",
            fontsize=10,
            fontweight="bold",
        )
        for axis in axes[row]:
            axis.text(
                0.015,
                0.96,
                chr(ord("A") + panel_index),
                transform=axis.transAxes,
                fontsize=12,
                fontweight="bold",
                ha="left",
                va="top",
                color="white" if axis is axes[row, 0] else "black",
                zorder=10,
            )
            panel_index += 1

    axes[-1, 0].set_xlabel(r"$m/z$ (Th)")
    axes[-1, 1].set_xlabel("Mass (kDa)")
    return figure


def make_tic_figure(time_axis: np.ndarray, tic: np.ndarray) -> plt.Figure:
    """Build the MRS TIC plot with all five selected windows highlighted."""
    figure, axis = plt.subplots(
        1,
        1,
        figsize=TIC_FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )
    axis.plot(time_axis, tic, color="gold", linewidth=1.5, label="DM_TIC")
    for (label, start, end), color in zip(PEAK_WINDOWS, PEAK_COLORS):
        axis.axvspan(start, end, color=color, alpha=0.22, label=label)
    axis.set_xlim(float(time_axis[0]), float(time_axis[-1]))
    axis.set_ylim(bottom=0.0)
    axis.set_xlabel("Time (min)")
    axis.set_ylabel("Ion Count")
    axis.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        frameon=False,
        ncol=6,
        fontsize=9,
        borderaxespad=0.0,
        columnspacing=1.2,
        handlelength=2.0,
    )
    style_axis(axis)
    return figure


def save_individual_peak_figures(
    mz_axis: np.ndarray,
    charge_axis: np.ndarray,
    peak_histograms: list[np.ndarray],
    mass_distributions: list[np.ndarray],
) -> None:
    """Save one paired PNG for each selected TIC peak."""
    for peak_number, ((label, start, end), histogram, mass_data) in enumerate(
        zip(PEAK_WINDOWS, peak_histograms, mass_distributions), start=1
    ):
        figure, axes = plt.subplots(
            1,
            2,
            figsize=INDIVIDUAL_FIGURE_SIZE_INCHES,
            constrained_layout=True,
        )
        plot_charge_mz(axes[0], figure, mz_axis, charge_axis, histogram)
        plot_mass_distribution(axes[1], mass_data)
        axes[0].set_xlabel(r"$m/z$ (Th)")
        axes[1].set_xlabel("Mass (kDa)")
        figure.suptitle(
            f"{label}: {start:.2f}-{end:.2f} min",
            fontsize=11,
            fontweight="bold",
        )
        output = OUTPUT_DIRECTORY / (
            f"Figure S25 Peak {peak_number} {start:.2f}-{end:.2f} min.png"
        )
        figure.savefig(output, dpi=OUTPUT_DPI, transparent=True)
        plt.close(figure)


def save_results(
    time_axis: np.ndarray,
    input_tic: np.ndarray,
    output_tic: np.ndarray,
    mz_axis: np.ndarray,
    charge_axis: np.ndarray,
    peak_histograms: list[np.ndarray],
    mass_distributions: list[np.ndarray],
) -> None:
    """Save interval summaries and all plotted numeric arrays."""
    rows: list[dict[str, float | int | str]] = []
    archive: dict[str, np.ndarray] = {
        "time_min": time_axis,
        "demultiplexed_tic": input_tic,
        "uccd_tic": output_tic,
        "mz_th": mz_axis,
        "charge": charge_axis,
    }
    for peak_number, ((label, start, end), histogram, mass_data) in enumerate(
        zip(PEAK_WINDOWS, peak_histograms, mass_distributions), start=1
    ):
        selected = (time_axis >= start) & (time_axis <= end)
        positive_mass = mass_data[:, 1] > 0
        rows.append(
            {
                "Peak": label,
                "Start Time (min)": start,
                "End Time (min)": end,
                "Time Points": int(np.count_nonzero(selected)),
                "Demultiplexed TIC Sum": float(np.sum(input_tic[selected])),
                "UCCD TIC Sum": float(np.sum(output_tic[selected])),
                "Charge-m/z Histogram Sum": float(np.sum(histogram)),
                "Mass Axis Minimum (Da)": float(mass_data[0, 0]),
                "Mass Axis Maximum (Da)": float(mass_data[-1, 0]),
                "Mass Distribution Maximum": float(np.max(mass_data[:, 1])),
                "Positive Mass Bins": int(np.count_nonzero(positive_mass)),
            }
        )
        key = f"peak_{peak_number}"
        archive[f"{key}_charge_mz_intensity"] = histogram
        archive[f"{key}_mass_da"] = mass_data[:, 0]
        archive[f"{key}_mass_intensity"] = mass_data[:, 1]

    with RESULTS_CSV.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(RESULTS_NPZ, **archive)


def parse_arguments() -> argparse.Namespace:
    """Parse the optional full-stack rerun switch."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rerun-stack",
        action="store_true",
        help='repeat the GUI "Run UniDec Demult. Stack" workflow first',
    )
    return parser.parse_args()


def main() -> None:
    """Load or rerun the stack, calculate all panels, and save Figure S25."""
    arguments = parse_arguments()
    if arguments.rerun_stack:
        rerun_demultiplexed_stack()
    validate_paths(require_stack=True)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10,
            "axes.labelsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    input_axes, input_tic, _ = read_sparse_stack(
        UCCD_INPUT_FILE, collect_peak_histograms=False
    )
    output_axes, output_tic, peak_histograms = read_sparse_stack(
        UCCD_DECON_FILE, collect_peak_histograms=True
    )
    assert_matching_axes(input_axes, output_axes)
    time_axis, mz_axis, charge_axis = output_axes
    mass_distributions = calculate_mass_distributions(
        mz_axis, charge_axis, peak_histograms
    )

    tic_figure = make_tic_figure(time_axis, input_tic)
    tic_figure.savefig(TIC_PNG, dpi=OUTPUT_DPI, transparent=True)
    tic_figure.savefig(TIC_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(tic_figure)

    figure = make_main_figure(
        mz_axis, charge_axis, peak_histograms, mass_distributions
    )
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)

    save_individual_peak_figures(
        mz_axis, charge_axis, peak_histograms, mass_distributions
    )
    save_results(
        time_axis,
        input_tic,
        output_tic,
        mz_axis,
        charge_axis,
        peak_histograms,
        mass_distributions,
    )

    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")
    print(f"Saved {TIC_PNG}")
    print(f"Saved {TIC_PDF}")
    print(f"Saved {RESULTS_CSV}")
    print(f"Saved {RESULTS_NPZ}")


if __name__ == "__main__":
    main()
