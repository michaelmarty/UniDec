import os
import sys
from pathlib import Path

import matplotlib
import numpy as np

UNIDEC_ROOT = Path(__file__).resolve().parents[3] / "public" / "UniDec"
sys.path.insert(0, str(UNIDEC_ROOT))

import unidec.engine as ud
import unidec.modules.MassSpecBuilder as msb

matplotlib.use(os.environ.get("MPLBACKEND", "TkAgg"))
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

matplotlib.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})

FIGURE_DIR = r"C:\Users\mm96978\The University of Texas at Austin\MartyLab - General\Papers\Artifact Suppression\Figures"

MASS_RANGE = (5000, 50000)
MZ_LIMITS = (1500, 3000)
COLORBAR_WIDTH_RATIO = 0.08
AXIS_LABEL_GAP_RATIO = 0.35
COLORBAR_TICK_GAP_RATIO = 0.28


def normalize(data):
    """Normalize one spectrum or deconvolution result to its own maximum."""
    maximum = np.max(data)
    if maximum == 0:
        return np.zeros_like(data)
    return data / maximum * 100


def run_deconvolution(spectrum, **suppression):
    engine = ud.UniDec()
    engine.pass_data_in(spectrum)
    engine.config.psig = 10
    engine.config.subbuff = 0
    engine.config.masslb, engine.config.massub = MASS_RANGE
    engine.config.suppression_satellite = suppression.get("satellite", 0)
    engine.config.suppression_harmonic = suppression.get("harmonic", 0)
    engine.config.suppression_topn = suppression.get("topn", 0)
    engine.config.suppression_topx = suppression.get("topx", 0)
    engine.config.beta = 0
    engine.autorun(auto_peak_width=True)

    mzgrid = engine.data.mzgrid
    mz = np.unique(mzgrid[:, 0])
    charge = np.unique(mzgrid[:, 1])
    intensity = mzgrid[:, 2].reshape(len(mz), len(charge))
    return mz, charge, normalize(intensity), engine.config.cmap


def label_panel(axis, panel, description, color="black", rotate_description=False,
                center_panel=False):
    axis.text(
        0.5 if center_panel else 0.03,
        0.95,
        panel,
        transform=axis.transAxes,
        va="top",
        ha="center" if center_panel else "left",
        color=color,
    )
    axis.text(
        0.5 if rotate_description else 0.95,
        0.88 if rotate_description else 0.95,
        description,
        transform=axis.transAxes,
        va="top",
        ha="center" if rotate_description else "right",
        rotation=90 if rotate_description else 0,
        color=color,
        fontsize=9,
    )


def label_primary_charge_states(axis, spectrum, params, charge_states):
    intensity = normalize(spectrum[:, 1])
    mass, _, charge_center, charge_width, _ = params[0]
    weights = np.exp(-0.5 * ((charge_states - charge_center) / charge_width) ** 2)
    primary_charges = charge_states[weights >= np.max(weights) * 0.1]

    for charge in primary_charges:
        mz = (mass + 1.00727647 * charge) / charge
        peak_index = np.argmin(np.abs(spectrum[:, 0] - mz))
        peak_intensity = intensity[peak_index]
        offset = -10 if peak_intensity > 90 else 3
        xoffset = -10 if peak_intensity > 90 else 0
        axis.annotate(
            f"{int(charge)}+",
            (spectrum[peak_index, 0], peak_intensity),
            xytext=(xoffset, offset),
            textcoords="offset points",
            ha="center",
            va="top" if offset < 0 else "bottom",
            fontsize=9,
        )


def plot_charge_vector(axis, charge, vector, panel, title, cmap, charge_limits):
    axis.imshow(
        vector[:, np.newaxis],
        origin="lower",
        aspect="auto",
        extent=(-0.5, 0.5, charge[0] - 0.5, charge[-1] + 0.5),
        cmap=cmap,
        vmin=0,
        vmax=100,
        interpolation="nearest",
    )
    axis.set_xlim(-0.5, 0.5)
    axis.set_ylim(charge_limits)
    axis.set_xticks([])
    axis.tick_params(axis="y", left=False, labelleft=False)
    label_panel(axis, panel, "", color="white", center_panel=True)
    axis.text(
        0.9,
        -0.02,
        title,
        transform=axis.transAxes,
        va="top",
        ha="right",
        rotation=45,
        color="black",
        fontsize=9,
        clip_on=False,
    )


def make_figure():
    # Use the same deterministic synthetic spectrum as Figure 2.
    np.random.seed(42)
    (spectrum, charge_states), params = msb.simple_spectrum2(
        [20000],
        resolution=50,
        psfun=0,
        baseline=0.03,
        noise=0.01,
        zwidth=0.75,
    )

    settings = [
        ("Normal", {}),
        ("Top 1", {"topn": 1}),
        ("Top 0.05", None),
        ("Harmonic", {"harmonic": 1}),
        ("Satellite 1", {"satellite": 1}),
    ]
    results = [
        run_deconvolution(spectrum, **options) if options is not None else None
        for _, options in settings
    ]

    mz, charge, normal_grid, cmap = results[0]
    selected_mz_index, selected_charge_index = np.unravel_index(
        np.argmax(normal_grid), normal_grid.shape
    )
    selected_mz = mz[selected_mz_index]

    # At the selected m/z, only these charges map into Figure 2's mass range.
    neutral_mz = selected_mz - 1.00727647
    min_charge = max(charge[0], np.ceil(MASS_RANGE[0] / neutral_mz))
    max_charge = min(charge[-1], np.floor(MASS_RANGE[1] / neutral_mz))
    charge_mask = (charge >= min_charge) & (charge <= max_charge)
    plotted_charge = charge[charge_mask]
    charge_limits = (min_charge - 0.5, max_charge + 0.5)

    fig = plt.figure(figsize=(6.5, 2.8))
    outer_grid = fig.add_gridspec(
        1,
        6,
        width_ratios=[
            1.0,
            AXIS_LABEL_GAP_RATIO,
            1.5,
            COLORBAR_WIDTH_RATIO,
            COLORBAR_TICK_GAP_RATIO,
            1.0,
        ],
        wspace=0.05,
    )
    vector_grid = outer_grid[0, 5].subgridspec(1, 5, wspace=0.08)
    raw_axis = fig.add_subplot(outer_grid[0, 0])
    grid_axis = fig.add_subplot(outer_grid[0, 2])
    colorbar_axis = fig.add_subplot(outer_grid[0, 3])
    vector_axes = [
        fig.add_subplot(vector_grid[0, i], sharey=grid_axis) for i in range(5)
    ]

    raw_axis.plot(spectrum[:, 0], normalize(spectrum[:, 1]), color="darkblue")
    raw_axis.set_xlim(MZ_LIMITS)
    raw_axis.set_xticks([1500, 2000, 2500, 3000])
    raw_axis.set_xlabel(r"$\mathit{m/z}$")
    raw_axis.set_ylim(0, 100)
    raw_axis.set_yticks([0, 50, 100], labels=["0", "%", "100"])
    raw_axis.spines["top"].set_visible(False)
    raw_axis.spines["right"].set_visible(False)
    label_panel(raw_axis, "A", "Data")
    label_primary_charge_states(raw_axis, spectrum, params, charge_states)

    contour = grid_axis.contourf(
        mz,
        charge,
        normal_grid.T,
        levels=np.linspace(0, 100, 101),
        cmap=cmap,
    )
    grid_axis.set_xlabel(r"$\mathit{m/z}$")
    grid_axis.set_ylabel("Charge")
    grid_axis.set_xlim(MZ_LIMITS)
    grid_axis.set_ylim(charge_limits)
    label_panel(grid_axis, "B", "Deconvolution", color="white")
    contour.set_rasterized(True)

    mz_step = np.median(np.diff(mz))
    selection_width = 40 * mz_step
    selection_box = Rectangle(
        (selected_mz - selection_width / 2, min_charge),
        selection_width,
        max_charge - min_charge,
        fill=False,
        edgecolor="red",
        linewidth=1.2,
        linestyle="--",
        zorder=10,
    )
    grid_axis.add_patch(selection_box)
    grid_axis.plot(
        selected_mz,
        charge[selected_charge_index],
        marker="o",
        markersize=3,
        color="red",
        zorder=11,
    )
    grid_axis.annotate(
        "Extracted\nColumn",
        (selected_mz + selection_width / 2, max_charge - 4),
        xytext=(4, 0),
        textcoords="offset points",
        color="red",
        fontsize=9,
        va="top",
        ha="left",
    )

    colorbar = fig.colorbar(contour, cax=colorbar_axis, ticks=[0, 50, 100])
    colorbar.set_ticklabels(["0", "%", "100"])

    normal_vector = normal_grid[selected_mz_index, charge_mask]
    for index, (axis, (title, _), result) in enumerate(
        zip(vector_axes, settings, results)
    ):
        if result is None:
            vector = normal_vector.copy()
            vector[vector < 0.05 * np.max(normal_vector)] = 0
            result_cmap = cmap
        else:
            result_mz, result_charge, result_grid, result_cmap = result
            if not np.array_equal(mz, result_mz) or not np.array_equal(charge, result_charge):
                raise ValueError("Suppression runs produced incompatible m/z-charge grids")
            vector = result_grid[selected_mz_index, charge_mask]
        plot_charge_vector(
            axis,
            plotted_charge,
            vector,
            chr(ord("C") + index),
            title,
            result_cmap,
            charge_limits,
        )
        if index == 1:
            axis.axhspan(
                9.5,
                charge_limits[1],
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

            axis.axhspan(
                8.5,
                charge_limits[0],
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

        if index == 2:
            axis.axhspan(
                10.5,
                charge_limits[1],
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

            axis.axhspan(
                7.5,
                charge_limits[0],
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

        if index == 3:
            axis.axhspan(
                17.5,
                18.5,
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

            axis.axhspan(
                20.5,
                19.5,
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

            axis.axhspan(
                15.5,
                16.5,
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

        if index == 4:
            axis.axhspan(
                17.5,
                16.5,
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

            axis.axhspan(
                19.5,
                18.5,
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

            axis.axhspan(
                9.5,
                10.5,
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

            axis.axhspan(
                7.5,
                8.5,
                facecolor="red",
                edgecolor="red",
                alpha=0.5,
                zorder=1,
            )

    # Match panel C's border to the extraction box in panel B.
    for spine in vector_axes[0].spines.values():
        spine.set_color("red")
        spine.set_linewidth(1.2)
        spine.set_linestyle("--")

    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.2, top=0.95)
    return fig


if __name__ == "__main__":
    figure = make_figure()
    figure.savefig(FIGURE_DIR + r"\Figure1.png", dpi=300)
    figure.savefig(FIGURE_DIR + r"\Figure1.pdf")
    plt.show()
