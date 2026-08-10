import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


def normalize_intensity(intensity, dataset_number):
    max_intensity = np.max(intensity)
    if max_intensity <= 0:
        raise ValueError(f"Cannot normalize dataset {dataset_number}; maximum intensity is {max_intensity}.")
    return intensity / max_intensity


def plot_mass_distributions(ax, files, legend_labels, normalize_data, stagger_plots, zoom_range):
    """Plot mass-distribution files on an existing Matplotlib axis."""
    if len(files) == 0:
        raise ValueError("files must contain at least one file path")
    if len(legend_labels) != len(files):
        raise ValueError("legend_labels must contain one label for each file")

    file_paths = [os.fspath(file) for file in files]
    data_list = [np.loadtxt(file_path) for file_path in file_paths]
    all_intensity = np.concatenate([data[:, 1] for data in data_list])
    stagger_offset_step = 0.1 if normalize_data else 0.1 * np.max(all_intensity)
    plotted_intensity_list = []

    for i, (data, legend_label) in enumerate(zip(data_list, legend_labels)):
        mass = data[:, 0] / 1000
        intensity = data[:, 1]

        if normalize_data:
            intensity = normalize_intensity(intensity, i + 1)
        if stagger_plots:
            intensity = intensity + i * stagger_offset_step

        plotted_intensity_list.append(intensity)
        ax.plot(mass, intensity, label=legend_label)

    ax.set_xlabel("Mass (kDa)")
    y_label = "Normalized Intensity" if normalize_data else "Intensity"
    if stagger_plots:
        y_label += " + Offset"
    ax.set_ylabel(y_label)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if stagger_plots:
        ax.set_ylim(0, np.max(np.concatenate(plotted_intensity_list)))
    elif normalize_data:
        ax.set_ylim(0, 1.0)
    else:
        ax.set_ylim(0, np.max(all_intensity))

    if zoom_range is None:
        all_mass = np.concatenate([data[:, 0] for data in data_list])
        ax.set_xlim(np.min(all_mass) / 1000, np.max(all_mass) / 1000)
    else:
        ax.set_xlim(*zoom_range)

    formatter = ScalarFormatter(useOffset=False, useMathText=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.legend()

    return ax


def plot_mass_distribution_grid(file_groups, legend_label_groups, nrows, ncols, normalize_data,
                                stagger_plots, zoom_range, figsize=None):
    """Plot groups of mass-distribution files in a row-major subplot grid.

    Each item in ``file_groups`` is the list of files to overlay in one panel.
    ``legend_label_groups`` must have the same nested structure.
    """
    if nrows <= 0 or ncols <= 0:
        raise ValueError("nrows and ncols must both be positive")
    if len(file_groups) == 0:
        raise ValueError("file_groups must contain at least one panel")
    if len(file_groups) != len(legend_label_groups):
        raise ValueError("legend_label_groups must contain one label list for each file group")

    panel_count = len(file_groups)
    grid_size = nrows * ncols
    if panel_count > grid_size:
        raise ValueError(f"The {nrows}x{ncols} grid has {grid_size} axes for {panel_count} panels")

    if figsize is None:
        figsize = (8 * ncols, 5 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    flat_axes = axes.ravel()

    for ax, files, legend_labels in zip(flat_axes, file_groups, legend_label_groups):
        plot_mass_distributions(
            ax,
            files,
            legend_labels=legend_labels,
            normalize_data=normalize_data,
            stagger_plots=stagger_plots,
            zoom_range=zoom_range,
        )

    for ax in flat_axes[panel_count:]:
        ax.set_visible(False)

    fig.tight_layout()
    return fig, axes


def plot_standalone_single_figure(files, legend_labels, normalize_data, stagger_plots, zoom_range, savepath=None, output_name=None):
    """Create, save, and display the standalone Figure 6 plot."""
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 28
    plt.rcParams["xtick.major.size"] = 12
    plt.rcParams["xtick.major.width"] = 2
    plt.rcParams["ytick.major.size"] = 12
    plt.rcParams["ytick.major.width"] = 2
    plt.rcParams["xtick.major.pad"] = 8
    plt.rcParams["ytick.major.pad"] = 8
    plt.rcParams["axes.labelpad"] = 18
    plt.rcParams["axes.titlepad"] = 18
    plt.rcParams["axes.linewidth"] = plt.rcParams["xtick.major.width"]
    plt.rcParams["pdf.fonttype"] = 42

    fig, ax = plt.subplots(figsize=(16, 9))
    plot_mass_distributions(
        ax,
        files,
        legend_labels=legend_labels,
        normalize_data=normalize_data,
        stagger_plots=stagger_plots,
        zoom_range=zoom_range,
    )
    fig.tight_layout()

    if output_name is None:
        output_name = "mass_distribution_comparison"
    if normalize_data:
        output_name += "_normalized"
    if stagger_plots:
        output_name += "_staggered"
    if zoom_range is not None:
        output_name += "_zoomed"

    if savepath is not None:
        output_directory = savepath
    else:
        output_directory = os.path.dirname(os.path.abspath(os.fspath(files[0])))
    output_path_png = os.path.join(output_directory, f"{output_name}.png")
    output_path_pdf = os.path.join(output_directory, f"{output_name}.pdf")
    fig.savefig(output_path_png, dpi=300, transparent=True)
    fig.savefig(output_path_pdf, dpi=300)
    plt.show()

    return output_path_png, output_path_pdf


if __name__ == "__main__":
    legend_labels = ("Replicate 1", "Replicate 2", "Replicate 3")
    normalize_data = True
    stagger_plots = False
    zoom_range = None


    topdir = r"Z:\Group Share\BHT\Temp\Paper"
    files = [
        "",
        "",
        "",
    ]
    plot_standalone_single_figure(
        files,
        legend_labels=legend_labels,
        normalize_data=normalize_data,
        stagger_plots=stagger_plots,
        zoom_range=zoom_range,
    )
