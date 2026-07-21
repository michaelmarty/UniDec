"""Interactive plotter for UniDec two-column text exports.

The script combines the original m/z spectra, mass distribution, and TIC
plotting workflows into one menu-driven tool.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


FIGURE_SIZE = (16, 9)
FONT_FAMILY = "Arial"
FONT_SIZE = 28
TICK_SIZE = 12
TICK_WIDTH = 2
STAGGER_STEP_NORMALIZED = 0.1
STAGGER_STEP_FRACTION = 0.1


@dataclass(frozen=True)
class PlotConfig:
    """Settings that differ between the three supported plot types."""

    menu_key: str
    name: str
    x_label: str
    output_stem: str
    single_output_stem: str
    zoom_prompt: str
    zoom_start_prompt: str
    zoom_end_prompt: str
    zoom_value_error: str
    zoom_order_error: str
    x_scale: float = 1.0


@dataclass(frozen=True)
class Dataset:
    path: Path
    values: np.ndarray
    legend_label: Optional[str] = None


@dataclass(frozen=True)
class PlotOptions:
    normalize: bool
    stagger: bool
    use_zoom: bool
    zoom_start: Optional[float] = None
    zoom_end: Optional[float] = None


# Only the axis text, mass scaling, zoom prompts, and output filenames differ.
PLOT_CONFIGS = (
    PlotConfig(
        menu_key="1",
        name="m/z Spectra",
        x_label=r"$m/z$",
        output_stem="mz_spectra_comparison",
        single_output_stem="mass_spectrum",
        zoom_prompt="Do you want to use a zoomed m/z x-axis range? (y/n): ",
        zoom_start_prompt="Enter the start m/z for zoomed x-axis: ",
        zoom_end_prompt="Enter the end m/z for zoomed x-axis: ",
        zoom_value_error="Please enter numeric values for the start and end m/z values.",
        zoom_order_error="End m/z must be greater than start m/z.",
    ),
    PlotConfig(
        menu_key="2",
        name="Mass Distribution",
        x_label=r"Mass (kDa)",
        output_stem="mass_distribution_comparison",
        single_output_stem="mass_distribution",
        zoom_prompt="Do you want to use a zoomed mass x-axis range? (y/n): ",
        zoom_start_prompt="Enter the start mass for zoomed x-axis (in kDa): ",
        zoom_end_prompt="Enter the end mass for zoomed x-axis (in kDa): ",
        zoom_value_error="Please enter numeric values for the start and end masses.",
        zoom_order_error="End mass must be greater than start mass.",
        x_scale=1 / 1000,
    ),
    PlotConfig(
        menu_key="3",
        name="Total Ion Chromatogram",
        x_label=r"Time (minutes)",
        output_stem="total_ion_chromatogram_comparison",
        single_output_stem="total_ion_chromatogram",
        zoom_prompt="Do you want to use a zoomed x-axis range? (y/n): ",
        zoom_start_prompt="Enter the start time for zoomed x-axis (in minutes): ",
        zoom_end_prompt="Enter the end time for zoomed x-axis (in minutes): ",
        zoom_value_error="Please enter numeric values for the start and end times.",
        zoom_order_error="End time must be greater than start time.",
    ),
)
PLOT_CONFIG_BY_KEY = {config.menu_key: config for config in PLOT_CONFIGS}


def prompt_yes_no(prompt):
    while True:
        response = input(prompt).strip().lower()
        if response in ("y", "yes"):
            return True
        if response in ("n", "no"):
            return False
        print("Please enter y or n.")


def prompt_positive_integer(prompt):
    while True:
        try:
            value = int(input(prompt))
        except ValueError:
            print("Please enter a whole number.")
            continue

        if value <= 0:
            print("Please enter a number greater than 0.")
            continue

        return value


def prompt_plot_config():
    print("Select the plot type:")
    for config in PLOT_CONFIGS:
        print(f"{config.menu_key}. {config.name}")

    valid_choices = ", ".join(config.menu_key for config in PLOT_CONFIGS)
    while True:
        response = input(f"Enter {valid_choices}: ").strip()
        if response in PLOT_CONFIG_BY_KEY:
            return PLOT_CONFIG_BY_KEY[response]
        print(f"Please enter {valid_choices}.")


def prompt_zoom_range(config):
    while True:
        try:
            zoom_start = float(input(config.zoom_start_prompt))
            zoom_end = float(input(config.zoom_end_prompt))
        except ValueError:
            print(config.zoom_value_error)
            continue

        if zoom_end <= zoom_start:
            print(config.zoom_order_error)
            continue

        return zoom_start, zoom_end


def prompt_plot_options(config, single_dataset):
    normalize = prompt_yes_no("Do you want to normalize the data? (y/n): ")
    stagger = False
    if not single_dataset:
        stagger = prompt_yes_no("Do you want to stagger the plots? (y/n): ")

    use_zoom = prompt_yes_no(config.zoom_prompt)

    if not use_zoom:
        return PlotOptions(normalize=normalize, stagger=stagger, use_zoom=False)

    zoom_start, zoom_end = prompt_zoom_range(config)
    return PlotOptions(
        normalize=normalize,
        stagger=stagger,
        use_zoom=True,
        zoom_start=zoom_start,
        zoom_end=zoom_end,
    )


def clean_file_path(raw_path):
    return Path(raw_path.strip().strip('"').strip("'"))


def load_dataset(path, dataset_number):
    """Load a UniDec-style text export with x values in column 1 and intensity in column 2."""
    data = np.loadtxt(path)

    if data.ndim == 1:
        if data.size == 2:
            data = data.reshape(1, 2)
        else:
            raise ValueError(
                f"Dataset {dataset_number} must have two columns: x-axis values and intensity."
            )

    if data.ndim != 2:
        raise ValueError(
            f"Dataset {dataset_number} must have two columns: x-axis values and intensity."
        )

    if data.shape[1] < 2:
        raise ValueError(
            f"Dataset {dataset_number} must have at least two columns: x-axis values and intensity."
        )

    return data


def prompt_datasets():
    num_datasets = prompt_positive_integer("Enter the number of datasets: ")
    datasets = []

    for index in range(1, num_datasets + 1):
        path = clean_file_path(input(f"Enter the file path for dataset {index} (quotes optional): "))
        values = load_dataset(path, index)
        legend_label = None
        if num_datasets > 1:
            legend_label = input(f"Enter the legend label for dataset {index}: ")
        datasets.append(Dataset(path=path, values=values, legend_label=legend_label))

    return datasets


def configure_matplotlib():
    """Apply the shared visual style from the original comparison scripts."""
    plt.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": FONT_SIZE,
            "xtick.major.size": TICK_SIZE,
            "xtick.major.width": TICK_WIDTH,
            "ytick.major.size": TICK_SIZE,
            "ytick.major.width": TICK_WIDTH,
            "xtick.major.pad": 8,
            "ytick.major.pad": 8,
            "axes.labelpad": 18,
            "axes.titlepad": 18,
            "axes.linewidth": TICK_WIDTH,
            "pdf.fonttype": 42,
        }
    )


def get_x_values(dataset, config):
    # Mass distribution exports are in Da; plot them in kDa.
    return dataset.values[:, 0] * config.x_scale


def get_intensity_values(dataset):
    return dataset.values[:, 1].astype(float, copy=True)


def normalize_intensity(intensity, dataset_number):
    max_intensity = np.max(intensity)
    if max_intensity <= 0:
        raise ValueError(f"Cannot normalize dataset {dataset_number}; maximum intensity is {max_intensity}.")
    return intensity / max_intensity


def get_stagger_offset_step(all_intensities, normalize):
    if normalize:
        return STAGGER_STEP_NORMALIZED
    return STAGGER_STEP_FRACTION * np.max(all_intensities)


def is_single_dataset(datasets):
    return len(datasets) == 1


def plot_datasets(config, datasets, options):
    all_raw_intensities = np.concatenate([get_intensity_values(dataset) for dataset in datasets])
    stagger_offset_step = get_stagger_offset_step(all_raw_intensities, options.normalize)
    plotted_intensities = []
    single_dataset = is_single_dataset(datasets)

    for dataset_number, dataset in enumerate(datasets, start=1):
        intensity = get_intensity_values(dataset)

        if options.normalize:
            intensity = normalize_intensity(intensity, dataset_number)
        if options.stagger:
            intensity = intensity + ((dataset_number - 1) * stagger_offset_step)

        plotted_intensities.append(intensity)
        if single_dataset:
            plt.plot(get_x_values(dataset, config), intensity, color="black")
        else:
            plt.plot(get_x_values(dataset, config), intensity, label=dataset.legend_label)

    return all_raw_intensities, plotted_intensities


def format_axes(ax, config, datasets, options, all_raw_intensities, plotted_intensities):
    ax.set_xlabel(config.x_label)

    y_label = "Normalized Intensity" if options.normalize else "Intensity"
    if options.stagger:
        y_label += " + Offset"
    ax.set_ylabel(y_label)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    set_y_limits(ax, options, all_raw_intensities, plotted_intensities)
    set_x_limits(ax, config, datasets, options)
    format_x_axis_without_scientific_notation(ax)


def set_y_limits(ax, options, all_raw_intensities, plotted_intensities):
    if options.stagger:
        y_max = np.max(np.concatenate(plotted_intensities))
    elif options.normalize:
        y_max = 1.0
    else:
        y_max = np.max(all_raw_intensities)

    ax.set_ylim(0, y_max)


def set_x_limits(ax, config, datasets, options):
    if options.use_zoom:
        ax.set_xlim(options.zoom_start, options.zoom_end)
        return

    all_x_values = np.concatenate([get_x_values(dataset, config) for dataset in datasets])
    ax.set_xlim(np.min(all_x_values), np.max(all_x_values))


def format_x_axis_without_scientific_notation(ax):
    formatter = ScalarFormatter(useOffset=False, useMathText=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)


def build_output_stem(config, datasets, options):
    if is_single_dataset(datasets):
        output_stem = config.single_output_stem
    else:
        output_stem = config.output_stem

    if options.normalize:
        output_stem += "_normalized"
    if options.stagger:
        output_stem += "_staggered"
    if options.use_zoom:
        output_stem += "_zoomed"
    return output_stem


def save_plot(config, datasets, options):
    output_stem = build_output_stem(config, datasets, options)
    output_dir = datasets[0].path.parent
    output_path_png = output_dir / f"{output_stem}.png"
    output_path_pdf = output_dir / f"{output_stem}.pdf"

    plt.savefig(output_path_png, dpi=300, transparent=True)
    plt.savefig(output_path_pdf, dpi=300)

    return output_path_png, output_path_pdf


def plot_comparison(config, datasets, options):
    configure_matplotlib()
    plt.figure(figsize=FIGURE_SIZE)

    all_raw_intensities, plotted_intensities = plot_datasets(config, datasets, options)
    ax = plt.gca()
    format_axes(ax, config, datasets, options, all_raw_intensities, plotted_intensities)

    if not is_single_dataset(datasets):
        plt.legend()
    plt.tight_layout()

    output_path_png, output_path_pdf = save_plot(config, datasets, options)
    print(f"Plot saved as PNG at: {output_path_png}")
    print(f"Plot saved as PDF at: {output_path_pdf}")

    plt.show()


def main():
    config = prompt_plot_config()
    datasets = prompt_datasets()
    options = prompt_plot_options(config, is_single_dataset(datasets))
    plot_comparison(config, datasets, options)


if __name__ == "__main__":
    main()
