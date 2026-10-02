"""
Mass-photometry histogram plotter with Excel export.

The same histogram data used to make the figure are exported to an Excel
workbook. This allows the results to be replotted later as:

    x = mass_kDa
    y = events

or:

    x = mass_kDa
    y = intensity

Required packages:
    pip install pandas numpy matplotlib openpyxl
"""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams


# -----------------------
# Plot style
# -----------------------
rcParams["ps.useafm"] = True
rcParams["ps.fonttype"] = 42
rcParams["pdf.fonttype"] = 42
rcParams["lines.linewidth"] = 1
rcParams["errorbar.capsize"] = 3
rcParams["patch.force_edgecolor"] = True
rcParams["lines.markersize"] = 7
rcParams["font.size"] = 12
rcParams["font.weight"] = 200
rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.titleweight"] = 700
rcParams["axes.titlesize"] = 16
rcParams["axes.linewidth"] = 1
rcParams["xtick.major.width"] = 1
rcParams["ytick.major.width"] = 1


# =====================================================
# FILES TO PLOT
# Put one file path on each active line.
# Keep the r before the quotation mark for Windows paths.
# =====================================================
files = [
    r"G:\EC\Mp data\MSP1D1_sfCherryS.csv",

    # Add more files here:
    # r"D:\UT Summer 2026\mass photometry\sample_2.csv",
]


# =====================================================
# SAMPLE LABELS
# Enter exactly one label for each file above, in the same order.
# =====================================================
sample_labels = [
    "Sample 1",

    # "Sample 2",
]


# =====================================================
# CUSTOM HEX COLORS
# Colors are reused if there are more files than colors.
# =====================================================
colors = [
    "#333333",
    "#00A9A5",
    "#C76BA0",
    "#2C7FB8",
]


# =====================================================
# MASS-PHOTOMETRY AND HISTOGRAM SETTINGS
# =====================================================

# Column containing calibrated mass values.
mass_column = "calibrated_values"

# Unit used in calibrated_values:
#   "kDa" = values are already in kDa
#   "Da"  = values will be divided by 1000
input_mass_unit = "kDa"

# Width of each histogram bin in kDa.
bin_width_kda = 5.0

# Adjustable x-axis limits in kDa.
# These limits also define the mass range exported to Excel.
# Use None for either value to determine that limit automatically.
x_axis_min_kda = 0.0
x_axis_max_kda = 1000.0

# Set to True to keep only rows marked as selected in the selections column.
selected_only = False

# Mark and label the histogram bin containing the most events.
annotate_tallest_bin = True

# Number of rows before the column-header row. Usually 0.
skip_rows = 0


# =====================================================
# OUTPUT SETTINGS
# =====================================================
outdir = r"G:\EC\Mp data"

output_figure_filename = "MP dat 2"
output_excel_filename = "mass_photometry_histogram_data.xlsx"

# Set to False when running without a graphical display.
show_plot = True


# -----------------------
# Read a mass-photometry file
# -----------------------
def read_mass_photometry_file(file_path):
    """Read a comma- or tab-separated file and return mass values in kDa."""
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File not found:\n{file_path}")

    dataframe = None
    last_error = None

    # sep=None detects comma-, tab-, and semicolon-delimited files.
    for encoding in ("utf-8", "utf-8-sig", "latin1", "utf-16"):
        try:
            dataframe = pd.read_csv(
                file_path,
                sep=None,
                engine="python",
                encoding=encoding,
                skiprows=skip_rows,
            )
            break
        except Exception as error:
            last_error = error

    if dataframe is None:
        raise RuntimeError(
            f"Could not read:\n{file_path}\n\nLast error: {last_error}"
        ) from last_error

    # Remove spaces and hidden byte-order marks from column names.
    dataframe.columns = (
        dataframe.columns
        .astype(str)
        .str.strip()
        .str.replace("\ufeff", "", regex=False)
    )

    # Match column names without requiring identical capitalization.
    column_lookup = {
        column.lower(): column
        for column in dataframe.columns
    }

    actual_mass_column = column_lookup.get(mass_column.lower())

    if actual_mass_column is None:
        available_columns = "\n".join(
            f"  - {column}"
            for column in dataframe.columns
        )
        raise KeyError(
            f"Column '{mass_column}' was not found in:\n{file_path}\n\n"
            f"Available columns are:\n{available_columns}"
        )

    # Optionally retain only events marked as selected.
    if selected_only:
        selection_column = column_lookup.get("selections")

        if selection_column is None:
            raise KeyError(
                "selected_only is True, but no 'selections' column was found."
            )

        accepted_values = {
            "1",
            "1.0",
            "true",
            "yes",
            "y",
            "selected",
            "keep",
        }

        selection_mask = (
            dataframe[selection_column]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(accepted_values)
        )

        dataframe = dataframe.loc[selection_mask].copy()

    masses = pd.to_numeric(
        dataframe[actual_mass_column],
        errors="coerce",
    )

    masses = (
        masses
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .to_numpy(dtype=float)
    )

    if len(masses) == 0:
        raise ValueError(
            f"No numeric values were found in '{mass_column}' for:\n{file_path}"
        )

    if input_mass_unit.lower() == "da":
        masses = masses / 1000.0
    elif input_mass_unit.lower() != "kda":
        raise ValueError("input_mass_unit must be either 'kDa' or 'Da'.")

    return masses


# -----------------------
# Create valid Excel worksheet names
# -----------------------
def make_sheet_name(label, used_names):
    """Return a valid, unique Excel worksheet name."""
    name = re.sub(r"[\\/*?:\[\]]", "_", str(label)).strip()
    name = name or "Sample"
    name = name[:31]

    candidate = name
    number = 2

    while candidate.lower() in used_names:
        suffix = f"_{number}"
        candidate = name[: 31 - len(suffix)] + suffix
        number += 1

    used_names.add(candidate.lower())
    return candidate


# -----------------------
# Apply simple formatting to Excel sheets
# -----------------------
def format_excel_sheet(worksheet):
    """Format headers, widths, filters, and numeric columns."""
    from openpyxl.styles import Alignment, Font, PatternFill

    header_fill = PatternFill("solid", fgColor="333333")
    header_font = Font(color="FFFFFF", bold=True)

    for cell in worksheet[1]:
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center")

    worksheet.freeze_panes = "A2"
    worksheet.auto_filter.ref = worksheet.dimensions

    for column_cells in worksheet.columns:
        column_letter = column_cells[0].column_letter
        maximum_length = 0

        for cell in column_cells:
            if cell.value is not None:
                maximum_length = max(maximum_length, len(str(cell.value)))

        worksheet.column_dimensions[column_letter].width = min(
            max(maximum_length + 2, 12),
            28,
        )

    header_map = {
        cell.value: cell.column
        for cell in worksheet[1]
    }

    decimal_columns = {
        "mass_bin_left_kDa",
        "mass_bin_right_kDa",
        "mass_kDa",
        "intensity",
        "normalized_intensity",
        "fraction_of_events",
    }

    for header in decimal_columns:
        column_number = header_map.get(header)
        if column_number is not None:
            for row_number in range(2, worksheet.max_row + 1):
                worksheet.cell(
                    row=row_number,
                    column=column_number,
                ).number_format = "0.000"

    events_column = header_map.get("events")
    if events_column is not None:
        for row_number in range(2, worksheet.max_row + 1):
            worksheet.cell(
                row=row_number,
                column=events_column,
            ).number_format = "0"


# -----------------------
# Main analysis
# -----------------------
def main():
    # Validate user settings.
    if len(files) == 0:
        raise ValueError("Add at least one file path to the files list.")

    if len(files) != len(sample_labels):
        raise ValueError(
            "There must be exactly one sample label for each input file."
        )

    if len(colors) == 0:
        raise ValueError("Add at least one color to the colors list.")

    if bin_width_kda <= 0:
        raise ValueError("bin_width_kda must be greater than zero.")

    if (
        x_axis_min_kda is not None
        and x_axis_max_kda is not None
        and x_axis_max_kda <= x_axis_min_kda
    ):
        raise ValueError(
            "x_axis_max_kda must be greater than x_axis_min_kda."
        )

    # Read and filter all samples.
    datasets = []
    raw_event_tables = []

    for file_path, sample_label in zip(files, sample_labels):
        all_masses_kda = read_mass_photometry_file(file_path)
        plotted_masses_kda = all_masses_kda.copy()

        if x_axis_min_kda is not None:
            plotted_masses_kda = plotted_masses_kda[
                plotted_masses_kda >= x_axis_min_kda
            ]

        if x_axis_max_kda is not None:
            plotted_masses_kda = plotted_masses_kda[
                plotted_masses_kda <= x_axis_max_kda
            ]

        if len(plotted_masses_kda) == 0:
            print(
                f"Warning: no events remain within the x-axis limits for "
                f"'{sample_label}'. This sample will be skipped."
            )
            continue

        datasets.append(
            {
                "file": file_path,
                "label": sample_label,
                "masses": plotted_masses_kda,
                "events_read": len(all_masses_kda),
            }
        )

        raw_event_tables.append(
            pd.DataFrame(
                {
                    "sample": sample_label,
                    "event_number": np.arange(
                        1,
                        len(plotted_masses_kda) + 1,
                    ),
                    "mass_kDa": plotted_masses_kda,
                }
            )
        )

    if len(datasets) == 0:
        raise ValueError(
            "No samples contained events within the selected x-axis limits."
        )

    # Build common histogram bins so all samples use the same mass scale.
    all_plotted_masses = np.concatenate(
        [dataset["masses"] for dataset in datasets]
    )

    if x_axis_min_kda is None:
        plot_minimum = (
            np.floor(np.min(all_plotted_masses) / bin_width_kda)
            * bin_width_kda
        )
    else:
        plot_minimum = float(x_axis_min_kda)

    if x_axis_max_kda is None:
        plot_maximum = (
            np.ceil(np.max(all_plotted_masses) / bin_width_kda)
            * bin_width_kda
        )
    else:
        plot_maximum = float(x_axis_max_kda)

    if plot_maximum <= plot_minimum:
        plot_maximum = plot_minimum + bin_width_kda

    bin_edges = np.arange(
        plot_minimum,
        plot_maximum,
        bin_width_kda,
        dtype=float,
    )

    if len(bin_edges) == 0 or not np.isclose(bin_edges[0], plot_minimum):
        bin_edges = np.insert(bin_edges, 0, plot_minimum)

    if not np.isclose(bin_edges[-1], plot_maximum):
        bin_edges = np.append(bin_edges, plot_maximum)

    if len(bin_edges) < 2:
        bin_edges = np.array(
            [plot_minimum, plot_maximum],
            dtype=float,
        )

    # Create the figure.
    figure_height = max(4.0, 3.3 * len(datasets))

    figure, axes = plt.subplots(
        len(datasets),
        1,
        figsize=(8, figure_height),
        sharex=True,
        dpi=300,
    )

    if len(datasets) == 1:
        axes = [axes]

    histogram_tables = []

    # Plot each sample and build its Excel table from the same histogram.
    for index, dataset in enumerate(datasets):
        axis = axes[index]
        sample_label = dataset["label"]
        masses_kda = dataset["masses"]
        sample_color = colors[index % len(colors)]

        counts, edges = np.histogram(
            masses_kda,
            bins=bin_edges,
        )

        bin_left = edges[:-1]
        bin_right = edges[1:]
        bin_center = (bin_left + bin_right) / 2.0
        bin_widths = bin_right - bin_left

        maximum_count = int(np.max(counts)) if len(counts) else 0
        total_events = int(np.sum(counts))

        if maximum_count > 0:
            normalized_intensity = counts.astype(float) / maximum_count
        else:
            normalized_intensity = np.zeros_like(counts, dtype=float)

        if total_events > 0:
            fraction_of_events = counts.astype(float) / total_events
        else:
            fraction_of_events = np.zeros_like(counts, dtype=float)

        # This table is the export-file maker. It uses the exact same bins
        # and counts as the histogram shown in the figure.
        histogram_table = pd.DataFrame(
            {
                "sample": sample_label,
                "mass_bin_left_kDa": bin_left,
                "mass_bin_right_kDa": bin_right,
                "mass_kDa": bin_center,
                "events": counts.astype(int),
                # "intensity" is normalized event count, with the tallest
                # histogram bin equal to 1.0.
                "intensity": normalized_intensity,
                "normalized_intensity": normalized_intensity,
                "fraction_of_events": fraction_of_events,
            }
        )

        histogram_tables.append((sample_label, histogram_table))

        # Plot Mass (kDa) versus Events.
        axis.bar(
            bin_left,
            counts,
            width=bin_widths,
            align="edge",
            color=sample_color,
            edgecolor="black",
            linewidth=0.4,
            alpha=0.9,
            label=sample_label,
        )

        if annotate_tallest_bin and maximum_count > 0:
            peak_index = int(np.argmax(counts))
            peak_mass = bin_center[peak_index]
            peak_count = int(counts[peak_index])

            axis.plot(
                peak_mass,
                peak_count,
                marker="o",
                color="black",
                markersize=4,
            )

            axis.axvline(
                peak_mass,
                color=sample_color,
                linestyle="--",
                linewidth=0.8,
                alpha=0.7,
            )

            axis.annotate(
                f"{peak_mass:.1f} kDa\n{peak_count:,} events",
                xy=(peak_mass, peak_count),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                bbox=dict(
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.75,
                    pad=1,
                ),
            )

        axis.text(
            0.98,
            0.95,
            f"N = {total_events:,} events",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=10,
        )

        axis.set_xlim(plot_minimum, plot_maximum)

        if maximum_count > 0:
            axis.set_ylim(0, maximum_count * 1.25)

        axis.margins(y=0)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(axis="y", length=3, width=0.8)
        axis.tick_params(axis="x", length=3, width=0.8)

        axis.legend(
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            handlelength=1.5,
        )

        print("----------------------------------------")
        print(f"Sample: {sample_label}")
        print(f"Events read: {dataset['events_read']:,}")
        print(f"Events plotted/exported: {total_events:,}")
        print(
            f"Plotted mass range: {plot_minimum:.1f} to "
            f"{plot_maximum:.1f} kDa"
        )

    # Shared labels and layout.
    figure.text(
        0.04,
        0.5,
        "Events",
        va="center",
        rotation="vertical",
    )

    axes[-1].set_xlabel("Mass (kDa)")

    plt.tight_layout(
        rect=[0.08, 0.05, 0.84, 0.97],
        h_pad=0.8,
    )

    plt.subplots_adjust(hspace=0.55)

    # Save the figure and the Excel export.
    os.makedirs(outdir, exist_ok=True)

    figure_output_path = os.path.join(
        outdir,
        output_figure_filename,
    )

    excel_output_path = os.path.join(
        outdir,
        output_excel_filename,
    )

    plt.savefig(
        figure_output_path,
        bbox_inches="tight",
        dpi=300,
    )

    combined_histogram_data = pd.concat(
        [table for _, table in histogram_tables],
        ignore_index=True,
    )

    raw_event_data = pd.concat(
        raw_event_tables,
        ignore_index=True,
    )

    settings_table = pd.DataFrame(
        {
            "setting": [
                "mass_column",
                "input_mass_unit",
                "bin_width_kDa",
                "x_axis_min_kDa",
                "x_axis_max_kDa",
                "selected_only",
            ],
            "value": [
                mass_column,
                input_mass_unit,
                bin_width_kda,
                x_axis_min_kda,
                x_axis_max_kda,
                selected_only,
            ],
        }
    )

    with pd.ExcelWriter(
        excel_output_path,
        engine="openpyxl",
    ) as writer:
        settings_table.to_excel(
            writer,
            sheet_name="settings",
            index=False,
        )

        combined_histogram_data.to_excel(
            writer,
            sheet_name="combined_data",
            index=False,
        )

        raw_event_data.to_excel(
            writer,
            sheet_name="raw_event_masses",
            index=False,
        )

        used_sheet_names = {
            "settings",
            "combined_data",
            "raw_event_masses",
        }

        for sample_label, histogram_table in histogram_tables:
            sheet_name = make_sheet_name(
                sample_label,
                used_sheet_names,
            )

            histogram_table.to_excel(
                writer,
                sheet_name=sheet_name,
                index=False,
            )

        for worksheet in writer.book.worksheets:
            format_excel_sheet(worksheet)

    print("----------------------------------------")
    print(f"Figure saved to:\n{figure_output_path}")
    print(f"Excel data saved to:\n{excel_output_path}")
    print(
        "For another plotting script, use 'mass_kDa' for x and either "
        "'events' or 'intensity' for y."
    )

    if show_plot:
        plt.show()
    else:
        plt.close(figure)


if __name__ == "__main__":
    main()