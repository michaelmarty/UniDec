import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import numpy as np

# -----------------------
# Plot settings
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
# =====================================================
files = [
    r"G:\EC\Mp data\MSP1D1.csv",

    # Add additional files here:
    # r"D:\UT Summer 2026\mass photometry\sample_2.csv",
    # r"D:\UT Summer 2026\mass photometry\sample_3.csv",
]


# =====================================================
# LEGEND LABELS
# One label is required for each file
# =====================================================
sample_labels = [
    "MSP1D1",

    # "Sample 2",
    # "Sample 3",
]


# =====================================================
# CUSTOM HEX COLORS
# One color may be supplied for each file.
# Colors are reused if there are more files than colors.
# =====================================================
colors = [
    "#333333",

    # "#00A9A5",
    # "#C76BA0",
    # "#2C7FB8",
]


# =====================================================
# MASS PHOTOMETRY SETTINGS
# =====================================================

# Column containing the calibrated mass values
mass_column = "calibrated_values"

# Unit currently used in the calibrated_values column.
# Use "kDa" if the values are already in kilodaltons.
# Use "Da" if the values need to be divided by 1000.
input_mass_unit = "kDa"

# Histogram bin width in kDa
bin_width_kda = 5.0

# Mass range to plot in kDa.
# Use None to determine a limit automatically from the data.
mass_min_kda = 0.0
mass_max_kda = None

# Set to True to plot only selected events.
# This uses the "selections" column.
selected_only = False

# Annotate the histogram bin containing the most events
annotate_tallest_bin = True

# Number of initial rows to skip before reading the header.
# Normally this should be 0 for a mass-photometer CSV.
skip_rows = 0


# =====================================================
# OUTPUT SETTINGS
# =====================================================
outdir = r"G:\EC\Mp data"
output_filename = "mass_photometry_histograms.pdf"


# -----------------------
# Detect likely delimiter
# -----------------------
def detect_delimiter(
    file_path,
    encodings=("utf-8", "utf-8-sig", "latin1", "utf-16")
):
    """
    Detect whether a file is comma-, tab-, or semicolon-separated.
    """

    for encoding in encodings:
        try:
            with open(
                file_path,
                "r",
                encoding=encoding,
                errors="strict"
            ) as file_handle:
                sample = file_handle.read(4096)

            if "\t" in sample:
                return "\t"

            if "," in sample:
                return ","

            if ";" in sample:
                return ";"

            return None

        except Exception:
            continue

    return None


# -----------------------
# Read mass-photometry CSV
# -----------------------
def read_mass_photometry_file(
    file_path,
    column_name="calibrated_values",
    unit="kDa",
    selected_events_only=False,
    rows_to_skip=0
):
    """
    Read a mass-photometry CSV or Excel file and return mass values in kDa.
    """

    _, extension = os.path.splitext(file_path)
    extension = extension.lower()

    if extension in (".xls", ".xlsx"):
        dataframe = pd.read_excel(
            file_path,
            skiprows=rows_to_skip
        )

    elif extension in (".csv", ".txt", ".tsv"):
        encodings_to_try = [
            "utf-8",
            "utf-8-sig",
            "latin1",
            "utf-16"
        ]

        separator = detect_delimiter(
            file_path,
            encodings=encodings_to_try
        )

        dataframe = None
        last_exception = None

        for encoding in encodings_to_try:
            try:
                if separator is None:
                    dataframe = pd.read_csv(
                        file_path,
                        sep=None,
                        engine="python",
                        encoding=encoding,
                        skiprows=rows_to_skip
                    )
                else:
                    dataframe = pd.read_csv(
                        file_path,
                        sep=separator,
                        encoding=encoding,
                        skiprows=rows_to_skip
                    )

                break

            except Exception as exception:
                last_exception = exception

        if dataframe is None:
            raise RuntimeError(
                f"Could not read file:\n{file_path}\n\n"
                f"Last error: {last_exception}"
            )

    else:
        raise ValueError(
            f"Unsupported file extension: {extension}\n"
            "Supported formats are .csv, .txt, .tsv, .xls, and .xlsx."
        )

    # Clean spaces and byte-order marks from column names
    dataframe.columns = (
        dataframe.columns
        .astype(str)
        .str.strip()
        .str.replace("\ufeff", "", regex=False)
    )

    # Find the mass column without requiring matching capitalization
    column_lookup = {
        column.lower(): column
        for column in dataframe.columns
    }

    actual_mass_column = column_lookup.get(
        column_name.lower()
    )

    if actual_mass_column is None:
        available_columns = "\n".join(
            f"  - {column}"
            for column in dataframe.columns
        )

        raise KeyError(
            f"The mass column '{column_name}' was not found in:\n"
            f"{file_path}\n\n"
            f"Available columns are:\n{available_columns}"
        )

    # Optionally keep only rows marked as selected
    if selected_events_only:
        selection_column = column_lookup.get("selections")

        if selection_column is None:
            raise KeyError(
                "selected_only is True, but the file does not contain "
                "a 'selections' column."
            )

        selection_text = (
            dataframe[selection_column]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
        )

        accepted_selection_values = {
            "1",
            "1.0",
            "true",
            "yes",
            "y",
            "selected",
            "keep"
        }

        selection_mask = selection_text.isin(
            accepted_selection_values
        )

        dataframe = dataframe.loc[selection_mask].copy()

    # Convert the mass column to numeric values
    masses = pd.to_numeric(
        dataframe[actual_mass_column],
        errors="coerce"
    )

    # Replace infinite values and remove missing values
    masses = masses.replace(
        [np.inf, -np.inf],
        np.nan
    )

    number_invalid = int(masses.isna().sum())

    masses = masses.dropna().to_numpy(dtype=float)

    if number_invalid > 0:
        print(
            f"Warning: ignored {number_invalid} missing or non-numeric "
            f"mass value(s) in:\n{file_path}"
        )

    if len(masses) == 0:
        raise ValueError(
            f"No numeric mass values were found in:\n{file_path}"
        )

    # Convert Da to kDa when needed
    if unit.lower() == "da":
        masses = masses / 1000.0

    elif unit.lower() != "kda":
        raise ValueError(
            "input_mass_unit must be either 'kDa' or 'Da'."
        )

    return masses


# -----------------------
# Check settings
# -----------------------
if len(files) == 0:
    raise ValueError("No files were entered in the files list.")

if len(files) != len(sample_labels):
    raise ValueError(
        "The number of files must equal the number of sample labels.\n"
        f"Number of files: {len(files)}\n"
        f"Number of labels: {len(sample_labels)}"
    )

if len(colors) == 0:
    raise ValueError("At least one color must be entered.")

if bin_width_kda <= 0:
    raise ValueError("bin_width_kda must be greater than zero.")


# -----------------------
# Read and filter all files
# -----------------------
datasets = []

for file_path, sample_label in zip(files, sample_labels):

    masses_kda = read_mass_photometry_file(
        file_path=file_path,
        column_name=mass_column,
        unit=input_mass_unit,
        selected_events_only=selected_only,
        rows_to_skip=skip_rows
    )

    number_before_filtering = len(masses_kda)

    # Apply the lower mass limit
    if mass_min_kda is not None:
        masses_kda = masses_kda[
            masses_kda >= mass_min_kda
        ]

    # Apply the upper mass limit
    if mass_max_kda is not None:
        masses_kda = masses_kda[
            masses_kda <= mass_max_kda
        ]

    if len(masses_kda) == 0:
        raise ValueError(
            f"No events remain inside the selected mass range for:\n"
            f"{file_path}"
        )

    number_removed = number_before_filtering - len(masses_kda)

    print("----------------------------------------")
    print(f"Sample: {sample_label}")
    print(f"File: {file_path}")
    print(f"Events read: {number_before_filtering:,}")
    print(f"Events inside plot range: {len(masses_kda):,}")
    print(f"Events excluded by mass limits: {number_removed:,}")
    print(
        f"Observed mass range: "
        f"{np.min(masses_kda):.2f}–{np.max(masses_kda):.2f} kDa"
    )

    datasets.append(
        {
            "file": file_path,
            "label": sample_label,
            "masses": masses_kda
        }
    )


# -----------------------
# Calculate common histogram bins
# -----------------------
all_masses = np.concatenate(
    [dataset["masses"] for dataset in datasets]
)

if mass_min_kda is None:
    plot_minimum = (
        np.floor(np.min(all_masses) / bin_width_kda)
        * bin_width_kda
    )
else:
    plot_minimum = mass_min_kda

if mass_max_kda is None:
    plot_maximum = (
        np.ceil(np.max(all_masses) / bin_width_kda)
        * bin_width_kda
    )
else:
    plot_maximum = mass_max_kda

# Make sure the range contains at least one bin
if plot_maximum <= plot_minimum:
    plot_maximum = plot_minimum + bin_width_kda

number_of_bins = int(
    np.ceil(
        (plot_maximum - plot_minimum)
        / bin_width_kda
    )
)

bin_edges = (
    plot_minimum
    + np.arange(number_of_bins + 1) * bin_width_kda
)


# -----------------------
# Plot setup
# -----------------------
figure_height = max(
    4.0,
    3.3 * len(datasets)
)

fig, axes = plt.subplots(
    len(datasets),
    1,
    figsize=(8, figure_height),
    sharex=True,
    dpi=300
)

if len(datasets) == 1:
    axes = [axes]


# -----------------------
# Plot each sample
# -----------------------
for index, dataset in enumerate(datasets):

    axis = axes[index]
    masses_kda = dataset["masses"]
    sample_label = dataset["label"]
    sample_color = colors[index % len(colors)]

    # Plot mass histogram
    counts, edges, patches = axis.hist(
        masses_kda,
        bins=bin_edges,
        color=sample_color,
        edgecolor="black",
        linewidth=0.4,
        alpha=0.9,
        label=sample_label
    )

    maximum_count = np.max(counts)

    # Annotate the histogram bin with the most events
    if annotate_tallest_bin and maximum_count > 0:

        tallest_bin_index = int(np.argmax(counts))

        tallest_bin_center = (
            edges[tallest_bin_index]
            + edges[tallest_bin_index + 1]
        ) / 2.0

        tallest_bin_count = int(
            counts[tallest_bin_index]
        )

        # Marker at top of tallest bin
        axis.plot(
            tallest_bin_center,
            tallest_bin_count,
            marker="o",
            color="black",
            markersize=4
        )

        # Vertical dashed line
        axis.axvline(
            tallest_bin_center,
            color=sample_color,
            linestyle="--",
            linewidth=0.8,
            alpha=0.7
        )

        # Text annotation
        axis.annotate(
            (
                f"{tallest_bin_center:.1f} kDa\n"
                f"{tallest_bin_count:,} events"
            ),
            xy=(
                tallest_bin_center,
                tallest_bin_count
            ),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            bbox=dict(
                facecolor="white",
                edgecolor="none",
                alpha=0.75,
                pad=1
            )
        )

    # Event-count annotation
    axis.text(
        0.98,
        0.95,
        f"N = {len(masses_kda):,} events",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=10
    )

    # Axis formatting
    axis.set_xlim(
        plot_minimum,
        plot_maximum
    )

    if maximum_count > 0:
        axis.set_ylim(
            0,
            maximum_count * 1.20
        )

    axis.margins(y=0)

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    axis.tick_params(
        axis="y",
        length=3,
        width=0.8
    )

    axis.tick_params(
        axis="x",
        length=3,
        width=0.8
    )

    axis.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        handlelength=1.5
    )


# -----------------------
# Shared axis labels
# -----------------------
fig.text(
    0.04,
    0.5,
    "Events",
    va="center",
    rotation="vertical"
)

axes[-1].set_xlabel("Mass (kDa)")


# -----------------------
# Layout and save
# -----------------------
plt.tight_layout(
    rect=[0.08, 0.05, 0.84, 0.97],
    h_pad=0.8
)

plt.subplots_adjust(
    hspace=0.55
)

os.makedirs(
    outdir,
    exist_ok=True
)

output_path = os.path.join(
    outdir,
    output_filename
)

plt.savefig(
    output_path,
    bbox_inches="tight"
)

print("----------------------------------------")
print(f"Plot saved to:\n{output_path}")

plt.show()