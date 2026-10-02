import os

import matplotlib.pyplot as plt
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
rcParams["patch.facecolor"] = "b"
rcParams["lines.markersize"] = 7
rcParams["font.size"] = 12
rcParams["font.weight"] = 200
rcParams["font.sans-serif"] = "Arial"
rcParams["axes.titleweight"] = 700
rcParams["axes.titlesize"] = 16
rcParams["axes.linewidth"] = 1
rcParams["xtick.major.width"] = 1
rcParams["ytick.major.width"] = 1


# -----------------------
# Files and labels
# -----------------------
files = [
    r"D:\UT Fall 2026\B2 akta files\MSPtit\09232026_B2_msp_tit_200Elut 001.csv",
    #r"Z:\Group Share\Group\Papers\FluoroND\PaperData\SEC\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep1.csv",
    #r"Z:\Group Share\Group\Papers\FluoroND\PaperData\SEC\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep3.csv",
]

sample_labels = [
    "B2 10X MSP",
    #"MSP1D1-sfGFP",
    #"MSP1D1-mCherry",
]

colors = [
    "#333333",
    "#00A9A5",
    "#C76BA0",
]


# -----------------------
# Plot settings
# -----------------------
standard_volumes = [13.81, 16.82, 18.09, 19.05, 20.79]
standard_color = "#8C8C8C"

x_min = 0.0
x_max = 24.0

# UV keeps the behavior of the original script. The analog channels stay raw.
normalize_uv = True

outdir = r"D:\UT Fall 2026\B2 akta files\MSPtit"
output_filename = "sec_uv_and_analog_signals_test1.pdf"


# Column positions in the UNICORN export (zero-based Python indices).
# A:B = UV, K:L = Analog In 1, M:N = Analog In 2.
SIGNAL_COLUMNS = {
    "uv": (0, 1),
    "analog_1": (10, 11),
    "analog_2": (12, 13),
}


def detect_delimiter(file_path, encodings=("utf-8", "utf-8-sig", "latin1", "utf-16")):
    """Detect a tab- or comma-delimited text export."""
    for encoding in encodings:
        try:
            with open(file_path, "r", encoding=encoding, errors="strict") as handle:
                sample = handle.read(4096)
            if "\t" in sample:
                return "\t"
            if "," in sample:
                return ","
            return None
        except (OSError, UnicodeError):
            continue
    return None


def _read_selected_columns(file_path, csv_sep=None):
    """Read only the six columns needed for the three signals."""
    _, extension = os.path.splitext(file_path)
    extension = extension.lower()
    selected_columns = [0, 1, 10, 11, 12, 13]

    if extension in (".xls", ".xlsx"):
        return pd.read_excel(
            file_path,
            skiprows=3,
            header=None,
            usecols=selected_columns,
        )

    if extension != ".csv":
        raise ValueError(
            f"Unsupported file extension: {extension}. "
            "Supported extensions are .xls, .xlsx, and .csv."
        )

    encodings = ("utf-8", "utf-8-sig", "latin1", "utf-16")
    if csv_sep is None:
        csv_sep = detect_delimiter(file_path, encodings)

    attempts = []
    if csv_sep is not None:
        attempts.extend(
            {"sep": csv_sep, "encoding": encoding} for encoding in encodings
        )
    attempts.extend(
        {"sep": None, "engine": "python", "encoding": encoding}
        for encoding in encodings
    )

    last_error = None
    for options in attempts:
        try:
            return pd.read_csv(
                file_path,
                skiprows=3,
                header=None,
                usecols=selected_columns,
                **options,
            )
        except Exception as error:
            last_error = error

    raise RuntimeError(
        f"Failed to read CSV file '{file_path}'. Last error: {last_error}"
    ) from last_error


def _clean_signal(volume, signal, max_volume=24.0):
    """Convert one x/y pair to numeric values and clean it independently."""
    clean = pd.DataFrame(
        {
            "ml": pd.to_numeric(volume, errors="coerce"),
            "signal": pd.to_numeric(signal, errors="coerce"),
        }
    )
    clean = clean.dropna(subset=["ml", "signal"])
    clean = clean[clean["ml"] <= max_volume]
    return clean.sort_values("ml").reset_index(drop=True)


def read_sec_file(file_path, csv_sep=None, max_volume=24.0):
    """Return independently cleaned UV and analog traces from one export."""
    raw = _read_selected_columns(file_path, csv_sep=csv_sep)
    raw.columns = [
        "uv_ml",
        "uv_mau",
        "analog_1_ml",
        "analog_1_mv",
        "analog_2_ml",
        "analog_2_mv",
    ]

    traces = {
        "uv": _clean_signal(raw["uv_ml"], raw["uv_mau"], max_volume),
        "analog_1": _clean_signal(
            raw["analog_1_ml"], raw["analog_1_mv"], max_volume
        ),
        "analog_2": _clean_signal(
            raw["analog_2_ml"], raw["analog_2_mv"], max_volume
        ),
    }

    if normalize_uv and not traces["uv"].empty:
        max_abs = traces["uv"]["signal"].abs().max()
        if pd.notna(max_abs) and max_abs != 0:
            traces["uv"]["signal"] = traces["uv"]["signal"] / max_abs

    return traces


def _format_axis(axis):
    axis.set_xlim(x_min, x_max)
    axis.margins(y=0.05)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(axis="both", length=3, width=0.8)


def plot_sec_chromatograms(file_paths, labels, sample_colors):
    """Create one UV/Analog 1/Analog 2 three-panel block per sample."""
    if not (len(file_paths) == len(labels)):
        raise ValueError("files and sample_labels must contain the same number of items.")
    if not file_paths:
        raise ValueError("The files list is empty.")

    panels_per_sample = 3
    row_count = len(file_paths) * panels_per_sample
    figure, axes = plt.subplots(
        row_count,
        1,
        figsize=(8, 2.0 * row_count),
        sharex=True,
        dpi=300,
        squeeze=False,
    )
    axes = axes[:, 0]

    for sample_index, (file_path, label) in enumerate(zip(file_paths, labels)):
        traces = read_sec_file(file_path, max_volume=x_max)
        color = sample_colors[sample_index % len(sample_colors)]
        uv_axis = axes[sample_index * panels_per_sample]
        analog_1_axis = axes[sample_index * panels_per_sample + 1]
        analog_2_axis = axes[sample_index * panels_per_sample + 2]

        panel_specs = [
            (uv_axis, traces["uv"], "Normalized UV absorbance", True),
            (analog_1_axis, traces["analog_1"], "Analog In 1 (mV)", False),
            (analog_2_axis, traces["analog_2"], "Analog In 2 (mV)", False),
        ]

        for axis, trace, y_label, fill_trace in panel_specs:
            if trace.empty:
                axis.text(
                    0.5,
                    0.5,
                    "No valid data",
                    transform=axis.transAxes,
                    ha="center",
                    va="center",
                    color="0.45",
                )
            else:
                x = trace["ml"].to_numpy()
                y = trace["signal"].to_numpy()
                axis.plot(x, y, color=color, linewidth=1.8)
                if fill_trace:
                    axis.fill_between(x, y, 0, color=color, alpha=0.9, linewidth=0)

            axis.set_ylabel(y_label)
            _format_axis(axis)

        uv_axis.set_title(label, loc="left", color=color, pad=7)
        if normalize_uv:
            uv_axis.set_ylim(0, 1.08)
            uv_axis.set_yticks([0, 1])

        # Add the SEC standard positions above the first panel only.
        if sample_index == 0:
            top_axis = uv_axis.secondary_xaxis("top")
            top_axis.set_xticks(standard_volumes)
            top_axis.set_xticklabels([])
            top_axis.tick_params(
                axis="x", length=5, width=1.0, colors=standard_color
            )
            top_axis.spines["top"].set_visible(False)
            top_axis.set_xlabel("SEC standards", color=standard_color, labelpad=6)

    axes[-1].set_xlabel("Retention Volume (mL)")
    figure.tight_layout(rect=[0.07, 0.03, 0.98, 0.98], h_pad=0.8)
    figure.subplots_adjust(hspace=0.45)
    return figure, axes


def main():
    figure, _ = plot_sec_chromatograms(files, sample_labels, colors)
    os.makedirs(outdir, exist_ok=True)
    output_path = os.path.join(outdir, output_filename)
    figure.savefig(output_path, bbox_inches="tight")
    print(f"Saved: {output_path}")
    plt.show()


if __name__ == "__main__":
    main()
