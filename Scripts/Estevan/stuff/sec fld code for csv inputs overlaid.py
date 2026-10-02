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
rcParams["font.size"] = 12
rcParams["font.weight"] = 200
rcParams["font.sans-serif"] = "Arial"
rcParams["axes.titleweight"] = 700
rcParams["axes.titlesize"] = 16
rcParams["axes.linewidth"] = 1
rcParams["xtick.major.width"] = 1
rcParams["ytick.major.width"] = 1


# -----------------------
# Files and sample labels
# -----------------------
files = [
    r"D:\UT Fall 2026\B2 akta files\MSPtit\09232026_B2_msp_tit_10_elut 001.csv",
    #r"Z:\Group Share\Group\Papers\FluoroND\PaperData\SEC\MSP1D1sfGFP DMPC Nanodiscs\MSP1D1 sfGFP DMPC rep1.csv",
    #r"Z:\Group Share\Group\Papers\FluoroND\PaperData\SEC\MSP1D1mCherry DMPC Nanodiscs\MSP1D1 mCherry DMPC rep3.csv",
]

sample_labels = [
    ".9 mmol MSP",
    #"MSP1D1-sfGFP",
    #"MSP1D1-mCherry",
]


# Each signal uses the same color in every sample panel.
signal_colors = {
    "uv": "#333333",
    "analog_1": "#00A9A5",
    "analog_2": "#C76BA0",
}

standard_volumes = [13.81, 16.82, 18.09, 19.05, 20.79]
standard_color = "#8C8C8C"

x_min = 0.0
x_max = 24.0
# Keep the UV trace in its original mAU units.
normalize_uv = False

outdir = r"D:\UT Fall 2026\B2 akta files\MSPtit"
output_filename = "sec_uv_and_analog_signals_overlaid.pdf"


def detect_delimiter(file_path, encodings=("utf-8", "utf-8-sig", "latin1", "utf-16")):
    """Detect whether a text export is tab- or comma-delimited."""
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
    """Read A:B, K:L, and M:N from a UNICORN export."""
    _, extension = os.path.splitext(file_path)
    extension = extension.lower()
    selected_columns = [0, 1, 14, 15, 12, 13]

    if extension in (".xls", ".xlsx"):
        return pd.read_excel(
            file_path,
            skiprows=3,
            header=None,
            usecols=selected_columns,
        )

    if extension != ".csv":
        raise ValueError(
            f"Unsupported extension: {extension}. Use .xls, .xlsx, or .csv."
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


def _clean_signal(volume, signal, max_volume):
    """Clean one x/y pair without dropping valid rows from other signals."""
    clean = pd.DataFrame(
        {
            "ml": pd.to_numeric(volume, errors="coerce"),
            "signal": pd.to_numeric(signal, errors="coerce"),
        }
    )
    clean = clean.dropna(subset=["ml", "signal"])
    clean = clean[clean["ml"] <= max_volume]
    return clean.sort_values("ml").reset_index(drop=True)


def read_sec_file(file_path, csv_sep=None, max_volume=20.0):
    """Return UV, Analog In 1, and Analog In 2 traces."""
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


def _style_y_axis(axis, color, label, side="left"):
    axis.set_ylabel(label, color=color)
    axis.tick_params(axis="y", colors=color, length=3, width=0.8)
    axis.spines[side].set_color(color)
    axis.spines["top"].set_visible(False)


def plot_overlaid_chromatograms(file_paths, labels):
    """Overlay all three signals in one panel for each sample."""
    if len(file_paths) != len(labels):
        raise ValueError("files and sample_labels must contain the same number of items.")
    if not file_paths:
        raise ValueError("The files list is empty.")

    figure, base_axes = plt.subplots(
        len(file_paths),
        1,
        figsize=(10, 3.2 * len(file_paths)),
        sharex=True,
        dpi=300,
        squeeze=False,
    )
    base_axes = base_axes[:, 0]

    for sample_index, (file_path, label) in enumerate(zip(file_paths, labels)):
        traces = read_sec_file(file_path, max_volume=x_max)
        uv_axis = base_axes[sample_index]
        analog_1_axis = uv_axis.twinx()
        analog_2_axis = uv_axis.twinx()

        # Move the third y-axis outward so both raw analog scales are readable.
        analog_2_axis.spines["right"].set_position(("outward", 70))
        analog_2_axis.spines["right"].set_visible(True)

        axis_specs = [
            (
                uv_axis,
                traces["uv"],
                signal_colors["uv"],
                "UV",
                "UV absorbance (mAU)",
            ),
            (
                analog_1_axis,
                traces["analog_1"],
                signal_colors["analog_1"],
                "Analog In 1",
                "Analog In 1 (mV)",
            ),
            (
                analog_2_axis,
                traces["analog_2"],
                signal_colors["analog_2"],
                "Analog In 2",
                "Analog In 2 (mV)",
            ),
        ]

        plotted_lines = []
        for axis, trace, color, legend_label, y_label in axis_specs:
            if trace.empty:
                continue
            line, = axis.plot(
                trace["ml"].to_numpy(),
                trace["signal"].to_numpy(),
                color=color,
                linewidth=1.8,
                label=legend_label,
            )
            plotted_lines.append(line)
            side = "left" if axis is uv_axis else "right"
            _style_y_axis(axis, color, y_label, side=side)
            axis.margins(y=0.08)

        uv_axis.set_xlim(x_min, x_max)
        uv_axis.tick_params(axis="x", length=3, width=0.8)
        uv_axis.spines["top"].set_visible(False)
        uv_axis.set_title(label, loc="left", pad=7)

        if normalize_uv and not traces["uv"].empty:
            uv_axis.set_ylim(0, 1.08)
            uv_axis.set_yticks([0, 1])

        if plotted_lines:
            uv_axis.legend(
                handles=plotted_lines,
                loc="upper right",
                frameon=False,
                handlelength=1.8,
            )

        if sample_index == 0:
            top_axis = uv_axis.secondary_xaxis("top")
            top_axis.set_xticks(standard_volumes)
            top_axis.set_xticklabels([])
            top_axis.tick_params(
                axis="x", length=5, width=1.0, colors=standard_color
            )
            top_axis.spines["top"].set_visible(False)
            top_axis.set_xlabel("SEC standards", color=standard_color, labelpad=6)

    base_axes[-1].set_xlabel("Retention Volume (mL)")

    # Leave extra room on the right for the two separate analog y-axes.
    figure.tight_layout(rect=[0.03, 0.03, 0.80, 0.98], h_pad=1.0)
    figure.subplots_adjust(hspace=0.50)
    return figure, base_axes


def main():
    figure, _ = plot_overlaid_chromatograms(files, sample_labels)
    os.makedirs(outdir, exist_ok=True)
    output_path = os.path.join(outdir, output_filename)
    figure.savefig(output_path, bbox_inches="tight")
    print(f"Saved: {output_path}")
    plt.show()


if __name__ == "__main__":
    main()
