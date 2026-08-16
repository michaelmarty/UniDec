"""Create the replicate TIC and ion-accumulation plots for Figure S1."""

from pathlib import Path
from string import ascii_uppercase
from xml.etree import ElementTree
from zipfile import ZipFile

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np


CHROMATOGRAMS = [
    (
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 1\merged"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_"
            r"set1_merged_21_unidecfiles\TIC_decon.txt"
        ),
        "Rep 1",
        "#0072B2",
    ),
    (
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 2\merged"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_"
            r"set2_merged_21_unidecfiles\TIC_decon.txt"
        ),
        "Rep 2",
        "#D55E00",
    ),
    (
        Path(
            r"Z:\Group Share\BHT\Temp\Paper\DMT"
            r"\Protein Mix 21 Injections Set 3\merged"
            r"\20260327_BHT_agilentc8column_protienmix_CDMS_0o1mgmLs_"
            r"set3_merged_21_unidecfiles\TIC_decon.txt"
        ),
        "Rep 3",
        "#009E73",
    ),
]

SPREADSHEET = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S1"
    r"\Protein Mix Total Ion Scatter Plot.xlsx"
)

TIME_LIMITS = (0.0, 5.5)
FIGURE_SIZE_INCHES = (6.5, 3.25)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S1"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S1.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S1.pdf"

DATA_COLOR = "#1F77B4"
FIT_COLOR = "#1F77B4"
ERROR_COLOR = "#000000"
XML_NAMESPACE = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
}


def load_chromatogram(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load and validate a two-column time/ion-count chromatogram."""
    if not path.is_file():
        raise FileNotFoundError(f"Chromatogram not found: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in: {path}")

    time = data[:, 0]
    ion_count = data[:, 1]
    valid = (
        np.isfinite(time)
        & np.isfinite(ion_count)
        & (time >= TIME_LIMITS[0])
        & (time <= TIME_LIMITS[1])
    )
    if not np.any(valid):
        raise ValueError(f"No finite data in the requested time range: {path}")

    # Remove only physically meaningless negative numerical noise.
    return time[valid], np.clip(ion_count[valid], 0.0, None)


def column_index(cell_reference: str) -> int:
    """Convert the alphabetic portion of an Excel cell reference to an index."""
    letters = "".join(character for character in cell_reference if character.isalpha())
    index = 0
    for character in letters.upper():
        if character not in ascii_uppercase:
            raise ValueError(f"Invalid Excel cell reference: {cell_reference}")
        index = index * 26 + ascii_uppercase.index(character) + 1
    return index - 1


def load_spreadsheet(path: Path) -> tuple[list[str], np.ndarray]:
    """Read the first six columns of the first XLSX sheet without openpyxl.

    Formula cells in the supplied workbook contain cached values, so the XLSX
    Open XML can be read directly with Python's standard library.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Spreadsheet not found: {path}")

    with ZipFile(path) as workbook:
        names = set(workbook.namelist())
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in names:
            shared_root = ElementTree.fromstring(
                workbook.read("xl/sharedStrings.xml")
            )
            for item in shared_root.findall("main:si", XML_NAMESPACE):
                text_parts = item.findall(".//main:t", XML_NAMESPACE)
                shared_strings.append("".join(part.text or "" for part in text_parts))

        worksheet_name = "xl/worksheets/sheet1.xml"
        if worksheet_name not in names:
            raise ValueError(f"First worksheet is missing from: {path}")
        worksheet_root = ElementTree.fromstring(workbook.read(worksheet_name))

        rows: list[list[str | float | None]] = []
        for row_element in worksheet_root.findall(
            ".//main:sheetData/main:row", XML_NAMESPACE
        ):
            row: list[str | float | None] = [None] * 6
            for cell in row_element.findall("main:c", XML_NAMESPACE):
                index = column_index(cell.get("r", ""))
                if index < 0 or index >= len(row):
                    continue

                value_element = cell.find("main:v", XML_NAMESPACE)
                if value_element is None or value_element.text is None:
                    continue
                raw_value = value_element.text
                if cell.get("t") == "s":
                    row[index] = shared_strings[int(raw_value)]
                else:
                    row[index] = float(raw_value)
            rows.append(row)

    if len(rows) < 2:
        raise ValueError(f"Spreadsheet contains no data rows: {path}")

    headers = [str(value) for value in rows[0]]
    try:
        data = np.asarray(rows[1:], dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Expected six numeric data columns in: {path}") from error

    if data.shape != (21, 6) or not np.all(np.isfinite(data)):
        raise ValueError(
            f"Expected 21 rows and 6 finite numeric columns; found {data.shape} in: {path}"
        )
    return headers, data


def style_axis(axis: plt.Axes, panel_letter: str) -> None:
    """Apply shared manuscript styling to one panel."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(1.5)
    axis.spines["bottom"].set_linewidth(1.5)
    axis.tick_params(direction="out", width=1.5, length=6)
    axis.text(
        0.015,
        0.97,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def make_figure() -> plt.Figure:
    """Build and return the complete two-panel Figure S1."""
    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )

    # A: raw chromatogram overlays for the three independent injection sets.
    maximum_ion_count = 0.0
    for path, label, color in CHROMATOGRAMS:
        time, ion_count = load_chromatogram(path)
        maximum_ion_count = max(maximum_ion_count, float(np.max(ion_count)))
        axes[0].plot(
            time,
            ion_count,
            color=color,
            linewidth=1.5,
            label=label,
        )

    style_axis(axes[0], "A")
    axes[0].set_xlim(*TIME_LIMITS)
    axes[0].set_ylim(0.0, maximum_ion_count * 1.18)
    axes[0].set_xticks(np.arange(0.0, 5.1, 1.0))
    axes[0].set_xlabel("Time (min)")
    axes[0].set_ylabel("Ion Count")
    axes[0].legend(
        loc="upper left",
        bbox_to_anchor=(0.12, 1.0),
        frameon=False,
        fontsize=10,
        handlelength=2.0,
        borderaxespad=0.0,
    )

    # B: averages and sample standard deviations from the three replicates.
    _, table = load_spreadsheet(SPREADSHEET)
    injection_number = table[:, 0]
    replicate_values = table[:, 1:4]
    supplied_average = table[:, 4]
    supplied_standard_deviation = table[:, 5]

    calculated_average = np.mean(replicate_values, axis=1)
    calculated_standard_deviation = np.std(replicate_values, axis=1, ddof=1)
    if not np.allclose(supplied_average, calculated_average):
        raise ValueError("Spreadsheet averages do not match replicate columns 2--4")
    if not np.allclose(supplied_standard_deviation, calculated_standard_deviation):
        raise ValueError("Spreadsheet standard deviations do not match columns 2--4")

    axes[1].errorbar(
        injection_number,
        supplied_average,
        yerr=supplied_standard_deviation,
        fmt="o",
        color=DATA_COLOR,
        markerfacecolor=DATA_COLOR,
        markeredgecolor=DATA_COLOR,
        markersize=4.5,
        ecolor=ERROR_COLOR,
        elinewidth=1.0,
        capsize=2.5,
        capthick=1.0,
        linestyle="none",
        label="Mean ± SD",
        zorder=3,
    )

    slope, intercept = np.polyfit(injection_number, supplied_average, 1)
    fit_x = np.array([injection_number.min(), injection_number.max()])
    fit_y = slope * fit_x + intercept
    axes[1].plot(
        fit_x,
        fit_y,
        color=FIT_COLOR,
        linewidth=1.5,
        label="Linear fit",
        zorder=2,
    )

    fitted_values = slope * injection_number + intercept
    residual_sum_squares = np.sum((supplied_average - fitted_values) ** 2)
    total_sum_squares = np.sum((supplied_average - supplied_average.mean()) ** 2)
    r_squared = 1.0 - residual_sum_squares / total_sum_squares
    axes[1].text(
        0.05,
        0.76,
        f"y = {slope:,.1f}x + {intercept:,.1f}\nR² = {r_squared:.5f}",
        transform=axes[1].transAxes,
        fontsize=10,
        ha="left",
        va="top",
    )

    style_axis(axes[1], "B")
    axes[1].set_xlim(0.5, 21.5)
    axes[1].set_ylim(bottom=0.0)
    axes[1].set_xticks([1, 5, 10, 15, 20])
    axes[1].set_xlabel("Number of injections")
    axes[1].set_ylabel("Total Number of Ions")
    axes[1].yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))

    return figure


def main() -> None:
    """Render Figure S1 as a transparent PNG and PDF."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.labelsize": 13,
            "legend.fontsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    figure = make_figure()
    # Avoid bbox_inches="tight" so the saved files remain exactly 7 inches wide.
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
