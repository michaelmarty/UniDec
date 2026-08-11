"""Create the high-throughput antibody ion scatter plots for Figure S11."""

from dataclasses import dataclass
from pathlib import Path
from string import ascii_uppercase
from xml.etree import ElementTree
from zipfile import ZipFile

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter
import numpy as np


@dataclass(frozen=True)
class SampleSource:
    """Spreadsheet and display title for one antibody sample."""

    title: str
    path: Path


SAMPLE_SOURCES = (
    SampleSource(
        title="Waters Intact mAb\nMass Check Standard",
        path=Path(
            r"C:\Users\bht442\The University of Texas at Austin"
            r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
            r"\Figures for Manuscript\Supporting Information Figures\Figure S11"
            r"\Waters mAb HighTotal Ion Scatter Plot.xlsx"
        ),
    ),
    SampleSource(
        title="Human IgG",
        path=Path(
            r"C:\Users\bht442\The University of Texas at Austin"
            r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
            r"\Figures for Manuscript\Supporting Information Figures\Figure S11"
            r"\Human IgG HighTotal Ion Scatter Plot.xlsx"
        ),
    ),
    SampleSource(
        title="Rabbit IgG",
        path=Path(
            r"C:\Users\bht442\The University of Texas at Austin"
            r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
            r"\Figures for Manuscript\Supporting Information Figures\Figure S11"
            r"\Rabbit IgG HighTotal Ion Scatter Plot.xlsx"
        ),
    ),
)

FIGURE_SIZE_INCHES = (7.0, 8.5)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_PNG = OUTPUT_DIRECTORY / "Figure S11.png"
OUTPUT_PDF = OUTPUT_DIRECTORY / "Figure S11.pdf"

DATA_COLOR = "#1F77B4"
FIT_COLOR = "#1F77B4"
ERROR_COLOR = "#000000"
AXIS_STROKE_WIDTH = 1.5
XML_NAMESPACE = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
}


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
    """Read the first six columns of the first XLSX sheet without openpyxl."""
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


def validate_statistics(table: np.ndarray, path: Path) -> None:
    """Verify the supplied means and sample standard deviations."""
    replicate_values = table[:, 1:4]
    supplied_average = table[:, 4]
    supplied_standard_deviation = table[:, 5]
    calculated_average = np.mean(replicate_values, axis=1)
    calculated_standard_deviation = np.std(replicate_values, axis=1, ddof=1)

    if not np.allclose(supplied_average, calculated_average):
        raise ValueError(f"Spreadsheet averages do not match replicates in: {path}")
    if not np.allclose(supplied_standard_deviation, calculated_standard_deviation):
        raise ValueError(
            f"Spreadsheet standard deviations do not match replicates in: {path}"
        )


def style_axis(axis: plt.Axes, panel_letter: str) -> None:
    """Apply the Figure S1 manuscript styling to one panel."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(AXIS_STROKE_WIDTH)
    axis.spines["bottom"].set_linewidth(AXIS_STROKE_WIDTH)
    axis.tick_params(
        direction="out",
        width=AXIS_STROKE_WIDTH,
        length=6,
    )
    axis.text(
        0.025,
        0.96,
        panel_letter,
        transform=axis.transAxes,
        fontsize=16,
        fontweight="bold",
        ha="left",
        va="top",
        zorder=10,
    )


def plot_scatter_panel(
    axis: plt.Axes,
    table: np.ndarray,
    source: SampleSource,
    panel_letter: str,
) -> None:
    """Plot replicate mean, sample SD, linear fit, and fit statistics."""
    injection_number = table[:, 0]
    average = table[:, 4]
    standard_deviation = table[:, 5]

    axis.errorbar(
        injection_number,
        average,
        yerr=standard_deviation,
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
        zorder=3,
    )

    slope, intercept = np.polyfit(injection_number, average, 1)
    fit_x = np.array([injection_number.min(), injection_number.max()])
    fit_y = slope * fit_x + intercept
    axis.plot(
        fit_x,
        fit_y,
        color=FIT_COLOR,
        linewidth=1.5,
        zorder=2,
    )

    fitted_values = slope * injection_number + intercept
    residual_sum_squares = np.sum((average - fitted_values) ** 2)
    total_sum_squares = np.sum((average - average.mean()) ** 2)
    r_squared = 1.0 - residual_sum_squares / total_sum_squares
    intercept_operator = "+" if intercept >= 0 else "-"
    axis.text(
        0.05,
        0.76,
        f"y = {slope:,.1f}x {intercept_operator} {abs(intercept):,.1f}\n"
        f"$R^2$ = {r_squared:.5f}",
        transform=axis.transAxes,
        fontsize=8.5,
        ha="left",
        va="top",
    )

    style_axis(axis, panel_letter)
    axis.set_xlim(0.5, 21.5)
    axis.set_ylim(bottom=0.0)
    axis.set_xticks([1, 5, 10, 15, 20])
    axis.set_xlabel("Number of injections")
    axis.set_ylabel("Total Number of Ions")
    axis.set_title(source.title, fontsize=12, fontweight="bold", pad=7)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    axis.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))


def make_figure(tables: tuple[np.ndarray, ...]) -> plt.Figure:
    """Build and return the complete three-panel Figure S11."""
    figure, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )
    for index, (axis, table, source) in enumerate(
        zip(axes, tables, SAMPLE_SOURCES)
    ):
        plot_scatter_panel(axis, table, source, chr(ord("A") + index))
    return figure


def main() -> None:
    """Load the spreadsheets and render Figure S11 as PNG and PDF."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 12,
            "axes.labelsize": 13,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.transparent": True,
        }
    )

    tables = []
    for source in SAMPLE_SOURCES:
        _, table = load_spreadsheet(source.path)
        validate_statistics(table, source.path)
        tables.append(table)

    figure = make_figure(tuple(tables))
    figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    figure.savefig(OUTPUT_PDF, dpi=OUTPUT_DPI, transparent=True)
    plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
