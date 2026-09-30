"""Recreate Figure S4's protein mixture bar charts for the presentation."""

from pathlib import Path
from xml.etree import ElementTree
from zipfile import ZipFile

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter
import numpy as np


SPREADSHEET = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Papers\Muliple Run Summing Reversed Phase LC-CD-MS"
    r"\Figures for Manuscript\Supporting Information Figures\Figure S4"
    r"\0.1 mz bin and Frequency Correction Test.xlsx"
)

PROTEINS = [
    ("Cytochrome C", 12360.08, "B"),
    ("Myoglobin", 16951.48, "C"),
    ("Carbonic Anhydrase", 29024.63, "D"),
]

FIGURE_SIZE_INCHES = (6.0, 3.5)
OUTPUT_DPI = 300
OUTPUT_DIRECTORY = Path(
    r"C:\Users\bht442\The University of Texas at Austin"
    r"\MartyLab - General\Presentations\Brian\4th Year Presentation\Figures"
)
OUTPUT_PNG = OUTPUT_DIRECTORY / "Protein Mixture Bar Charts.png"

BAR_COLOR = "#1F77B4"
BAR_EDGE_COLOR = "#000000"
THEORETICAL_COLOR = "#000000"
AXIS_STROKE_WIDTH = 0.8
XML_NAMESPACE = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
}


def read_shared_strings(workbook: ZipFile) -> list[str]:
    """Read plain and rich-text strings from an XLSX archive."""
    if "xl/sharedStrings.xml" not in workbook.namelist():
        return []
    root = ElementTree.fromstring(workbook.read("xl/sharedStrings.xml"))
    strings = []
    for item in root.findall("main:si", XML_NAMESPACE):
        text_parts = item.findall(".//main:t", XML_NAMESPACE)
        strings.append("".join(part.text or "" for part in text_parts))
    return strings


def read_first_sheet_cells(path: Path) -> dict[str, str | float]:
    """Return cached values from the first worksheet without openpyxl."""
    if not path.is_file():
        raise FileNotFoundError(f"Spreadsheet not found: {path}")

    with ZipFile(path) as workbook:
        shared_strings = read_shared_strings(workbook)
        worksheet_name = "xl/worksheets/sheet1.xml"
        if worksheet_name not in workbook.namelist():
            raise ValueError(f"First worksheet is missing from: {path}")
        root = ElementTree.fromstring(workbook.read(worksheet_name))

    cells: dict[str, str | float] = {}
    for cell in root.findall(".//main:sheetData/main:row/main:c", XML_NAMESPACE):
        reference = cell.get("r")
        value_element = cell.find("main:v", XML_NAMESPACE)
        if reference is None or value_element is None or value_element.text is None:
            continue
        raw_value = value_element.text
        if cell.get("t") == "s":
            cells[reference] = shared_strings[int(raw_value)]
        else:
            cells[reference] = float(raw_value)
    return cells


def load_figure_data(path: Path) -> tuple[list[str], np.ndarray]:
    """Load the two header rows and four treatment rows from columns A--D."""
    cells = read_first_sheet_cells(path)
    treatments = [str(cells[f"A{row}"]) for row in range(3, 7)]
    values = np.asarray(
        [
            [float(cells[f"{column}{row}"]) for column in ("B", "C", "D")]
            for row in range(3, 7)
        ],
        dtype=float,
    )

    if values.shape != (4, 3) or not np.all(np.isfinite(values)):
        raise ValueError(
            f"Expected four finite treatment rows in columns B--D; found {values.shape}"
        )
    for protein, theoretical_mass, column in PROTEINS:
        header = str(cells.get(f"{column}2", ""))
        if protein not in header or f"{theoretical_mass:.2f}" not in header:
            raise ValueError(
                f"Unexpected column {column} header for {protein}: {header!r}"
            )
    return treatments, values


def format_treatment_label(treatment: str) -> str:
    """Wrap a workbook treatment name for a compact category axis."""
    normalized = treatment.lower()
    bin_width = "0.1" if normalized.startswith("0.1") else "1"
    if "frequency corrected" in normalized:
        return f"{bin_width} $m/z$ bin + FC"
    return f"{bin_width} $m/z$ bin"


def set_mass_limits(
    axis: plt.Axes,
    values: np.ndarray,
    theoretical_mass: float,
) -> None:
    """Set a local mass scale that includes all bars and the reference line."""
    all_values = np.append(values, theoretical_mass)
    lower_value = float(np.min(all_values))
    upper_value = float(np.max(all_values))
    span = upper_value - lower_value
    lower_padding = max(span * 0.18, 0.5)
    upper_padding = max(span * 0.40, 0.5)
    axis.set_ylim(lower_value - lower_padding, upper_value + upper_padding)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.yaxis.set_major_formatter(StrMethodFormatter("{x:,.1f}"))
    axis.yaxis.get_offset_text().set_visible(False)


def style_axis(axis: plt.Axes) -> None:
    """Apply the manuscript styling shared by the SI figures."""
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(AXIS_STROKE_WIDTH)
    axis.spines["bottom"].set_linewidth(AXIS_STROKE_WIDTH)
    axis.tick_params(direction="out", width=AXIS_STROKE_WIDTH, length=4)


def make_figure(treatments: list[str], values: np.ndarray) -> plt.Figure:
    """Build and return the three-panel treatment bar-chart figure."""
    figure, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=FIGURE_SIZE_INCHES,
        constrained_layout=True,
    )
    x_positions = np.arange(len(treatments))
    treatment_labels = [format_treatment_label(name) for name in treatments]

    for index, (axis, (protein, theoretical_mass, _)) in enumerate(
        zip(axes, PROTEINS)
    ):
        protein_values = values[:, index]
        axis.bar(
            x_positions,
            protein_values,
            width=0.72,
            color=BAR_COLOR,
            edgecolor=BAR_EDGE_COLOR,
            linewidth=AXIS_STROKE_WIDTH,
            zorder=2,
        )
        axis.axhline(
            theoretical_mass,
            color=THEORETICAL_COLOR,
            linestyle=(0, (4, 3)),
            linewidth=1.2,
            label="Theo. mass",
            zorder=3,
        )

        set_mass_limits(axis, protein_values, theoretical_mass)
        style_axis(axis)
        axis.set_title(
            f"{protein}\n({theoretical_mass:,.2f} Da)",
            fontsize=10,
            fontweight="bold",
            pad=8,
            x=0.40 if index == 2 else 0.5,
        )
        axis.set_xticks(x_positions, labels=treatment_labels)
        axis.tick_params(axis="x", labelrotation=35, labelsize=10)
        for label in axis.get_xticklabels():
            label.set_horizontalalignment("right")
            label.set_rotation_mode("anchor")
        if index == 0:
            axis.set_ylabel("Mass Centroid (Da)")

        if index == 0:
            axis.legend(
                loc="upper left",
                bbox_to_anchor=(0.03, 0.90),
                frameon=False,
                fontsize=10,
                handlelength=2.2,
            )

    return figure


def main() -> None:
    """Render a 6 by 3.5 inch, transparent 300 DPI PNG."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "savefig.transparent": True,
        }
    )

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    treatments, values = load_figure_data(SPREADSHEET)
    figure = make_figure(treatments, values)
    try:
        figure.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, transparent=True)
    finally:
        plt.close(figure)
    print(f"Saved {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
