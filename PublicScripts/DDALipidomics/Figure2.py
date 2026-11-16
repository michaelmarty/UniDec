"""Save Figure 2A's confirmed-lipid Venn diagram as a standalone figure.

Use the same workbook, dataset order, and sum-composition names as Figure 2.
Save Figure2.png and Figure2.pdf in Final Figures beside the workbook. Set MPLBACKEND=Agg
to render without an interactive window.
"""

import os
from pathlib import Path

import matplotlib

matplotlib.use(os.environ.get("MPLBACKEND", "WxAgg"))
import matplotlib.pyplot as plt
import pandas as pd

from Figure3 import DATASET_LABELS, datasets, input_file, prepare_overlap, topdir
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import venn_diagram


INPUT_FILE = Path(topdir) / input_file
OUTPUT_DIR = INPUT_FILE.parent / "Final Figures"


def plot_overlap(source):
    overlap = prepare_overlap(source)
    if overlap.empty:
        raise ValueError("No confirmed lipids in the five manuscript datasets.")
    fig, ax = plt.subplots(figsize=(3.25, 3.6))
    venn_diagram(
        overlap,
        sets=[DATASET_LABELS[name] for name in datasets],
        ax=ax,
        number_fontsize=8,
        legend_fontsize=8,
        legend_loc="upper center",
        legend_bbox_to_anchor=(0.5, 0.0),
        legend_ncol=2,
    )
    # The Venn patches fit within the unit square; remove the helper's margins.
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.tight_layout(pad=0.3)
    return fig


if __name__ == "__main__":
    source = pd.read_excel(INPUT_FILE, sheet_name="Sheet1")
    figure = plot_overlap(source)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf"]:
        output_path = OUTPUT_DIR / f"Figure2.{extension}"
        figure.savefig(output_path, dpi=600, bbox_inches="tight", pad_inches=0.04)
        print(f"Saved {output_path}")
    plt.show(block=True)
