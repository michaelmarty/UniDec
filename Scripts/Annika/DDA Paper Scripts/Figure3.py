import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

input_file = "Combined_Fragment_Assignment_Colored.xlsx"
output_basename = "Figure 3"
datasets = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]

# Muted sequential colors emphasize that coverage is ordered from one to five
# datasets. This is deliberately not viridis because viridis identifies
# individual datasets in the other figures.
coverage_colors = {
    1: "#D9EAF4",
    2: "#B8D8E8",
    3: "#8FC1DA",
    4: "#5FA6C8",
    5: "#347FA8",
}

# Collapse related subclasses in the same way as the Figure 2 class plot.
class_groups = {
    "Cer_AP": "Cer",
    "Cer_HDS": "Cer",
    "Cer_NDS": "Cer",
    "Cer_NP": "Cer",
    "Cer_NS": "Cer",
    "Cer_EODS": "Cer",
    "Cer_EOS": "Cer",
    "Cer_HS": "Cer",
    "Cer_AS": "Cer",
    "Cer_ADS": "Cer",
    "HexCer_AP": "HexCer",
    "HexCer_NDS": "HexCer",
    "HexCer_HS": "HexCer",
    "HexCer_NS": "HexCer",
    "HexCer_EOS": "HexCer",
    "HexCer_HDS": "HexCer",
    "AHexCer": "HexCer",
    "Hex2Cer": "HexCer",
    "Hex3Cer": "HexCer",
    "SHexCer": "HexCer",
    "ASM": "SM",
    "LNAPE": "PE",
    "LNAPS": "PS",
    "NAPE": "PE",
    "NAPS": "PS",
    "PE_Cer": "Cer",
    "PI_Cer": "Cer",
    "ADGGA": "GL",
    "GM3": "GL",
    "DGDG": "GL",
    "DGGA": "GL",
    "MGDG": "GL",
    "SQDG": "GL",
}


def simplify_lipid_class(value):
    lipid_class = str(value).strip()
    lipid_class = class_groups.get(lipid_class, lipid_class)
    lipid_class = lipid_class.replace("Ox", "").replace("Ether", "")
    if lipid_class.startswith("L"):
        lipid_class = lipid_class[1:]
    return class_groups.get(lipid_class, lipid_class)


df = pd.read_excel(input_file, sheet_name="Sheet1")

required_columns = {
    "Dataset",
    "Metabolite name",
    "Ontology",
    "Overall Match Status",
}
missing_columns = required_columns.difference(df.columns)
if missing_columns:
    raise ValueError(
        "The input workbook is missing required columns: "
        + ", ".join(sorted(missing_columns))
    )

df["Dataset"] = (
    df["Dataset"]
    .astype(str)
    .str.replace("\u00a0", " ", regex=False)
    .str.strip()
)
df["Ontology"] = (
    df["Ontology"]
    .astype(str)
    .str.replace("\u00a0", " ", regex=False)
    .str.strip()
)

true_values = [True, "TRUE", "True", "true", 1, "1"]
confirmed_df = df[
    df["Dataset"].isin(datasets)
    & df["Overall Match Status"].isin(true_values)
].copy()

# Use sum-composition names so coverage is comparable with Figures 1 and 2.
confirmed_df["Lipid Name"] = (
    confirmed_df["Metabolite name"]
    .astype(str)
    .str.split("|", n=1)
    .str[0]
    .str.strip()
)
confirmed_df["Lipid Class"] = confirmed_df["Ontology"].apply(
    simplify_lipid_class
)

# Count a lipid once per dataset, then assign every unique class/name pair to
# the number of datasets in which it was confirmed.
confirmed_df = confirmed_df.drop_duplicates(
    subset=["Dataset", "Lipid Class", "Lipid Name"]
)
unique_lipids = (
    confirmed_df
    .groupby(["Lipid Class", "Lipid Name"])["Dataset"]
    .nunique()
    .rename("Coverage")
    .reset_index()
)

coverage_groups = list(range(1, len(datasets) + 1))
coverage_counts = pd.crosstab(
    unique_lipids["Lipid Class"],
    unique_lipids["Coverage"],
).reindex(columns=coverage_groups, fill_value=0).sort_index()

class_totals = coverage_counts.sum(axis=1)
coverage_percentages = coverage_counts.div(class_totals, axis=0) * 100

fig, ax = plt.subplots(figsize=(6.5, 4.875))
y_positions = np.arange(len(coverage_percentages.index))
left = np.zeros(len(coverage_percentages.index))

for coverage in coverage_groups:
    values = coverage_percentages[coverage].to_numpy()
    ax.barh(
        y_positions,
        values,
        left=left,
        height=0.72,
        color=coverage_colors[coverage],
        edgecolor="white",
        linewidth=0.7,
        label=f"{coverage} Dataset{'s' if coverage != 1 else ''}",
    )
    left += values

ax.set_yticks(y_positions, labels=coverage_percentages.index)
ax.invert_yaxis()
ax.set_xlim(0, 100)
ax.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
ax.set_xlabel("Confirmed Lipids in Class (%)", labelpad=12)
ax.set_ylabel("Lipid Class")
ax.set_title("Dataset Coverage of Confirmed Lipids by Class")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(
    title="Dataset Coverage",
    title_fontsize=9,
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    frameon=False,
    ncol=1,
)

fig.tight_layout()
fig.savefig(
    f"{output_basename}.png",
    dpi=600,
)
fig.savefig(
    f"{output_basename}.pdf",
)

print("\nConfirmed lipids in each class:")
print(class_totals.to_string())
print("\nCoverage percentages within each lipid class:")
print(coverage_percentages.round(1).to_string())

plt.show(block=True)
