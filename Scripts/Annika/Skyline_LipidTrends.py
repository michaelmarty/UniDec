import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# ============================================================
# USER SETTINGS
# ============================================================

input_file = r"Z:\Group Share\Estevan\Lipid Stuffs\Molecule Transition Results.csv"

# Folder where the plots and processed tables will be saved.
output_folder = r"Z:\Group Share\Estevan\Lipid Stuffs\LipidTrends.csv"

# Set to True if you want replicate names shortened automatically.
shorten_replicate_names = True

color_map_file = r"Z:\Group Share\Annika\Lipid_ColorMap.xlsx"

color_df = pd.read_excel(color_map_file)

color_df["Lipid Class"] = color_df["Lipid Class"].astype(str).str.strip()
color_df["Color"] = color_df["Color"].astype(str).str.strip()

lipid_colors = dict(
    zip(color_df["Lipid Class"], color_df["Color"])
)

other_classes = [
    "LPA",
    "LPG",
    "LPI",
    "LPS",
    "MAG",
    "NAE",
    "DLCL",
    "MGDG",
    "LPE-N",
]

# ============================================================
# READ THE SKYLINE FILE
# ============================================================

output_folder = Path(output_folder)
output_folder.mkdir(parents=True, exist_ok=True)

# Automatically detect comma-separated versus tab-separated files.
df = pd.read_csv(input_file, sep=None, engine="python")

required_columns = [
    "Molecule",
    "Molecule List Name",
    "Replicate Name",
    "Fragment Ion",
    "Area",
]

missing_columns = [column for column in required_columns if column not in df.columns]

if missing_columns:
    raise ValueError(
        "The input file is missing these required columns:\n"
        + "\n".join(missing_columns)
    )


# Make sure Area is numeric.
df["Area"] = pd.to_numeric(df["Area"], errors="coerce").fillna(0)

# ============================================================
# MANUALLY ADD CHOLESTEROL PEAK AREAS
# ============================================================

cholesterol_data = pd.DataFrame({
    "Molecule": [
        "Cholesterol",
        "Cholesterol",
    ],
    "Molecule List Name": [
        "Cholesterol",
        "Cholesterol",
    ],
    "Replicate Name": [
        "20260723_AS_05_HEK_Discs_EC",
        "20260723_AS_06_HEKExtr_EC",
    ],
    "Fragment Ion": [
        "manual",
        "manual",
    ],
    "Area": [
        2194922.72,   # HEK Discs cholesterol area
        201383787.61,   # HEK Extract cholesterol area
    ],
})

df = pd.concat(
    [df, cholesterol_data],
    ignore_index=True
)
# ============================================================
# USE ALL TRANSITION AREAS
# ============================================================

# Keep precursor and fragment-ion transitions.
precursor_df = df.copy()
# ============================================================
# OPTIONAL: SHORTEN REPLICATE NAMES
# ============================================================

def shorten_name(name):
    """
    Removes a leading date, initials, and injection number when present.

    Example:
    20260723_AS_05_HEK_Discs_EC
    becomes:
    HEK_Discs_EC
    """

    name = str(name)

    shortened = re.sub(
        r"^\d{8}_[A-Za-z]+_\d+_",
        "",
        name
    )

    return shortened

precursor_df["Lipid Class"] = (
    precursor_df["Molecule List Name"]
    .astype(str)
    .str.strip()
)

precursor_df["Lipid Class"] = precursor_df["Lipid Class"].replace({
    "CerP": "Cer",
    "Cer_ADS": "Cer",
    "Cer_AP": "Cer",
    "Cer_EOS": "Cer",
    "Cer_NDS": "Cer",
    "Cer_NP": "Cer",
    "Cer_NS": "Cer",
    "Cer_AS": "Cer",
    "Cer_EODS": "Cer",

    "ASM": "SM",

    "AHexCer": "HexCer",
    "SHexCer": "HexCer",
    "Hex2Cer": "HexCer",
    "Hex3Cer": "HexCer",
    "HexCer_HDS": "HexCer",
    "HexCer_HS": "HexCer",
    "HexCer_NDS": "HexCer",
    "HexCer_NS": "HexCer",
})

other_classes = [
    "NAE",
    "LNAPS",
    "LNAPE",
    "HBMP",
    "Sph",
    "GM3",
    "DLCL",
    "MLCL",
    "MGDG",
    "LPE-N",
    "EtherDG",
    "EtherMGDG",
    "OxPE",
    "DHSph",
    "FAHFA"
]

precursor_df.loc[
    precursor_df["Lipid Class"].isin(other_classes),
    "Lipid Class"
] = "Other"
if shorten_replicate_names:
    precursor_df["Sample"] = precursor_df["Replicate Name"].apply(shorten_name)
else:
    precursor_df["Sample"] = precursor_df["Replicate Name"].astype(str)


# ============================================================
# EXTRACT LIPID SUM COMPOSITION
# ============================================================

def get_sum_composition(molecule):
    """
    Returns the lipid annotation before the vertical bar.

    Example:
    ADGGA 52:0|ADGGA 16:0_18:0_18:0
    becomes:
    ADGGA 52:0
    """

    return str(molecule).split("|")[0].strip()


def extract_chain_length(sum_composition):
    """
    Extracts the total carbon number from a sum composition.

    Examples:
    ADGGA 52:0      -> 52
    AHexCer 53:2;O3 -> 53
    PC O-36:2       -> 36
    """

    match = re.search(r"(?<!\d)(\d+):(\d+)", str(sum_composition))

    if match:
        return int(match.group(1))

    return pd.NA


def extract_unsaturation(sum_composition):
    """
    Extracts the total number of double bonds from a sum composition.

    Examples:
    ADGGA 52:0      -> 0
    AHexCer 53:2;O3 -> 2
    PC O-36:2       -> 2
    """

    match = re.search(r"(?<!\d)(\d+):(\d+)", str(sum_composition))

    if match:
        return int(match.group(2))

    return pd.NA


precursor_df["Sum Composition"] = precursor_df["Molecule"].apply(
    get_sum_composition
)

precursor_df["Chain Length"] = precursor_df["Sum Composition"].apply(
    extract_chain_length
)

precursor_df["Unsaturation"] = precursor_df["Sum Composition"].apply(
    extract_unsaturation
)


# ============================================================
# COLLAPSE DUPLICATE ROWS
# ============================================================

# This sums rows belonging to the same lipid annotation within a sample.
# Different molecular species listed as separate Molecule rows remain separate.
lipid_df = (
    precursor_df.groupby(
        [
            "Replicate Name",
            "Sample",
            "Molecule",
            "Sum Composition",
            "Lipid Class",
            "Chain Length",
            "Unsaturation",
        ],
        dropna=False,
        as_index=False,
    )["Area"]
    .sum()
)


# ============================================================
# CALCULATE mTIC AND NORMALIZE EACH LIPID
# ============================================================

# mTIC is the total precursor area in each individual replicate.
mtic_df = (
    lipid_df.groupby(
        ["Replicate Name", "Sample"],
        as_index=False
    )["Area"]
    .sum()
    .rename(columns={"Area": "mTIC"})
)

lipid_df = lipid_df.merge(
    mtic_df,
    on=["Replicate Name", "Sample"],
    how="left"
)

if (lipid_df["mTIC"] <= 0).any():
    bad_samples = lipid_df.loc[
        lipid_df["mTIC"] <= 0,
        "Replicate Name"
    ].unique()

    raise ValueError(
        "The following replicates have an mTIC of zero:\n"
        + "\n".join(map(str, bad_samples))
    )

# Fraction of the replicate mTIC.
lipid_df["mTIC Normalized Area"] = (
    lipid_df["Area"] / lipid_df["mTIC"]
)

# Percentage of the replicate mTIC.
lipid_df["Percent of mTIC"] = (
    lipid_df["mTIC Normalized Area"] * 100
)

# ============================================================
# CHOLESTEROL VS POPC AREA PERCENT
# ============================================================

chol_popc_df = lipid_df[
    (lipid_df["Molecule"].str.startswith("ST 27:1", na=False))
    |
    (lipid_df["Molecule"].str.startswith("PC 34:1", na=False))
].copy()

chol_popc_df["Component"] = "POPC"

chol_popc_df.loc[
    chol_popc_df["Molecule"].str.startswith("ST 27:1", na=False),
    "Component"
] = "Cholesterol"

chol_popc_summary = (
    chol_popc_df.groupby(
        ["Replicate Name", "Sample", "Component"],
        as_index=False
    )["Area"]
    .sum()
)

chol_popc_summary["Total Area"] = (
    chol_popc_summary.groupby(
        ["Replicate Name", "Sample"]
    )["Area"]
    .transform("sum")
)

chol_popc_summary["Area Percent"] = (
    chol_popc_summary["Area"]
    / chol_popc_summary["Total Area"]
    * 100
)
# ============================================================
# CREATE SUMMARY TABLES
# ============================================================

def make_percent_summary(data, category_column):
    """
    Sums the normalized lipid percentages for each category within
    each replicate.
    """

    summary = (
        data.dropna(subset=[category_column])
        .groupby(
            ["Replicate Name", "Sample", category_column],
            as_index=False
        )["Percent of mTIC"]
        .sum()
    )

    # Recalculate within each replicate to protect against minor
    # rounding differences or lipids lacking a parsed category.
    category_total = summary.groupby(
        "Replicate Name"
    )["Percent of mTIC"].transform("sum")

    summary["Percent"] = (
        summary["Percent of mTIC"] / category_total * 100
    )

    return summary


class_summary = make_percent_summary(
    lipid_df[
        lipid_df["Lipid Class"] != "Cholesterol"
    ],
    "Lipid Class"
)

unsaturation_summary = make_percent_summary(
    lipid_df[
        lipid_df["Unsaturation"].between(0, 6, inclusive="both")
    ],
    "Unsaturation"
)

chain_length_summary = make_percent_summary(
    lipid_df[
        lipid_df["Chain Length"].between(24, 60, inclusive="both")
    ],
    "Chain Length"
)

cholesterol_summary = lipid_df[
    lipid_df["Lipid Class"] == "Cholesterol"
][
    [
        "Replicate Name",
        "Sample",
        "Area",
        "mTIC",
        "Percent of mTIC",
    ]
].copy()

# ============================================================
# STACKED BAR PLOTTING FUNCTION
# ============================================================

def plot_stacked_percent(
    summary,
    category_column,
    title,
    legend_title,
    output_name,
):
    plot_df = summary.pivot_table(
        index="Sample",
        columns=category_column,
        values="Percent",
        aggfunc="sum",
        fill_value=0,
    )

    # Sort categories.
    if category_column in ["Unsaturation", "Chain Length"]:
        plot_df = plot_df.reindex(
            sorted(plot_df.columns),
            axis=1
        )
    else:
        plot_df = plot_df.reindex(
            sorted(plot_df.columns, key=str.lower),
            axis=1
        )

    # Lipid class stays as a stacked bar chart.
    if category_column == "Lipid Class":

        plot_colors = [
            lipid_colors.get(str(lipid_class).strip(), "#B3B3B3")
            for lipid_class in plot_df.columns
        ]

        figure_width = max(8, len(plot_df.index) * 1.2)

        ax = plot_df.plot(
            kind="bar",
            stacked=True,
            figsize=(figure_width, 7),
            width=0.75,
            edgecolor="black",
            linewidth=0.3,
            color=plot_colors,
        )

        ax.set_xlabel("")
        ax.set_ylabel("Percent of total mTIC")
        ax.set_ylim(0, 100)

        ax.tick_params(
            axis="x",
            rotation=45
        )

        ax.legend(
            title=legend_title,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            ncol=2,
            frameon=False,
        )

    # Unsaturation and chain length become line plots.
    else:

        ax = plot_df.T.plot(
            kind="line",
            marker="o",
            figsize=(8, 6),
            linewidth=2,
            color=["lightseagreen", "slateblue"],
        )

        ax.set_xlabel(legend_title)
        ax.set_ylabel("Percent")
        ax.set_ylim(bottom=0)

        ax.legend(
            title="Sample",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
        )

    ax.set_title(title, fontsize=15)

    plt.tight_layout()
    plt.show()
    plt.close()


# ============================================================
# GENERATE THE THREE PLOTS
# ============================================================
def plot_cholesterol(cholesterol_summary):

    ax = cholesterol_summary.plot(
        x="Sample",
        y="Percent of mTIC",
        kind="bar",
        figsize=(6, 6),
        width=0.65,
        color=lipid_colors.get("Cholesterol", "#B40808"),
        edgecolor="black",
        linewidth=0.5,
        legend=False,
    )

    ax.set_title("Cholesterol Relative to Total mTIC", fontsize=15)
    ax.set_xlabel("")
    ax.set_ylabel("Cholesterol (% of total mTIC)")
    ax.tick_params(axis="x", rotation=45)

    # Add percentage labels above the bars.
    for container in ax.containers:
        ax.bar_label(
            container,
            fmt="%.2f%%",
            padding=3
        )

    plt.tight_layout()
    plt.show()
    plt.close()

plot_stacked_percent(
    summary=class_summary,
    category_column="Lipid Class",
    title="Lipid Class Composition",
    legend_title="Lipid Class",
    output_name="01_lipid_class_percent",
)

plot_stacked_percent(
    summary=unsaturation_summary,
    category_column="Unsaturation",
    title="Lipid Unsaturation Distribution",
    legend_title="Double Bonds",
    output_name="02_unsaturation_percent",
)

plot_stacked_percent(
    summary=chain_length_summary,
    category_column="Chain Length",
    title="Lipid Chain-Length Distribution",
    legend_title="Total Carbons",
    output_name="03_chain_length_percent",
)

# ============================================================
# PLOT CHOLESTEROL VS POPC
# ============================================================

chol_popc_plot = chol_popc_summary.pivot_table(
    index="Sample",
    columns="Component",
    values="Area Percent",
    fill_value=0
)

chol_popc_plot = chol_popc_plot.reindex(
    columns=["Cholesterol", "POPC"],
    fill_value=0
)

ax = chol_popc_plot.plot(
    kind="bar",
    stacked=True,
    figsize=(6, 6),
    width=0.65,
    color=[
        lipid_colors.get("ST", "#B40808"),
        lipid_colors.get("PC", "#788BFF"),
    ],
    edgecolor="black",
    linewidth=0.5,
)

ax.set_title("Cholesterol vs POPC Area Composition", fontsize=15)
ax.set_xlabel("")
ax.set_ylabel("Area Percent")
ax.set_ylim(0, 100)
ax.tick_params(axis="x", rotation=45)

ax.legend(
    title="Component",
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False
)

plt.tight_layout()
plt.show()
plt.close()
# ============================================================
# SAVE PROCESSED DATA
# ============================================================

lipid_df.to_csv(
    output_folder / "mTIC_normalized_lipids.csv",
    index=False
)

mtic_df.to_csv(
    output_folder / "mTIC_by_replicate.csv",
    index=False
)

class_summary.to_csv(
    output_folder / "percent_by_lipid_class.csv",
    index=False
)

unsaturation_summary.to_csv(
    output_folder / "percent_by_unsaturation.csv",
    index=False
)

chain_length_summary.to_csv(
    output_folder / "percent_by_chain_length.csv",
    index=False
)


# ============================================================
# PRINT BASIC CHECKS
# ============================================================

print("\nCalculated mTIC values:")
print(mtic_df.to_string(index=False))

print("\nPercent total by sample and plot type:")

print("\nLipid class totals:")
print(
    class_summary.groupby("Sample")["Percent"]
    .sum()
    .round(5)
)

print("\nUnsaturation totals:")
print(
    unsaturation_summary.groupby("Sample")["Percent"]
    .sum()
    .round(5)
)

print("\nChain-length totals:")
print(
    chain_length_summary.groupby("Sample")["Percent"]
    .sum()
    .round(5)
)

unparsed = lipid_df[
    lipid_df["Chain Length"].isna()
    | lipid_df["Unsaturation"].isna()
]

if not unparsed.empty:
    print(
        "\nWARNING: The following sum compositions could not be parsed "
        "for chain length or unsaturation:"
    )

    print(
        unparsed[
            ["Lipid Class", "Sum Composition"]
        ]
        .drop_duplicates()
        .to_string(index=False)
    )

# print(f"\nFinished. Results were saved to:\n{output_folder}")