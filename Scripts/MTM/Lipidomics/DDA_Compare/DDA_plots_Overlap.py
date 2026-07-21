import os
from PublicScripts.Lipids.LipidFunctions import *
from venn import venn
import matplotlib.pyplot as plt
import matplotlib as mpl
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *

mpl.use("WxAgg")

topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

df = pd.read_excel("Combined_DDA_Full.xlsx", sheet_name="Full Data")
df = pd.read_excel("Combined_DDA_OutlierAnalysis_Results.xlsx", sheet_name="All Data")
# df = df[df["Outlier_Prediction"] == 0]
df = pd.read_excel("Combined_DDA_OutlierAnalysis_Results.xlsx", sheet_name="All Data")


# # Drop any with class "BRSE"
# df = df[df["Ontology"] != "BRSE"]
df = df[df["Ontology"] != "HBMP"]


def plot_total_ids_by_dataset(df, datasets, dataset_col="Dataset"):
    """
    Makes a separate bar plot showing total IDs per dataset.
    Counts total rows, not unique metabolite names.
    """

    plot_df = df.copy()

    plot_df[dataset_col] = (
        plot_df[dataset_col]
        .astype(str)
        .str.replace("\u00a0", " ", regex=False)
        .str.strip()
    )

    print("\nDataset counts available to plotting function:")
    print(plot_df[dataset_col].value_counts(dropna=False).to_string())

    counts = (
        plot_df[plot_df[dataset_col].isin(datasets)]
        .groupby(dataset_col)
        .size()
        .reindex(datasets, fill_value=0)
    )

    print("\nCounts being plotted:")
    print(counts.to_string())

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.bar(counts.index, counts.values, edgecolor="black")

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Total IDs")
    ax.set_title("Total Lipid IDs by Dataset")

    for i, value in enumerate(counts.values):
        ax.text(
            i,
            value,
            str(int(value)),
            ha="center",
            va="bottom"
        )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()


def export_lipid_dataset_presence(
    df,
    datasets,
    lipid_classes,
    lipid_col="Lipid",
    class_col="Ontology",
    dataset_col="Dataset",
    output_file="Lipid_Dataset_Presence.xlsx"
):
    """
    Exports dataset presence/absence information for specified lipid classes.

    Output includes:
        - One row per unique lipid
        - 1/0 presence columns for each dataset
        - Number of datasets containing the lipid
        - Number of datasets missing the lipid
        - Names of datasets containing the lipid
        - Names of datasets missing the lipid
        - Separate Excel sheets grouped by number of datasets present
    """

    presence_df = df.copy()

    # Clean relevant text columns
    for col in [lipid_col, class_col, dataset_col]:
        presence_df[col] = (
            presence_df[col]
            .astype(str)
            .str.replace("\u00a0", " ", regex=False)
            .str.strip()
        )

    # Keep only requested datasets and lipid classes
    presence_df = presence_df[
        presence_df[dataset_col].isin(datasets)
        & presence_df[class_col].isin(lipid_classes)
    ].copy()

    # Remove rows without a meaningful lipid name
    presence_df = presence_df[
        presence_df[lipid_col].notna()
        & ~presence_df[lipid_col].isin(["", "nan", "None"])
    ]

    # One row per lipid/class/dataset before making the table
    presence_df = presence_df[
        [class_col, lipid_col, dataset_col]
    ].drop_duplicates()

    # Build presence/absence matrix
    presence_matrix = (
        presence_df
        .assign(Present=1)
        .pivot_table(
            index=[class_col, lipid_col],
            columns=dataset_col,
            values="Present",
            aggfunc="max",
            fill_value=0
        )
        .reindex(columns=datasets, fill_value=0)
        .reset_index()
    )

    # Count datasets present and missing
    presence_matrix["Datasets Present"] = (
        presence_matrix[datasets].sum(axis=1)
    )

    presence_matrix["Datasets Missing"] = (
        len(datasets) - presence_matrix["Datasets Present"]
    )

    # List dataset names
    presence_matrix["Present In"] = presence_matrix.apply(
        lambda row: ", ".join(
            dataset for dataset in datasets if row[dataset] == 1
        ),
        axis=1
    )

    presence_matrix["Missing From"] = presence_matrix.apply(
        lambda row: ", ".join(
            dataset for dataset in datasets if row[dataset] == 0
        ),
        axis=1
    )

    # Sort by class, then by greatest coverage, then lipid name
    presence_matrix = presence_matrix.sort_values(
        by=[class_col, "Datasets Present", lipid_col],
        ascending=[True, False, True]
    )

    # Summary by class and number of datasets present
    summary = (
        presence_matrix
        .groupby([class_col, "Datasets Present"])
        .size()
        .reset_index(name="Number of Lipids")
        .sort_values([class_col, "Datasets Present"])
    )

    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:

        # Complete table
        presence_matrix.to_excel(
            writer,
            sheet_name="All Lipids",
            index=False
        )

        # Summary counts
        summary.to_excel(
            writer,
            sheet_name="Summary",
            index=False
        )

        # Separate sheet for each presence count
        for number_present in range(len(datasets) + 1):
            subset = presence_matrix[
                presence_matrix["Datasets Present"] == number_present
            ]

            if not subset.empty:
                subset.to_excel(
                    writer,
                    sheet_name=f"Present in {number_present}",
                    index=False
                )

        # Separate sheets for each requested class
        for lipid_class in lipid_classes:
            class_subset = presence_matrix[
                presence_matrix[class_col] == lipid_class
            ]

            if not class_subset.empty:
                safe_sheet_name = str(lipid_class)[:31]

                class_subset.to_excel(
                    writer,
                    sheet_name=safe_sheet_name,
                    index=False
                )

    print(f"\nPresence/absence spreadsheet saved as:")
    print(os.path.abspath(output_file))

    print("\nNumber of lipids by class and dataset coverage:")
    print(summary.to_string(index=False))

    return presence_matrix



    print("\nCounts being plotted:")
    print(counts.to_string())

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.bar(counts.index, counts.values, edgecolor="black")

    ax.set_xlabel("Dataset")
    ax.set_ylabel("Total IDs")
    ax.set_title("Total Lipid IDs by Dataset")

    for i, value in enumerate(counts.values):
        ax.text(
            i,
            value,
            str(int(value)),
            ha="center",
            va="bottom"
        )

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.show()


# datasets = ["ITIT", "OTIT", "OTOT", "Stellar", "Stellarboth", "AX", "Apex_AX"]
datasets = ["ITIT", "OTIT", "OTOT", "Stellar", "Stellarboth"]

plot_total_ids_by_dataset(df, datasets=datasets)

compare_classes_plot(df, datasets=datasets, use_simple_names=False, drop_low_quality=True)

# presence_results = export_lipid_dataset_presence(
#     df=df,
#     datasets=datasets,
#     lipid_classes=["PC", "PE"],
#     lipid_col="Metabolite name",
#     class_col="Ontology",
#     dataset_col="Dataset",
#     output_file="PC_PE_Dataset_Presence.xlsx"
# )


