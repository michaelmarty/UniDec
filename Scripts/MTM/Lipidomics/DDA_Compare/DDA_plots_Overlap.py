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

# # Drop any with class "BRSE"
# df = df[df["Ontology"] != "BRSE"]
df = df[df["Ontology"] != "HBMP"]

def plot_total_ids_by_dataset(df, datasets, dataset_col="Dataset"):
    """
    Makes a separate bar plot showing total IDs per dataset.
    Counts total rows, not unique metabolite names.
    """

    plot_df = df[df[dataset_col].isin(datasets)].copy()

    counts = (
        plot_df
        .groupby(dataset_col)
        .size()
        .reindex(datasets, fill_value=0)
    )

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
datasets = ["ITIT", "OTIT", "OTOT", "Stellar", "Stellarboth", "AX"]

plot_total_ids_by_dataset(df, datasets=datasets)

compare_classes_plot(df, datasets=datasets, use_simple_names=False, drop_low_quality=True)


