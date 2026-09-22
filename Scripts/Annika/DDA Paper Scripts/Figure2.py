import os
from PublicScripts.Lipids.LipidFunctions import *
import matplotlib.pyplot as plt
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
    }
)


topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

# Use the final comparison workbook so accepted IDs are defined directly by
# Overall Match Status rather than inferred from a prefiltered export.
df = pd.read_excel("Combined_Fragment_Assignment_Colored.xlsx", sheet_name="Sheet1")

stellar_sm = df[
    (df["Dataset"] == "Stellar")
    & (df["Ontology"] == "SM")
    & (
        df["Overall Match Status"].isin(
            [True, "TRUE", "True", "true", 1, "1"]
        )
    )
].copy()

stellar_sm["Lipid Match Key"] = (
    stellar_sm["Metabolite name"]
    .astype(str)
    .str.split("|", n=1)
    .str[0]
    .str.strip()
)

removed_df = pd.read_excel("Removed_Identifications.xlsx")

# # Drop any with class "BRSE"
# df = df[df["Ontology"] != "BRSE"]
# df = df[df["Ontology"] != "HBMP"]

df["Overlap Name"] = (
    df["Metabolite name"]
    .astype(str)
    .str.split("|")
    .str[0]
    .str.strip()
)


datasets = ["OTOT", "OTIT", "ITIT", "Stellar", "Stellarboth"]

# Clean dataset and class names
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

# Keep only the five selected datasets
df = df[df["Dataset"].isin(datasets)].copy()


confirmed_df = df[
    df["Overall Match Status"].isin(
        [True, "TRUE", "True", "true", 1, "1"]
    )
].copy()

# Make the overlap analysis use the same confirmed-ID population as Figure 1
overlap_df = confirmed_df.copy()

# Convert detailed names to sum-composition names
overlap_df["Metabolite name"] = (
    overlap_df["Metabolite name"]
    .astype(str)
    .str.split("|", n=1)
    .str[0]
    .str.strip()
)




overlap_df = overlap_df.drop_duplicates(
    subset=[
        "Dataset",
        "Ontology",
        "Metabolite name"
    ]
).copy()

dataset_display_names = {
    "OTOT": "OT/OT",
    "OTIT": "OT/IT",
    "ITIT": "IT/IT",
    "Stellarboth": "Stellar (+/-)",
}
overlap_df["Dataset"] = overlap_df["Dataset"].replace(dataset_display_names)
plot_datasets = [dataset_display_names.get(name, name) for name in datasets]

compare_classes_plot(
    overlap_df,
    datasets=plot_datasets,
    use_simple_names=False,
    drop_low_quality=False,
    output_basename="Figure 2",
    figure_size=(6.5, 3.61),
    save_bbox_inches=None,
    overlap_style="venn",
    class_legend_fontsize=8,
    class_legend_bbox_to_anchor=(1.02, 1.0),
    class_legend_loc="upper left",
    hide_class_top_and_right=True,
    rounded_six_set_venn=True,
    venn_number_fontsize=9,
    venn_legend_fontsize=8,
    venn_legend_loc="upper center",
    venn_legend_bbox_to_anchor=(0.5, -0.03),
    venn_legend_ncol=2
)


# Load the full fragment-validation table
fragment_df = pd.read_excel("Combined_Fragment_Assignment_Colored.xlsx", sheet_name="Sheet1")



# Display all figures at the end
plt.show(block=True)
