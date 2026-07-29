import os
from PublicScripts.Lipids.LipidFunctions import *
import matplotlib.pyplot as plt
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *


topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

df = pd.read_excel("Combined_DDA_OutlierAnalysis.xlsx", sheet_name="Sheet1")

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

datasets = ["ITIT", "OTIT", "OTOT", "Stellar", "Stellarboth"]

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


# Make a separate dataframe for overlap analysis
overlap_df = df.copy()

# Convert detailed names to sum-composition names
overlap_df["Metabolite name"] = (
    overlap_df["Metabolite name"]
    .astype(str)
    .str.split("|", n=1)
    .str[0]
    .str.strip()
)
plot_total_ids_by_dataset(df, datasets=datasets, bar_color="seagreen")


overlap_df = overlap_df.drop_duplicates(
    subset=[
        "Dataset",
        "Ontology",
        "Metabolite name"
    ]
).copy()

compare_classes_plot(overlap_df, datasets=datasets,  use_simple_names=False, drop_low_quality=True)

class_dataset_counts = plot_class_dataset_heatmap(df=df, datasets=datasets)

# Load the full fragment-validation table
fragment_df = pd.read_excel("Combined_Fragment_Assignment_Colored.xlsx", sheet_name="Sheet1")

# # Plot True versus False quality distributions
# quality_distributions = plot_true_false_quality_distributions(candidate_df=fragment_df)

# # Plot which quality variables best separate True and False IDs
# feature_separation = plot_quality_feature_separation(candidate_df=fragment_df)

fragment_df = pd.read_excel("Combined_Fragment_Assignment_Colored.xlsx", sheet_name="Sheet1")

rejected_counts = plot_rejected_ids_by_dataset(candidate_df=fragment_df, datasets=datasets, bar_color="darkmagenta")

# Display all figures at the end
plt.show(block=True)


# presence_results = export_lipid_dataset_presence(
#     df=df,
#     datasets=datasets,
#     lipid_classes=["DG"],
#     lipid_col="Metabolite name",
#     class_col="Ontology",
#     dataset_col="Dataset",
#     output_file="Dataset_Presence.xlsx"
# )