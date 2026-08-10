import os
from PublicScripts.Lipids.LipidFunctions import *
import matplotlib.pyplot as plt
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *


topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

df = pd.read_excel("Combined_DDA_OutlierAnalysis.xlsx", sheet_name="Sheet1")

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

print("\nRaw confirmed Stellar SM rows:")
print(len(stellar_sm))

print("\nUnique Stellar SM sum compositions:")
print(stellar_sm["Lipid Match Key"].nunique())

print("\nStellar SM names being counted:")
print(
    stellar_sm[
        ["Lipid Match Key", "Metabolite name"]
    ]
    .drop_duplicates(subset=["Lipid Match Key"])
    .sort_values("Lipid Match Key")
    .to_string(index=False)
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

#Figure 1, all IDs across all data
plot_total_ids_by_dataset(confirmed_df, datasets=datasets, bar_color="steelblue")


overlap_df = overlap_df.drop_duplicates(
    subset=[
        "Dataset",
        "Ontology",
        "Metabolite name"
    ]
).copy()

compare_classes_plot(overlap_df, datasets=datasets,  use_simple_names=False, drop_low_quality=True)

class_dataset_counts = plot_class_dataset_heatmap(df=confirmed_df, datasets=datasets)

# Load the full fragment-validation table
fragment_df = pd.read_excel("Combined_Fragment_Assignment_Colored.xlsx", sheet_name="Sheet1")

confirmed_counts, rejected_counts = (
    plot_confirmed_rejected_separate_bars(
        confirmed_df=confirmed_df,
        candidate_df=fragment_df,
        datasets=datasets,
        confirmed_color="#076E3D",
        rejected_color="#6E0739",
        title_size=14,
        label_size=14,
        tick_size=10,
        border_width=1
    )
)

scan_counts = {
    "OTOT":17477,
    "OTIT":61534,
    "ITIT": 63625,
    "Stellar": 187148,
    "Stellarboth": 70975

}
scan_summary = plot_scan_counts(
        scan_counts=scan_counts,
        bar_color="#035D99",
        title_size=14,
        label_size=14,
        tick_size=10,
        border_width=1

)

rejected_counts = plot_rejected_ids_by_dataset(candidate_df=fragment_df, datasets=datasets, bar_color="mediumvioletred")

# Display all figures at the end
plt.show(block=True)