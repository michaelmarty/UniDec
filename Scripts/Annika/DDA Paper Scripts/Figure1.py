import os
from PublicScripts.Lipids.LipidFunctions import *
import matplotlib.pyplot as plt
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *


plt.rcParams.update(
    {
        "font.family": "Arial",
        "axes.titlesize": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
    }
)


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


overlap_df = overlap_df.drop_duplicates(
    subset=[
        "Dataset",
        "Ontology",
        "Metabolite name"
    ]
).copy()


# Load the full fragment-validation table
fragment_df = pd.read_excel("Combined_Fragment_Assignment_Colored.xlsx", sheet_name="Sheet1")

dataset_display_names = {"Stellarboth": "Stellar (+/-)"}
confirmed_df["Dataset"] = confirmed_df["Dataset"].replace(dataset_display_names)
fragment_df["Dataset"] = fragment_df["Dataset"].replace(dataset_display_names)
plot_datasets = [dataset_display_names.get(name, name) for name in datasets]

confirmed_counts, rejected_counts = (
    plot_confirmed_rejected_separate_bars(
        confirmed_df=confirmed_df,
        candidate_df=fragment_df,
        datasets=plot_datasets,
        output_basename="Figure1",
        confirmed_color="#076E3D",
        rejected_color="#6E0739",
        title_size=16,
        label_size=16,
        tick_size=14,
        border_width=1
    )
)
# Display all figures at the end
plt.show(block=True)
