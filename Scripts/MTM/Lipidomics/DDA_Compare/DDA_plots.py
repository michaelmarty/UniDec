import os
from PublicScripts.Lipids.LipidFunctions import *
from venn import venn
import matplotlib.pyplot as plt
import matplotlib as mpl
from PublicScripts.Lipids.SkylineResultsPlotter import class_color_map

mpl.use("WxAgg")

def venn_diagram(df, sets, namecol="Metabolite name", title="", ax=None):
    vdata = {}
    for s in sets:
        vdata[s + " " + str(len(df[df["Dataset"] == s][namecol].unique()))] = set(df[df["Dataset"] == s][namecol])
        print("Length of set", s, ":", len(df[df["Dataset"] == s][namecol].unique()))

    # Print number of unique names across all sets
    unique_names = set()
    for s in sets:
        unique_names.update(df[df["Dataset"] == s][namecol].unique())
    print("Number of unique names across all sets:", len(unique_names))

    if ax is None:
        plt.figure(figsize=(8, 8))
        ax = plt.gca()
    venn(vdata, ax=ax)
    plt.title(title)

def count_number_of_datasets(df, namecol="Metabolite name", dataset_col="Dataset", class_col="Ontology"):
    # Drop duplicate rows based on dataset and class to avoid counting the same metabolite multiple times
    df = df.drop_duplicates(subset=[dataset_col, class_col, namecol])
    # Count the number of unique datasets each metabolite is found in
    dataset_counts = df.groupby(namecol)[dataset_col].nunique()
    # Addback into the dataframe
    df = df.merge(dataset_counts.rename("Dataset Count"), left_on=namecol, right_index=True)
    return df

def get_class_counts(df, name_col="Metabolite name", dataset_col="Dataset", class_col="Ontology"):
    datasets = sorted(df[dataset_col].unique())

    counts = []

    for i in range(2, len(datasets) + 1):
        subdf = df[df["Dataset Count"] == i]
        # Drop repeated names
        subdf = subdf.drop_duplicates(subset=[class_col, name_col])
        class_counts = subdf.groupby([dataset_col, class_col]).size().unstack(fill_value=0)
        # Sum across datasets to get total counts per class for coloring
        class_counts = class_counts.sum(axis=0)
        print(f"Number of unique names in {i} datasets:", subdf[name_col].nunique())
        cdf = pd.DataFrame(class_counts, columns=[str(i)]).transpose()
        counts.append(cdf)

    counts = counts[::-1]  # Reverse to have higher counts on bottom of stack
    class_counts_1 = df[df["Dataset Count"] == 1].groupby([dataset_col, class_col]).size().unstack(fill_value=0)
    counts.append(class_counts_1)
    # Combine counts for 1, 2, and 3 datasets
    class_counts = pd.concat(counts, axis=1).fillna(0)
    # Merge duplicate classes together by summing across them
    class_counts = class_counts.groupby(class_counts.columns, axis=1).sum()
    return class_counts

def bar_chart_of_classes(df, dataset_col="Dataset", class_col="Ontology", name_col="Simple name", title="", ax=None):
    # Drop duplicate rows based on dataset and class to avoid counting the same metabolite multiple times
    df = df.drop_duplicates(subset=[dataset_col, class_col, name_col])
    # Create counts of classes in all datasets
    # Print number of unique classes
    print("Number of unique classes:", df[class_col].nunique())
    df = count_number_of_datasets(df)
    class_counts = get_class_counts(df, name_col=name_col, dataset_col=dataset_col, class_col=class_col)

    if ax is None:
        plt.figure(figsize=(10, 6))
    else:
        plt.sca(ax)
    colors = [class_color_map.get(cls, "#333333") for cls in class_counts.columns]
    class_counts.plot(kind="bar", stacked=True, ax=ax, color=colors)

    plt.title(title)
    plt.ylabel("Count")
    plt.xlabel(dataset_col)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, reverse=True)
    return df


topdir = r"Z:\Group Share\Shivam\DDA"
os.chdir(topdir)

datasets = ["ITIT", "OTIT", "OTOT", "OTOT_AqX", "Stellar_RT_1", "Stellar"]

df = pd.read_csv("Combined_DDA.csv")

# Remove low quality ids where comment includes "Low quality"
df = df[~df["Comment"].str.contains("Low quality", case=False, na=False)]

simpnames = df["Metabolite name"].apply(lambda x: x.split("|")[0] if "|" in x else x)
df["Simple name"] = simpnames

# Simplify classes by applying hg_simplifier dict to Ontology column
df["Ontology"] = df["Ontology"].replace(hg_simplifier)

plt.figure(figsize=(18, 10))
plt.subplot(1, 2, 1)
venn_diagram(df, datasets, namecol="Simple name", title="Overlap of Identified Metabolites", ax=plt.gca())
plt.subplot(1, 2, 2)
df = bar_chart_of_classes(df, dataset_col="Dataset", class_col="Ontology", title="Class Distribution by Dataset", ax=plt.gca())
# Write out df with dataset counts for reference

# Drop duplicates based on class and name to get unique metabolites for each class
# df = df.drop_duplicates(subset=["Ontology", "Metabolite name"])
# Sort by ontology and metabolite name
df = df.sort_values(by=["Ontology", "Metabolite name", "Adduct type"])
df.to_csv("Combined_DDA_with_Dataset_Counts.csv", index=False)
plt.tight_layout()
plt.show()

