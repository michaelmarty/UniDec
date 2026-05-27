import os
from PublicScripts.Lipids.LipidFunctions import *
import PublicScripts.Lipids.LipidOutlierAnalysis as ot
from venn import venn
import matplotlib.pyplot as plt
import matplotlib as mpl

def prep_datasets(topdir, datasets, min_count=None, drop_d7=True, write_output=True, dtype="msdial", rttol=0.2):
    os.chdir(topdir)
    dfs = []
    # Data sets
    for d in datasets:
        if dtype=="msdial":
            mergeddf, posdf, negdf = cleanup_and_merge_msdial(f"{d}_pos.csv", f"{d}_neg.csv", f"{d}_pos.xml",
                                                              f"{d}_neg.xml", rttol=rttol)
        elif dtype=="mzmine":
            mergeddf, posdf, negdf = cleanup_and_merge_mzmine(f"{d}_Pos_features_data_msdial.csv", f"{d}_Neg_features_data_msdial.csv", rttol=rttol)
        # Write out merged dataframes for reference
        # mergeddf.to_csv(f"{d}_Merged.csv", index=False)
        mergeddf["Dataset"] = d
        dfs.append(mergeddf)

    # Combine all datasets
    combineddf = pd.concat(dfs)
    # Sort by class, then by name
    combineddf = combineddf.sort_values(by=["Ontology", "Metabolite name"])

    if drop_d7:
        # Drop any with (d7) in name
        combineddf = combineddf[~combineddf["Metabolite name"].str.contains(r"\(d7\)", na=False)]

    # Drop any where the class has less than min_count total entries across all datasets, to focus on more common classes
    if min_count is not None:
        class_counts = combineddf["Ontology"].value_counts()
        classes_to_keep = class_counts[class_counts >= 2].index
        combineddf = combineddf[combineddf["Ontology"].isin(classes_to_keep)]

    if write_output:
        # Write out combined dataframe for reference
        combineddf.to_csv("Combined_DDA_Full1.csv", index=False)
    return combineddf

def full_name_assignment(df, rttol=0.1, rtcol="Average Rt(min)"):
    df = df.copy()
    for i, row in df.iterrows():
        name = row["Metabolite name"]
        if "|" in name:
            simple_name = name.split("|")[0]
        else:
            simple_name = name

        candidates = df[df["Metabolite name"].str.startswith(simple_name)]
        # Take only candidates not equal to the current row
        candidates = candidates[candidates["Metabolite name"] != name]
        if len(candidates) == 0:
            # print(f"Warning: No candidates found for {row['Metabolite name']}")
            continue
        for j, candidate in candidates.iterrows():
            # For cls and other partially named candidates, check if the candidate name is longer
            if len(candidate["Metabolite name"]) > len(name):
                # Check if RTs are within rttol
                if abs(candidate[rtcol] - row[rtcol]) < rttol:
                    df.at[i, "Metabolite name"] = candidate["Metabolite name"]
                    break
                else:
                    continue
            else:
                # print(f"Warning: Candidate {candidate['Metabolite name']} is not longer than {row['Metabolite name']}")
                # print(f"Warning: Candidate {candidate['Metabolite name']} does not match simple name {simple_name}")
                pass

    return df

def get_number_of_datasets(df, namecol="Metabolite name"):
    df = df.copy()
    for i, row in df.iterrows():
        name = row[namecol]
        count = len(df[df[namecol] == name]["Dataset"].unique())
        df.at[i, "Dataset Count"] = count
    return df

def split_dfs(df):
    max_num_datasets = int(df["Dataset Count"].max())
    dfs = {}
    # For everything with 2, 3, 4, etc. datasets, split into separate dfs
    for i in range(2, max_num_datasets + 1):
        subdf = df[df["Dataset Count"] == i]
        dfs[i] = subdf
    # For everything with only 1 dataset, split by dataset
    subdf = df[df["Dataset Count"] == 1]
    for dataset in subdf["Dataset"].unique():
        dataset_subdf = subdf[subdf["Dataset"] == dataset]
        dfs[dataset] = dataset_subdf
    return dfs

def write_dfs_to_excel(df, filename="Combined_DDA_Full.xlsx"):
    subdfs = split_dfs(df)

    # Write df to first sheet, then each subdf to separate sheets in the same Excel file
    with pd.ExcelWriter(filename) as writer:
        df.to_excel(writer, sheet_name="Full Data", index=False)
        for key, subdf in subdfs.items():
            if isinstance(key, int):
                sheet_name = f"{key} Datasets"
            else:
                sheet_name = f"{key} Only"
            subdf.to_excel(writer, sheet_name=sheet_name, index=False)
        # Split df by dataset and write each to a separate sheet
        for dataset in df["Dataset"].unique():
            dataset_subdf = df[df["Dataset"] == dataset]
            sheet_name = f"{dataset} Total"
            dataset_subdf.to_excel(writer, sheet_name=sheet_name, index=False)

def dataset_split(df, rttol=0.2, rtcol="Average Rt(min)", namecol="Metabolite name", write_output=True):
    df = full_name_assignment(df, rttol=rttol, rtcol=rtcol)
    df = get_number_of_datasets(df, namecol=namecol)
    if write_output:
        write_dfs_to_excel(df)
    return df

def outlier_setup_dda(df, mztol=0.3, write_output=True, do_heads=True, do_tails=False):
    df = smiles_analysis(df)

    posdf = df[df["Polarity"] == "Positive"]
    negdf = df[df["Polarity"] == "Negative"]

    outdf = ot.outlier_analysis(posdf, negdf, tol=mztol, contribs=False, add_all_cols_heads=True,
                                add_all_cols_tails=False, drop_cols=False, do_tails=do_tails, do_heads=do_heads)

    if write_output:
        # Write out combined dataframe for reference
        outdf.to_excel("Combined_DDA_OutlierAnalysis.xlsx", index=False)

    return outdf


def venn_diagram(df, sets, namecol="Metabolite name", title="", ax=None):
    vdata = {}
    for s in sets:
        subdf = df[df["Dataset"] == s]
        unique_names = subdf[namecol].unique()
        l = len(unique_names)
        vdata[s + " " + str(l)] = set(unique_names)

    if ax is None:
        plt.figure(figsize=(8, 8))
        ax = plt.gca()
    venn(vdata, ax=ax)
    plt.title(title)

def count_number_of_datasets(df, namecol="Metabolite name", dataset_col="Dataset", class_col="Ontology"):
    # If Dataset Count column already exists, drop it to avoid confusion
    if "Dataset Count" in df.columns:
        df = df.drop(columns=["Dataset Count"])
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
        cdf = pd.DataFrame(class_counts, columns=[str(i) + " Datasets"]).transpose()
        counts.append(cdf)

    counts = counts[::-1]  # Reverse to have higher counts on bottom of stack
    class_counts_1 = df[df["Dataset Count"] == 1].groupby([dataset_col, class_col]).size().unstack(fill_value=0)
    counts.append(class_counts_1)
    # Combine counts for 1, 2, and 3 datasets
    class_counts = pd.concat(counts, axis=1).fillna(0)
    # Merge duplicate classes together by summing across them
    class_counts = class_counts.T.groupby(class_counts.columns).sum().T
    # Rename rows that have a dataset in them to say dataset + only
    newrows = []
    for idx in class_counts.index:
        if any(dataset in idx for dataset in datasets):
            newrows.append(idx + " Only")
        else:
            newrows.append(idx)
    class_counts.index = newrows
    return class_counts

def bar_chart_of_classes(df, datasets, dataset_col="Dataset", class_col="Ontology", name_col="Metabolite name", title="", ax=None):
    class_counts = get_class_counts(df, name_col=name_col, dataset_col=dataset_col, class_col=class_col)
    # Sort order of class counts rows to match datasets order with only at the end, add the other classes at the beginning in the same order they appear in the dataframe
    new_order = []
    for idx in class_counts.index:
        if "Only" not in idx:
            new_order.append(idx)
    for dataset in datasets:
        only_row = dataset + " Only"
        if only_row in class_counts.index:
            new_order.append(only_row)
    class_counts = class_counts.reindex(new_order)

    if ax is None:
        plt.figure(figsize=(10, 6))
    else:
        plt.sca(ax)
    colors = [class_color_map.get(cls, "#333333") for cls in class_counts.columns]
    class_counts.plot(kind="bar", stacked=True, ax=ax, color=colors)
    # rotate x tick labels by 45 degrees
    plt.xticks(rotation=45, ha="right")

    plt.title(title)
    plt.ylabel("Count")
    # plt.xlabel(dataset_col)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, reverse=True)

def compare_classes_plot(df, datasets=None, use_simple_names=True, use_simple_classes=True, drop_low_quality=True):
    df = df.copy()
    if datasets is None:
        datasets = sorted(df["Dataset"].unique())

    if drop_low_quality:
        # Remove low quality ids where comment includes "Low quality"
        df = df[~df["Comment"].str.contains("Low quality", case=False, na=False)]

    if use_simple_names:
        simpnames = df["Metabolite name"].apply(lambda x: x.split("|")[0] if "|" in x else x)
        df["Simple name"] = simpnames
        namecol = "Simple name"
    else:
        namecol = "Metabolite name"

    # Drop duplicate rows based on dataset and class to avoid counting the same metabolite multiple times
    df = df.drop_duplicates(subset=["Dataset", namecol])

    # Simplify classes by applying hg_simplifier dict to Ontology column
    df["Ontology"] = df["Ontology"].replace(hg_simplifier)

    if use_simple_classes:
        # Simplify classes by applying class_simplifier dict to Ontology column
        df = simplify_class_df(df)

    # Count number of datasets for each metabolite name and add as a column to the dataframe
    df = count_number_of_datasets(df, namecol=namecol)
    # Rename datasets to add / to OTOT, OTIT, and ITIT
    df["Dataset"] = df["Dataset"].replace({"OTOT": "OT/OT", "OTIT": "OT/IT", "ITIT": "IT/IT"})
    datasets = [d.replace("OTOT", "OT/OT").replace("OTIT", "OT/IT").replace("ITIT", "IT/IT") for d in datasets]

    plt.figure(figsize=(18, 10))
    plt.subplot(1, 2, 1)
    venn_diagram(df, datasets, namecol=namecol, title="Overlap of Identified Metabolites", ax=plt.gca())
    plt.subplot(1, 2, 2)
    bar_chart_of_classes(df, datasets, dataset_col="Dataset", class_col="Ontology", title="Class Distribution by Dataset",
                         name_col=namecol, ax=plt.gca())
    plt.tight_layout()
    plt.savefig("DDA_Overlap_and_Class_Distribution.png", dpi=600, bbox_inches="tight")
    plt.savefig("DDA_Overlap_and_Class_Distribution.pdf", bbox_inches="tight")
    plt.show()