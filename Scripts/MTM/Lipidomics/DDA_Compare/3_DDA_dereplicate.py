import os
import pandas as pd
import matplotlib as mpl

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

mpl.use("WxAgg")

topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Annika\IQ-X\Feb 2026 UntDDA\MSDial\Marked Low Score\RTP2"
os.chdir(topdir)

df = pd.read_csv("Combined_DDA_RT_Corrected.csv")

df = full_name_assignment(df, rttol=0.2, rtcol="Average Rt(min)")
df = get_number_of_datasets(df, namecol="Metabolite name")

subdfs = split_dfs(df)

# Write df to first sheet, then each subdf to separate sheets in the same Excel file
with pd.ExcelWriter("Combined_DDA_Full.xlsx") as writer:
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
