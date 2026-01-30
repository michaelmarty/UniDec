import os
import pandas as pd
import PublicScripts.Lipids.LipidFragPrediction as lfp

def rename_overlapping_species(df):
    # Find all species that share the same name and adduct but have different retention times, rename them as _A _B, etc
    grouped = df.groupby(["Molecule Name", "Precursor Adduct"])
    rename_map = {}
    for (name, adduct), group in grouped:
        if len(group) > 1:
            # Check if all retention times are different and continue if they are all the same
            rts = group["Explicit Retention Time"].values
            if len(set(rts)) == 1:
                continue

            # If not, rename all with the same retention time with suffixes
            rt_groups = group.groupby("Explicit Retention Time")
            suffixes = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            i=0
            for rt, rt_group in rt_groups:
                if len(rt_group) >= 1:
                    for j, idx in enumerate(rt_group.index):
                        new_name = f"{name}_{suffixes[i]}"
                        rename_map[idx] = new_name
                    i+=1

    df = df.copy()
    for idx, new_name in rename_map.items():
        df.at[idx, "Molecule Name"] = new_name
    return df

filepath = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Full Library Master Lists\Negative_updated_PeakMaster_1_2026_01_23_14_07_18_TransitionList.csv"
outfile = filepath.replace(".csv", "_NoPrecursorAssigned.csv")

df = pd.read_csv(filepath)

# Remove any where "Product Mz" is greater than 5 less than "Precursor Mz"
df_filtered = df[df["Product Mz"] < (df["Precursor Mz"] - 5)]

# Only take PC class for now
# df_filtered = df_filtered[df_filtered["Molecule List Name"] == "PC"]

df_filtered = rename_overlapping_species(df_filtered)

# Remove any where Precursor Adduct is [M-2H]
# df_filtered = df_filtered[df_filtered["Precursor Adduct"] != "[M-2H]"]

# Remove any where Class is PG and Product Mz is close to 227.0326 (common precursor fragment)
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PG") & (df_filtered["Product Mz"].between(227.0326 - 0.1, 227.0326 + 0.1)))]
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PG") & (df_filtered["Product Mz"].between(171.006 - 0.1, 171.006 + 0.1)))]

# Remove PE fragment near 152.9958
# df_filtered = df_filtered[~((df_filtered["Molecule List Name"] == "PE") & (df_filtered["Product Mz"].between(152.9958 - 0.1, 152.9958 + 0.1)))]

df_filtered = lfp.skyline_tail_namer(df_filtered)

# Drop any fragments where Product Name contains +HG
# df_filtered = df_filtered[~df_filtered["Product Name"].str.contains(r"\+HG")]
# df_filtered = df_filtered[~df_filtered["Product Name"].str.startswith("HG")]
# df_filtered = df_filtered[~df_filtered["Product Name"].str.startswith("NOT")]

# Drop molecular formula column if it exists
if "Precursor Mz" in df_filtered.columns:
    df_filtered = df_filtered.drop(columns=["Precursor Mz"])

# If charge is 0 set to 1
# df_filtered.loc[df_filtered["Precursor Charge"] == 0, "Precursor Charge"] = 1
# If charge is 0 drop those rows
df_filtered = df_filtered[df_filtered["Precursor Charge"] != 0]



df_filtered.to_csv(outfile, index=False)