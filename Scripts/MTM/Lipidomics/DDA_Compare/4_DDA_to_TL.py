import pandas as pd
import os
from PublicScripts.Lipids.LipidFunctions import msdial_to_transitionlist, tl_to_il, derplicate_il

topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Annika\IQ-X\Feb 2026 UntDDA\MSDial\Marked Low Score\RTP2"
os.chdir(topdir)

filepath = r"Combined_DDA_Full.xlsx"
# Read in the combined DDA data with different sheets for each dataset
dfs = pd.read_excel(filepath, sheet_name=None)

tldfs = {}
ildfs = {}
for key, df in dfs.items():
    # print(f"Processing sheet: {key}")
    # Drop any that have "Low Quality" in the "Comment" column, to focus on good IDs
    df = df[~df["Comment"].str.contains("Low Quality", case=False, na=False)]
    tl_df = msdial_to_transitionlist(df)
    tldfs[key] = tl_df
    ildf = tl_to_il(tl_df)
    ildf = derplicate_il(ildf)
    ildfs[key] = ildf
    print("Length of isolation list for {}: {}".format(key, len(ildf)))

# Write out each transition list to a separate sheet in the same Excel file
with pd.ExcelWriter("Combined_DDA_TransitionLists.xlsx") as writer:
    for key, tl_df in tldfs.items():
        sheet_name = f"{key} TL"
        tl_df.to_excel(writer, sheet_name=sheet_name, index=False)


# Write out each isolation list to a separate sheet in the same Excel file
with pd.ExcelWriter("Combined_DDA_IsolationLists.xlsx") as writer:
    for key, ildf in ildfs.items():
        sheet_name = f"{key} IL"
        ildf.to_excel(writer, sheet_name=sheet_name, index=False)