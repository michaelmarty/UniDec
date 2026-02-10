import os
from PublicScripts.Lipids.LipidFunctions import *

topdir = r"Z:\Group Share\Shivam\DDA"
os.chdir(topdir)

datasets = ["ITIT", "OTIT", "OTOT", "OTOT_AqX", "Stellar_RT_1", "Stellar"]

dfs = []
# Data sets
for d in datasets:
    mergeddf, posdf, negdf = cleanup_and_merge_msdial(f"{d}_pos.csv", f"{d}_neg.csv", f"{d}_pos.xml", f"{d}_neg.xml", rttol=0.2)
    # Write out merged dataframes for reference
    # mergeddf.to_csv(f"{d}_Merged.csv", index=False)
    mergeddf["Dataset"] = d
    dfs.append(mergeddf)

# Combine all datasets
combineddf = pd.concat(dfs)
# Sort by class, then by name
combineddf = combineddf.sort_values(by=["Ontology", "Metabolite name"])

# Drop any with (d7) in name
combineddf = combineddf[~combineddf["Metabolite name"].str.contains(r"\(d7\)", na=False)]

# Drop any where the class has less than 2 total entries across all datasets, to focus on more common classes
# class_counts = combineddf["Ontology"].value_counts()
# classes_to_keep = class_counts[class_counts >= 2].index
# combineddf = combineddf[combineddf["Ontology"].isin(classes_to_keep)]

# Write out combined dataframe for reference
combineddf.to_csv("Combined_DDA.csv", index=False)