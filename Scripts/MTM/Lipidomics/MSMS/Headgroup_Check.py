import os
import pandas as pd
import numpy as np

from PublicScripts.Lipids.LipidFunctions import check_headgroup_fragments



os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

fragmentslib = "LipidHeadgroupFragments.csv"
posfile = "PosIDsRTP2.csv"
negfile = "NegIDsRTP2.csv"

fragdf = pd.read_csv(fragmentslib)

posdf = pd.read_csv(posfile)
negdf = pd.read_csv(negfile)

# Sort posdf and negdf by Ontology and Adduct type
posdf = posdf.sort_values(by=["Ontology", "Adduct type"]).reset_index(drop=True)
negdf = negdf.sort_values(by=["Ontology", "Adduct type"]).reset_index(drop=True)

posdf = check_headgroup_fragments(posdf, fragdf)
negdf = check_headgroup_fragments(negdf, fragdf)

# write out to modified csv
posdf.to_csv("PosIDsRTP2_HeadgroupCheck.csv", index=False)
negdf.to_csv("NegIDsRTP2_HeadgroupCheck.csv", index=False)


import matplotlib.pyplot as plt

plt.figure(figsize=(16,12))

# add categorical column for NL match
# posdf["HG_NL_Category"] = posdf["HG_NL_Match"].apply(categorize_nl)
# negdf["HG_NL_Category"] = negdf["HG_NL_Match"].apply(categorize_nl)

# Set viridis colormap
cm = plt.get_cmap("viridis_r")


# Positive: HG_Fragment_Match
plt.subplot(2,2,1)
posdf_frag = posdf[posdf["HG_Fragment_Match"] != -1]
pos_group_frag = posdf_frag.groupby(["Ontology", "Adduct type", "HG_Fragment_Match"]).size().unstack(fill_value=0)
pos_group_frag.plot(kind="bar", stacked=True, ax=plt.gca(), legend=True, color=cm(np.linspace(0,1,pos_group_frag.shape[1])))
plt.title("Positive Mode HG Fragment Match by Class and Adduct")
plt.xlabel("Compound Class / Adduct")
plt.ylabel("Count")
plt.xticks(rotation=45, ha="right")

# Positive: HG_NL_Match
plt.subplot(2,2,2)
posdf_nl = posdf[posdf["HG_NL_Match"] != -1]
pos_group_nl = posdf_nl.groupby(["Ontology", "Adduct type", "HG_NL_Match"]).size().unstack(fill_value=0)
pos_group_nl.plot(kind="bar", stacked=True, ax=plt.gca(), legend=True, color=cm(np.linspace(0,1,pos_group_nl.shape[1])))
plt.title("Positive Mode HG NL Match by Class and Adduct")
plt.xlabel("Compound Class / Adduct")
plt.ylabel("Count")
plt.xticks(rotation=45, ha="right")

# Negative: HG_Fragment_Match
plt.subplot(2,2,3)
negdf_frag = negdf[negdf["HG_Fragment_Match"] != -1]

neg_group_frag = negdf_frag.groupby(["Ontology", "Adduct type", "HG_Fragment_Match"]).size().unstack(fill_value=0)
neg_group_frag.plot(kind="bar", stacked=True, ax=plt.gca(), legend=True, color=cm(np.linspace(0,1,neg_group_frag.shape[1])))
plt.title("Negative Mode HG Fragment Match by Class and Adduct")
plt.xlabel("Compound Class / Adduct")
plt.ylabel("Count")
plt.xticks(rotation=45, ha="right")

# Negative: HG_NL_Match
plt.subplot(2,2,4)
negdf_nl = negdf[negdf["HG_NL_Match"] != -1]
print(negdf_nl)
neg_group_nl = negdf_nl.groupby(["Ontology", "Adduct type", "HG_NL_Match"]).size().unstack(fill_value=0)
neg_group_nl.plot(kind="bar", stacked=True, ax=plt.gca(), legend=True, color=cm(np.linspace(0,1,neg_group_nl.shape[1])))
plt.title("Negative Mode HG NL Match by Class and Adduct")
plt.xlabel("Compound Class / Adduct")
plt.ylabel("Count")
plt.xticks(rotation=45, ha="right")

plt.tight_layout()
plt.show()



