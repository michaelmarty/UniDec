import os
import pandas as pd
import numpy as np

def categorize_nl(x):
    if pd.isna(x) or x == -1:
        return "No NL"
    if x == 1:
        return "All"
    if x == 0:
        return "None"
    return "Partial"

def check_headgroup_fragments(df, fragdf, tol=0.01, intthresh=0.01):
    good = []
    goodnl = []
    for index, row in df.iterrows():
        class_name = row["Ontology"]
        adduct = row["Adduct type"]
        if "+" in adduct:
            ionmode = "positive"
        else:
            ionmode = "negative"

        fragments = fragdf[
            (fragdf["CompoundClass"] == class_name) &
            (fragdf["Adduct"] == adduct) &
            (fragdf["IonMode"] == ionmode)
        ]
        if fragments.empty:
            good.append(-1)
            goodnl.append(-1)
            # print("Not Found:", class_name, adduct, ionmode)
            continue

        msms_mz = row["MS/MS spectrum"]
        msms_mz = np.array([m.split(":") for m in msms_mz.split(" ")]).astype(float)
        b1 = msms_mz[:, 1] >= intthresh * np.max(msms_mz[:, 1])
        msms_mz = msms_mz[b1, :]

        prodmz = fragments["HeadgroupFragment_mz"].values[0]
        try:
            prodmz = np.array([prodmz]).astype(float)
        except:
            prodmz = np.array([float(f) for f in prodmz.split(",") if f])

        if prodmz is None or np.any(np.isnan(prodmz)) or len(prodmz) == 0:
            good.append(-1)
            # print("No Fragments:", class_name, adduct, ionmode)
        else:
            match_count = sum(1 for mz in msms_mz[:,0] if any(abs(mz - f) < tol for f in prodmz))
            good.append(match_count/len(prodmz))
            # if match_count == len(prodmz):
            #     good.append(1)
            # else:
            #     good.append(0)

        nlmz = fragments["HeadgroupNL_mz"].values[0]
        try:
            nlmz = np.array([nlmz]).astype(float)
        except:
            nlmz = np.array([float(f) for f in nlmz.split(",") if f])

        if nlmz is None or np.any(np.isnan(nlmz)) or len(nlmz) == 0:
            goodnl.append(-1)
            # print("No NL Fragments:", class_name, adduct, ionmode)
        else:
            nl_mz = row["Average Mz"] - msms_mz[:,0]
            match_count_nl = sum(1 for mz in nl_mz if any(abs(mz - f) < tol for f in nlmz))
            goodnl.append(match_count_nl/len(nlmz))
            # if match_count_nl == len(prodmz):
            #     goodnl.append(1)
            # else:
            #     goodnl.append(0)

        print(f"Row {index}: Class={class_name} Adduct={adduct}, IonMode={ionmode}, Good={good[-1]} GoodNL={goodnl[-1]}")

    df["HG_Fragment_Match"] = np.array(good, dtype=float).round(2)
    df["HG_NL_Match"] = np.array(goodnl, dtype=float).round(2)

    # Set match counts > 1 to 1
    df.loc[df["HG_Fragment_Match"] > 1, "HG_Fragment_Match"] = 1.0
    df.loc[df["HG_NL_Match"] > 1, "HG_NL_Match"] = 1.0
    return df


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
posdf["HG_NL_Category"] = posdf["HG_NL_Match"].apply(categorize_nl)
negdf["HG_NL_Category"] = negdf["HG_NL_Match"].apply(categorize_nl)

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



