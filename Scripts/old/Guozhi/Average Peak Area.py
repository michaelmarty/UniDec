import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from matplotlib import rcParams


rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5




def get_overall_tails(name):
    try:
        front = name.split("|")[0][-4:]
        tl = front.split(":")[0]
        tu = front.split(":")[1]
    except:
        return name, -1, -1
    try:
        tl = int(tl)
    except:
        tl= -1
    try:
        tu = int(tu)
    except:
        tu = -1
    return front, tl, tu

def set_tails(df):
    names = df["Metabolite name"]
    tus = []
    tls = []
    tnames = []
    for n in names:
        total, tl, tu = get_overall_tails(n)
        tus.append(tu)
        tls.append(tl)
        tnames.append(total)
    df["Tail Lengths"] = tls
    df["Tail Unsaturation"] = tus
    df["Tail Names"] = tnames
    return df

def norm_from_areas(df, files1, norm=True):
    avals = []
    for f in files1:
        areas = df[f].to_numpy()
        if norm:
            norm = np.mean(areas)
            areas = areas/norm
        avals.append(areas)

    m1 = np.mean(avals, axis=0)
    s1 = np.std(avals, axis=0, ddof=1)

    n1 = len(files1)

    return df


file = "Z:\\mtmarty\\Group Share\\Guozhi Zhang\\2D_LC\\20220713 DIA and DDA comparison\\Area_1_20227251518.txt"

df = pd.read_csv(file, skiprows=4, delimiter="\t")

df = df[df["MS/MS matched"]]
df["Ontology"].fillna("No MS2", inplace=True)

print(np.unique(df["Ontology"]))
onts = ["PE", "PG", "CL", "EtherPE", "LPE", "LPG", "OxPE", "OxPG"]
colors = ["orange", "g", "blue", "r", "purple", "blk", "lightblue", "lightgreen"]



df = df.loc[df["Ontology"].isin(onts)]
df = df.loc[~df["Metabolite name"].str.contains("Unsettled")]
df = df.loc[~df["Metabolite name"].str.contains("Unknown")]
df = set_tails(df)

#df = df[df["Tail Lengths"]>5]
df = df[df["S/N average"]>10]
df = df[df["Dot product"]>20]


print("Number of Species:", len(df))


files1 = ["20220713_Neg_DDA_10Lipid20ISD_newbuffer_Quad440_run1",
          "20220713_Neg_DDA_10Lipid20ISD_newbuffer_Quad440_run2",
          "20220713_Neg_DDA_10Lipid20ISD_newbuffer_Quad440_run3"]



df = norm_from_areas(df, files1)

print(df)

df.to_excel("Z:\\mtmarty\\Group Share\\Guozhi Zhang\\2D_LC\\20220713 DIA and DDA comparison\\output_true.xlsx")


