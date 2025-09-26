import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from matplotlib import rcParams

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

def filter_df(df, onts):
    df = df[df["MS/MS matched"]]
    df["Ontology"].fillna("No MS2", inplace=True)

    df = df.loc[df["Ontology"].isin(onts)]
    df = set_tails(df)

    df = df[df["Tail Lengths"] > 24]
    df = df[df["S/N average"] > 10]
    df = df[df["Dot product"] > 50]
    return df

#rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5
#rcParams['font.size']=8


posfile = "Z:\\Group Share\\Guozhi Zhang\\SynaptXS\\BrianDecRerun\\2\\processed.txt"
negfile = "Z:\\Group Share\\Guozhi Zhang\\SynaptXS\\BrianDecRerun\\2\\processedneg.txt"

df = pd.read_csv(posfile)
ndf = pd.read_csv(negfile)

onts = ["PE", "PG", "CL", "LPE", "PS", "PI", "PC"]
colors = ["orange", "g", "black", "r", "purple", "blue", "brown"]

posdf = filter_df(df, onts)
negdf = filter_df(ndf, onts)

df = pd.concat([posdf, negdf])

plt.figure(figsize=(7,5))
ax3=plt.gca()
tlabels = []
filter = True
pcutoff = 0.05
fcutoff = 1.1
sncutoff = 10
sort = df.sort_values("Log2Fold")
for index, row in sort.iterrows():
    met = row["Metabolite name"]
    met = met.split("|")[0]
    logp = row["Log10P"]
    logfold = row["Log2Fold"]
    fold = row["Fold"]
    sd = row["FoldSD"]
    if (logp <= -np.log10(pcutoff) or row["S/N average"] <= sncutoff or np.abs(logfold) <= np.log2(fcutoff)) and filter:
        # Comment out to remove the grey bars
        ax3.bar(len(tlabels), fold, yerr=sd, color="grey")
        tlabels.append(met)
        pass
    else:
        ont = row["Ontology"]
        for x in range(0, len(onts)):
            if onts[x]==ont:
                break
        ax3.bar(len(tlabels), fold, yerr=sd, color=colors[x])
        tlabels.append(met)
ax3.set_xticks(np.arange(len(tlabels)))
ax3.hlines(1, 0, len(tlabels), linestyle="--")
ax3.set_xticklabels(tlabels, rotation="vertical")
ax3.set_ylabel("Fold Change")

plt.tight_layout()
plt.show()
