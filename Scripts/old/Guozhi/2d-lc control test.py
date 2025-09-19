import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from matplotlib import rcParams


rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5

def get_pvalues(df, mc1, mc2, n1=4, n2=1, norm=True):
    m1 = df[mc1].to_numpy()
    s1 = df[mc1 + '.1'].to_numpy()
    m2 = df[mc2].to_numpy()
    s2 = df[mc2 + '.1'].to_numpy()
    if norm:
        norm1 = np.mean(m1)
        m1 = m1/norm1
        s1 = s1/norm1
        norm2 = np.mean(m2)
        m2 = m2/norm2
        s2 = s2/norm2

    fold = m2 / m1
    df["Fold"] = fold.tolist()
    df["Log2Fold"] = pd.Series(np.log2(df["Fold"]))
    pvals = [stats.ttest_ind_from_stats(m1[i], s1[i], n1, m2[i], s2[i], n2, equal_var=True)[1] for i in range(len(m1))]
    df["Pval"] = np.array(pvals).tolist()
    df["Log10P"] = pd.Series(-np.log10(df["Pval"]))
    df["FoldSD"] = fold * np.sqrt((s1/m1)**2 + (s2/m2)**2)
    return df

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

def norm_from_areas(df, files1, files2, norm=True):
    avals = []
    for f in files1:
        areas = df[f].to_numpy()
        if norm:
            norm = np.mean(areas)
            areas = areas/norm
        avals.append(areas)

    m1 = np.mean(avals, axis=0)
    s1 = np.std(avals, axis=0, ddof=1)

    avals = []
    for f in files2:
        areas = df[f].to_numpy()
        if norm:
            norm = np.mean(areas)
            areas = areas / norm
        avals.append(areas)
    m2 = np.mean(avals, axis=0)
    s2 = np.std(avals, axis=0, ddof=1)

    n1 = len(files1)
    n2 = len(files2)

    fold = m2 / m1
    df["Fold"] = fold.tolist()
    df["Log2Fold"] = pd.Series(np.log2(df["Fold"]))
    pvals = [stats.ttest_ind_from_stats(m1[i], s1[i], n1, m2[i], s2[i], n2, equal_var=True)[1] for i in range(len(m1))]
    df["Pval"] = np.array(pvals).tolist()
    df["Log10P"] = pd.Series(-np.log10(df["Pval"]))
    df["FoldSD"] = fold * np.sqrt((s1 / m1) ** 2 + (s2 / m2) ** 2)

    return df


file = "Z:\\mtmarty\\Group Share\\Guozhi Zhang\\2D_LC\\04260427and0428 cobined\\Area_0_202253163.txt"
df = pd.read_csv(file, skiprows=4, delimiter="\t")

df = df[df["MS/MS matched"]]
df["Ontology"].fillna("No MS2", inplace=True)

print(np.unique(df["Ontology"]))
onts = ['EtherLPE', 'EtherPE', 'EtherPG', 'LPE', 'LPG', 'OxPE', 'PE', 'PEtOH', 'PG']
colors = ["orange", "g", "blue", "r", "purple", "black", "lightblue", "lightgreen", "pink"]



df = df.loc[df["Ontology"].isin(onts)]
df = df.loc[~df["Metabolite name"].str.contains("Unsettled")]#exclude ontology
df = df.loc[~df["Metabolite name"].str.contains("Unknown")]
df = set_tails(df)

df = df[df["S/N average"]>10]
df = df[df["Dot product"]>50]


print("Number of Species:", len(df))
print(df.columns)

files1 = ["20220426 GZ 2D LC Test_flattrap_Ecoli_noLock_Quad440_newMPAnewMPB_newlipid_Capillary2_Desol250_scan10and5_excl35_50pTIC_increaseCE_swoff50p",
          "20220427 GZ 2D LC Test_flattrap_Ecoli_noLock_Quad440_newMPAoldMPB_newlipid_Capillary1_Desol250_scan10and5_excl35_increaseCE_swoff50p",
          "20220427 GZ 2D LC Test_flattrap_Ecoli_noLock_Quad440_newMPAoldMPB_newlipid_Capillary1_Desol400_scan10and10_excl35_50pTICincreaseCE_swoff50p_2",
          "20220427 GZ 2D LC Test_sharptrap_Ecoli_withLock_Quad440_newMPAoldMPB_newlipid_Capillary1_Desol400_scan10and10_exclall_50pTICincreaseCE_swoff50p_4_DDA9_2ions"]

files2 = ["20220428_EcoliND"]


df = norm_from_areas(df, files1, files2)


psorted = np.sort(df["Pval"])


print(psorted)



len_pval=len(psorted)
print("len_pval:", len_pval)

rank =1
p_adj=[]

for p in psorted:
    fdr_adj_p=p*len_pval/rank
    rank += 1
    p_adj.append(
        {
            "p_val": p, "fdr_adj_p": fdr_adj_p
        }
    )
print(p_adj)

print(fdr_adj_p)

pcutoff = 0.01
fcutoff = 1.1
sncutoff = 10

plt.figure(figsize=(7,5))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)


ax1.scatter(df["Log2Fold"], df["Log10P"], color="gray")
ontlabels = []
filter = True
for i, o in enumerate(onts):
    if o == "nan":
        continue

    df2 = df[df["Ontology"] == o]
    if len(df2) < 1 and filter:
        continue

    b1 = df2["Log10P"] > -np.log10(pcutoff)
    df2 = df2[b1]
    if len(df2) < 1 and filter:
        continue

    b1 = df2["S/N average"] > sncutoff
    df2 = df2[b1]
    if len(df2) < 1 and filter:
        continue

    folds = df2["Log2Fold"]

    b2 = np.abs(df2["Log2Fold"]) < np.log2(fcutoff)
    b2 = np.logical_not(b2)
    df2 = df2[b2]

    ax1.scatter(df2["Log2Fold"], df2["Log10P"], label=o, color=colors[i])
    ax2.bar(len(ontlabels), np.mean(folds), yerr=np.std(folds), color=colors[i])
    ontlabels.append(o)

ax2.set_xticks(np.arange(len(ontlabels)))
ax2.set_xticklabels(ontlabels, rotation="vertical")
#ax1.vlines([-0.5,0.5],0, 1, transform=ax1.get_xaxis_transform(), linestyle="--", color= "blue")
#ax1.hlines([1], 0,1, transform=ax1.get_yaxis_transform(), linestyle="--", color= "blue")
ax1.set_xlabel("Log2Fold")
ax1.set_ylabel("Log10P")
try:
    ax1.legend()
except:
    pass

tlabels = []
filter = True
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
ax3.hlines(1, 0, len(tlabels), linestyle="--", color="blue")
ax3.set_xticklabels(tlabels, rotation="vertical")
ax3.set_ylabel("Fold Change")


tls = np.unique(df["Tail Lengths"])

tlabels = []
filter = True
for i, t in enumerate(tls):
    df2 = df[df["Tail Lengths"] == t]
    if len(df2) <= 1 and filter:
        continue

    b1 = df2["Log10P"] > -np.log10(pcutoff)
    df2 = df2[b1]
    if len(df2) <= 1 and filter:
        continue

    b1 = df2["S/N average"] > sncutoff
    df2 = df2[b1]
    if len(df2) <= 1 and filter:
        continue

    folds = df2["Log2Fold"]

    if len(folds) < 1:
        continue

    b2 = np.abs(df2["Log2Fold"]) < np.log2(fcutoff)
    b2 = np.logical_not(b2)
    df2 = df2[b2]

    ax4.bar(len(tlabels), np.mean(folds), yerr=np.std(folds))
    tlabels.append(t)

tls = np.unique(df["Tail Unsaturation"])

filter = True
for i, t in enumerate(tls):
    df2 = df[df["Tail Unsaturation"] == t]
    if len(df2) <= 1 and filter:
        continue

    b1 = df2["Log10P"] > -np.log10(pcutoff)
    df2 = df2[b1]
    if len(df2) <= 1 and filter:
        continue

    b1 = df2["S/N average"] > sncutoff
    df2 = df2[b1]
    if len(df2) <= 1 and filter:
        continue

    folds = df2["Log2Fold"]

    if len(folds) < 1:
        continue

    b2 = np.abs(df2["Log2Fold"]) < np.log2(fcutoff)
    b2 = np.logical_not(b2)
    df2 = df2[b2]

    ax4.bar(len(tlabels), np.mean(folds), yerr=np.std(folds))
    tlabels.append(t)


ax4.set_xticks(np.arange(len(tlabels)))
ax4.set_xticklabels(tlabels, rotation="vertical")
plt.tight_layout()
plt.show()
