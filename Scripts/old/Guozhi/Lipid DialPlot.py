import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

def get_pvalues(df, mc1, sc1, mc2, sc2, n1, n2):
    m1 = df[mc1]#mean
    s1 = df[sc1]#standard d
    m2 = df[mc2]
    s2 = df[sc2]
    df["Fold"] = pd.Series(m1 / m2)#here why?
    df["Log2Fold"] = pd.Series(np.log2(df["Fold"]))
    pvals = [stats.ttest_ind_from_stats(m1[i], s1[i], n1, m2[i], s2[i], n2, equal_var=False)[1] for i in range(len(m1))]
    df["Pval"] = pd.Series(pvals)
    df["Log10P"] = pd.Series(-np.log10(df["Pval"]))
    return df


file = "Z:\\Group Share\\Guozhi Zhang\\MSDial Data\\Normalized_BrainPos.txt"
file = "Z:\\Group Share\\Guozhi Zhang\\SynaptXS\\0813_0818\\Time Exchange\\Normalized_2_20218191220.txt"


df = pd.read_csv(file, skiprows=4, delimiter="\t")
df = get_pvalues(df, "1", "1.1", "2", "2.1", 3, 3)#class name

#df = df[df["MS/MS matched"]]
#df = df[df["Ontology"] != "null"] #include==


df = df[df["S/N average"]>20]
df = df[df["Average Mz"]>450]

includelist=["PG", "PE"]
b1 = np.zeros_like(df["Ontology"])# == includelist[0]
for o in includelist:
    b1 = np.logical_or(b1, df["Ontology"] == o)
df = df[b1]

#excludelist=["FA", "Unknown"]
#for e in excludelist:
    #df = df[df["Ontology"] != e]




ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
onts = np.unique(df["Ontology"].astype(np.str))
ax1.scatter(df["Log2Fold"], df["Log10P"], color="gray")
for i, o in enumerate(onts):
    if o != "null":
        df2 = df[df["Ontology"] == o]

        b1 = df2["Log10P"] > 2
        b2 = np.abs(df2["Log2Fold"]) > 0.5
        #b2 = np.logical_not(b2)
        b2 = np.logical_and(b1, b2)
        df2 = df2[b2]

        folds = df2["Log2Fold"]
        print(np.mean(folds), o)


        ax1.scatter(df2["Log2Fold"], df2["Log10P"], label=o)
        ax2.bar(i, np.mean(folds), yerr=np.std(folds))
ax2.set_xticks(np.arange(len(onts)))
ax2.set_xticklabels(onts, rotation="vertical")
ax1.legend()
ax1.set_xlabel("Log2Fold")
ax1.set_ylabel("Log10P")
plt.show()
