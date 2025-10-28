import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def subset_tukey(df2, out="OX", input="LIPID"):
    tukey = pairwise_tukeyhsd(df2[out], df2[input], alpha=0.05)
    a = tukey.summary()
    html = a.as_html()
    df3 = pd.read_html(html, header=0)[0]
    df4 = df3[df3["reject"]]
    if False:
        print(df3)
        if len(df4) > 0:
            print(df4)
        else:
            print("Nothing significant found...")

    return df4


def run_tukey(df, out="OX", lname="LIPID", rname="RATIO"):
    ratios = np.unique(df[rname])
    lipids = np.unique(df[lname])

    dfout = pd.DataFrame()

    for r in ratios:
        print("RATIO:", r)
        df2 = df[df[rname] == r]
        df4 = subset_tukey(df2, out=out, input=lname)
        df4[rname] = r
        dfout = pd.concat([dfout, df4])

    for l in lipids:
        print("LIPID TYPE:", l)
        df2 = df[df[lname] == l]
        df4 = subset_tukey(df2, out=out, input=rname)
        df4[lname] = l
        dfout = pd.concat([dfout, df4])

    return dfout


def check_tukey(df, l1, r1, l2, r2, lname="LIPID", rname="RATIO"):
    if l1 == "E.coli":
        l1 = "Ecoli"
    if l2 == "E.coli":
        l2 = "Ecoli"
    if l1 == l2:
        df2 = df[df[lname] == l1]
        df2 = df2[pd.to_numeric(df2["group1"]) == float(r1)]
        df2 = df2[pd.to_numeric(df2["group2"]) == float(r2)]
        df3 = []
    elif r1 == r2:
        df2 = df[df[rname] == r1]
        df2 = df2[df2["group1"] == l1]
        df2 = df2[df2["group2"] == l2]

        df3 = df[df[rname] == r1]
        df3 = df3[df3["group1"] == l2]
        df3 = df3[df3["group2"] == l1]
    else:
        return False

    if len(df2) > 0 or len(df3) > 0:
        return True

    else:
        return False


if __name__ == "__main__":
    path = "Z:\\Group Share\\Deseree\\FPOP Paper 2\\4R_Raw_June14.csv"

    df = pd.read_csv(path)

    dfout = run_tukey(df)

    print(dfout)

    outpath = "Z:\\Group Share\\Deseree\\FPOP Paper 2\\4R_Raw_June14_Tukey.csv"
    dfout.to_csv(outpath)
