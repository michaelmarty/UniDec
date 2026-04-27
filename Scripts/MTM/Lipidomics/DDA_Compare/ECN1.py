import pandas as pd
from PublicScripts.Lipids.LipidFunctions import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.use("WxAgg")

def fit_rt(x, y):
    fit = np.polyfit(x, y, 5)
    fitdata = np.polyval(fit, x)
    return fit, fitdata

def calc_ecn(df, tlname="Tail Length", tuname="Tail Unsaturation", k=2):
    # If TL is 0 or TU is -1 set ECN to -1
    df["ECN"] = df.apply(lambda row: -1 if row[tlname] == 0 or row[tuname] == -1 else row[tlname] - k * row[tuname],
                         axis=1)
    return df

def plot_ecn(df, ecn_col="ECN", rt_col="Average Rt(min)", hue_col="Tail Unsaturation"):
    plt.figure(figsize=(10, 6))
    for dataset in df[hue_col].unique():
        subset = df[df[hue_col] == dataset]
        plt.scatter(subset[ecn_col], subset[rt_col], label=dataset)

        xdata = subset[ecn_col]
        ydata = subset[rt_col]
        # Sort data by x values for fitting
        sorted_indices = np.argsort(xdata)
        xdata_sorted = xdata.iloc[sorted_indices]
        ydata_sorted = ydata.iloc[sorted_indices]

        fit , fitdata = fit_rt(xdata_sorted, ydata_sorted)
        plt.plot(xdata_sorted, fitdata, label=f"{dataset} Fit", linestyle='--')

    plt.title("ECN vs RT")
    plt.xlabel("ECN")
    plt.ylabel("Average RT (min)")
    plt.legend(title=hue_col)
    plt.show()

def align_ecn(df, ecn_col="ECN", rt_col="Average Rt(min)", ref_ecn=0):
    return df


file = r"Z:\Group Share\Annika\DDA_Compare\mzmine\Combined_DDA_OutlierAnalysis.xlsx"

df = pd.read_excel(file)
df = apply_ontology_mapping(df)
# print(df.columns)

# Drop all "[M+Na]+" adducts
df = df[df["Adduct type"] != "[M+Na]+"]

# Filter to PCs
df = df[df["Ontology"] == "PG"]
# Filter dataset to "OTOT"
df = df[df["Dataset"] == "OTOT"]

df = set_tails(df, "Metabolite name")

df = calc_ecn(df)

# Filter for unsaturation of 0
# df = df[df["Tail Unsaturation"] == 0]

print(df.to_string())

# Drop rows with ECN -1
df = df[df["ECN"] > 0]

plot_ecn(df)
