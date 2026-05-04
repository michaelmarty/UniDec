import os
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

    ecn_range = np.sort(df[ecn_col].unique())
    ecn_avg = []
    for r in ecn_range:
        ecn_subset = df[df[ecn_col] == r]
        ecn_avg.append(ecn_subset[rt_col].mean())

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
    plt.plot(ecn_range, ecn_avg, label="Average RT", linestyle='-')
    plt.title("ECN vs RT")
    plt.xlabel("ECN")
    plt.ylabel("Average RT (min)")
    plt.legend(title=hue_col)
    plt.show()

def plot_rtmatch(df, rt_col="Average Rt(min)", ref_rt_col="Reference RT", hue_col="Tail Unsaturation"):
    plt.figure(figsize=(10, 6))
    plt.subplot(121)
    for dataset in df[hue_col].unique():
        subset = df[df[hue_col] == dataset]
        plt.scatter(subset[ref_rt_col], subset[rt_col], label=dataset)

    plt.title("RT Match")
    plt.xlabel("Reference RT (min)")
    plt.ylabel("Average RT (min)")
    plt.legend(title=hue_col)
    plt.subplot(122)
    # Plot RTdiff vs Reference RT
    for dataset in df[hue_col].unique():
        subset = df[df[hue_col] == dataset]
        plt.scatter(subset[ref_rt_col], subset["RTdiff"], label=dataset)
    plt.show()

def calc_delta_rt_ecn(df, ecn_col="ECN", rt_col="Average Rt(min)"):
    # Calculate average RT for each ECN value
    ecn_rt_means = df.groupby(ecn_col)[rt_col].mean()

    # Map average RT back to each row based on ECN and calculate delta RT
    df["ECN_RT_Mean"] = df[ecn_col].map(ecn_rt_means)
    df["Delta_RT"] = df[rt_col] - df["ECN_RT_Mean"]
    # Drop the intermediate ECN_RT_Mean column
    df.drop(columns=["ECN_RT_Mean"], inplace=True)

    return df




file = r"Z:\Group Share\Annika\DDA_Compare\mzmine\Combined_DDA_OutlierAnalysis.xlsx"
os.chdir(os.path.dirname(file))

df = pd.read_excel(file)
# df = apply_ontology_mapping(df)
# print(df.columns)

# Print all unique ontologies
print("Unique ontologies:", np.sort(df["Ontology"].unique()))

# Drop all "[M+Na]+" adducts
df = df[df["Adduct type"] != "[M+Na]+"]

# Filter to PCs
# df = df[df["Ontology"] == "SM"]
# Filter dataset to "OTOT"
# df = df[df["Dataset"] == "OTOT"]

# print(df[["Metabolite name", "Ontology", "Adduct type"]].to_string())


df = set_tails(df, "Metabolite name")

df = calc_ecn(df)
df = calc_delta_rt_ecn(df)
# Filter for unsaturation of 0
# df = df[df["Tail Unsaturation"] == 0]

# print(df.to_string())

# Drop rows with ECN -1
# df = df[df["ECN"] > 0]

# Write output
df.to_excel("ECN_Analysis.xlsx", index=False)

plot_ecn(df)
plot_rtmatch(df)
