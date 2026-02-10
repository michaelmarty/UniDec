import os
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("WxAgg")

def correct_rts(df, rtcol, ref_dataset="Stellar", namecol="Metabolite name", adductcol="Adduct type"):
    df = df.copy()
    refdf = df[df["Dataset"] == ref_dataset]
    # Plot 1-1 RT line for comparison
    plt.figure(figsize=(8, 8))
    plt.subplot(1,2,1)
    plt.plot([0, 20], [0, 20], "k--", label="1:1 line")
    plt.subplot(1,2,2)
    plt.plot([0, 20], [0, 20], "k--", label="1:1 line")
    # Take only unique names in the reference dataset, and their RTs
    for dataset in df["Dataset"].unique():
        if dataset == ref_dataset:
            continue
        # Get subset of data for this dataset
        subset = df[df["Dataset"] == dataset]

        xmatch = []
        ymatch = []
        for y in subset.iterrows():
            match = refdf[(refdf[namecol] == y[1][namecol]) & (refdf[adductcol] == y[1][adductcol])]
            if len(match) > 0:
                xmatch.append(match.iloc[0][rtcol])
                ymatch.append(y[1][rtcol])

        print(len(xmatch), len(ymatch))

        # Sort the matches by xmatch
        xmatch = np.array(xmatch)
        ymatch = np.array(ymatch)
        sort_idx = np.argsort(xmatch)
        xmatch = xmatch[sort_idx]
        ymatch = ymatch[sort_idx]
        # Drop any duplicates in xmatch and ymatch
        _, unique_idx = np.unique(xmatch, return_index=True)
        xmatch = xmatch[unique_idx]
        ymatch = ymatch[unique_idx]

        # Remove outliers that are more than 1 min away from the 1-1 line
        median_shift = np.median(ymatch - xmatch)
        print(f"Median RT shift for {dataset}: {median_shift:.2f} minutes")
        # Remove points that are more than 1 min away from the median shift
        mask = np.abs((ymatch - xmatch) - median_shift) < 1
        xmatch = xmatch[mask]
        ymatch = ymatch[mask]

        plt.subplot(1,2,1)
        # Plot these RTs to check for linearity
        plt.scatter(xmatch, ymatch, label=dataset)
        plt.xlabel("Reference RT")
        plt.ylabel("Dataset RT")
        print("SSE:", np.sum((xmatch - ymatch) ** 2))

        # cs = CubicSpline(xmatch, ymatch, bc_type='natural')
        fit = np.polyfit(ymatch, xmatch, deg=3)
        correction = np.poly1d(fit)
        # Define a function to apply the correction
        newrts = correction(df.loc[df["Dataset"] == dataset, rtcol])
        deltart = newrts - df.loc[df["Dataset"] == dataset, rtcol]
        print(f"Mean RT shift for {dataset}: {np.mean(deltart):.2f} minutes")
        df.loc[df["Dataset"] == dataset, rtcol] = newrts
        df.loc[df["Dataset"] == dataset, "RT Shift"] = deltart

        # Plot corrected values
        plt.subplot(1,2,2)
        plt.scatter(xmatch, correction(ymatch))
        print("SSE after correction:", np.sum((xmatch - correction(ymatch)) ** 2))
        plt.xlabel("Reference RT")
        plt.ylabel("Corrected Dataset RT")
    plt.subplot(1,2,1)
    plt.legend()
    plt.show()
    return df

topdir = r"Z:\Group Share\Shivam\DDA"
os.chdir(topdir)

df = pd.read_csv("Combined_DDA.csv")

# Drop all but Stellar datasets
# df = df[df["Dataset"].isin(["ITIT", "Stellar"])]

correct_rts(df, "Average Rt(min)", ref_dataset="Stellar")

# Write out corrected dataframe for reference
df.to_csv("Combined_DDA_RT_Corrected.csv", index=False)