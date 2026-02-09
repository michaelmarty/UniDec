import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def plot_attributes_by_datasetcount(df, attrlist = []):
    # Plot attributes by dataset count
    plt.figure(figsize=(12, 6))
    for attr in attrlist:
        plt.subplot(1, len(attrlist), attrlist.index(attr) + 1)
        for count in sorted(df["Dataset Count"].unique()):
            subset = df[df["Dataset Count"] == count]
            # Box wisker plot of attribute for this dataset count
            plt.boxplot(subset[attr], positions=[count], widths=0.6)
        plt.xlabel(attr)
        plt.ylabel("Dataset Count")
        plt.legend()
        plt.title(f"{attr} by Dataset Count")
    plt.tight_layout()
    plt.show()

mpl.use("WxAgg")

topdir = r"Z:\Group Share\Shivam\DDA"
os.chdir(topdir)

df = pd.read_csv("Combined_DDA_RT_Corrected_FullNames.csv")
print(df.columns)
# Create logSN column
df["logSN"] = df["S/N average"].apply(lambda x: np.log10(x) if x > 0 else 0)
plot_attributes_by_datasetcount(df, attrlist=["S/N average", "logSN", "Total score", "Fill %"])

