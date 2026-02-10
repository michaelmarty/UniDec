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

def good_percentage_plot(df):
    plt.figure(figsize=(8, 6))
    dataset_counts = sorted(df["Dataset Count"].unique())
    datasets = sorted(df["Dataset"].unique())
    percentages = []
    labels = []
    for count in dataset_counts:
        subset = df[df["Dataset Count"] == count]
        # Bad ones are where "Comment" contaings "Low Quality"
        bad_count = len(subset[subset["Comment"].str.contains("Low Quality", case=False, na=False)])
        total_count = len(subset)
        good_percentage = 100 * (total_count - bad_count) / total_count if total_count > 0 else 0
        percentages.append(good_percentage)
        labels.append(f"{count} datasets")

    for dataset in datasets:
        # Isolate all where only this dataset is present
        subset = df[df["Dataset Count"] == 1 & (df["Dataset"] == dataset)]
        bad_count = len(subset[subset["Comment"].str.contains("Low Quality", case=False, na=False)])
        total_count = len(subset)
        good_percentage = 100 * (total_count - bad_count) / total_count if total_count > 0 else 0
        percentages.append(good_percentage)
        print(f"{dataset} only: {good_percentage:.2f}% good IDs")
        labels.append(f"{dataset} only")

    # Repeat for any in the dataset
    for dataset in datasets:
        subset = df[df["Dataset"].str.contains(dataset)]
        bad_count = len(subset[subset["Comment"].str.contains("Low Quality", case=False, na=False)])
        total_count = len(subset)
        good_percentage = 100 * (total_count - bad_count) / total_count if total_count > 0 else 0
        percentages.append(good_percentage)
        print(f"{dataset} present: {good_percentage:.2f}% good IDs")
        labels.append(f"{dataset} Total")

    # plt.bar(dataset_counts, percentages, color="skyblue")
    plt.bar(labels, percentages)
    plt.xlabel("Dataset Count")
    plt.ylabel("Percentage of Good IDs")
    plt.title("Percentage of Good IDs by Dataset Count")
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.show()

mpl.use("WxAgg")

topdir = r"Z:\Group Share\Shivam\DDA"
os.chdir(topdir)

# df = pd.read_csv("Combined_DDA_RT_Corrected_FullNames.csv")
df = pd.read_excel("Combined_DDA_OutlierAnalysis_Results.xlsx", sheet_name="All Data")
# Drop all where Outlier_Result is 1
df = df[df["Outlier_Prediction"] == 0]
good_percentage_plot(df)
exit()
print(df.columns)
# Create logSN column
df["logSN"] = df["S/N average"].apply(lambda x: np.log10(x) if x > 0 else 0)
plot_attributes_by_datasetcount(df, attrlist=["S/N average", "logSN", "Total score", "Fill %"])

