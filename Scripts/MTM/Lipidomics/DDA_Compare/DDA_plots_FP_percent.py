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

def good_percentage_plot(df, datasets=None, type="All", figname="Good_Percentage_by_Dataset", use_simple_names=False):
    if use_simple_names:
        simpnames = df["Metabolite name"].apply(lambda x: x.split("|")[0] if "|" in x else x)
        df["Simple name"] = simpnames
        namecol = "Simple name"
    else:
        namecol = "Metabolite name"

    if datasets is None:
        dataset_counts = sorted(df["Dataset Count"].unique())
        datasets = sorted(df["Dataset"].unique())
    else:
        dataset_counts = []
    percentages = []
    labels = []
    totalcounts = []
    for count in dataset_counts:
        subset = df[df["Dataset Count"] == count]
        # Drop duplicate names based on the name column
        subset = subset.drop_duplicates(subset=[namecol])
        # Bad ones are where "Comment" contaings "Low Quality"
        bad_count = len(subset[subset["Comment"].str.contains("Low Quality", case=False, na=False)])
        total_count = len(subset)
        good_percentage = 100 * (total_count - bad_count) / total_count if total_count > 0 else 0
        percentages.append(good_percentage)
        totalcounts.append(total_count)
        labels.append(f"{count} datasets")
        print(f"{count} datasets: {good_percentage:.2f}% good IDs. Total Number: {total_count}")

    if type == "All" or type == "Only":
        for dataset in datasets:
            # Isolate all where only this dataset is present
            subset = df[df["Dataset Count"] == 1 & (df["Dataset"] == dataset)]
            # Drop duplicate names based on the name column
            subset = subset.drop_duplicates(subset=[namecol])

            bad_count = len(subset[subset["Comment"].str.contains("Low Quality", case=False, na=False)])
            total_count = len(subset)
            good_percentage = 100 * (total_count - bad_count) / total_count if total_count > 0 else 0
            percentages.append(good_percentage)
            totalcounts.append(total_count)
            print(f"{dataset} only: {good_percentage:.2f}% good IDs. Total Number: {total_count}")
            labels.append(f"{dataset} only")

    if type == "All" or type == "Total":
        # Repeat for any in the dataset
        for dataset in datasets:
            subset = df[df["Dataset"].str.contains(dataset)]
            # Drop duplicate names based on the name column
            subset = subset.drop_duplicates(subset=[namecol])

            bad_count = len(subset[subset["Comment"].str.contains("Low Quality", case=False, na=False)])
            total_count = len(subset)
            good_percentage = 100 * (total_count - bad_count) / total_count if total_count > 0 else 0
            percentages.append(good_percentage)
            totalcounts.append(total_count)
            print(f"{dataset} present: {good_percentage:.2f}% good IDs. Total Number: {total_count}")
            labels.append(f"{dataset}")

    # Add slash in labels for OTOT, OTIT, and ITIT
    labels = [label.replace("OTOT", "OT/OT").replace("OTIT", "OT/IT").replace("ITIT", "IT/IT") for label in labels]

    # plt.bar(dataset_counts, percentages, color="skyblue")
    plt.figure(figsize=(10, 6))
    plt.subplot(122)
    plt.bar(labels, percentages)
    plt.xlabel("Dataset")
    plt.ylabel("Percentage of True Positives")
    # plt.title("Percentage of Good IDs by Dataset Count")
    # plt.xticks(rotation=45, ha="right")
    pmin = np.floor(np.amin(percentages) * 10) / 10 - 5
    pmin = 50
    plt.ylim(pmin, 100)

    # plt.subplot(122)
    # plt.bar(labels, totalcounts, color="lightgreen")
    # plt.xlabel("Dataset Count")
    # plt.ylabel("Total Count of IDs")
    # plt.title("Total Count of IDs by Dataset Count")
    # plt.xticks(rotation=45, ha="right")

    plt.subplot(121)
    # Number of good and bad IDs in different colors on stacked bar chart
    plt.bar(labels, totalcounts, color="forestgreen", label="True")
    bad_counts = [total - (total * good / 100) for total, good in zip(totalcounts, percentages)]
    plt.bar(labels, bad_counts, color="red", label="False")
    plt.xlabel("Dataset")
    plt.ylabel("Count of IDs")
    # plt.title("Good vs Bad IDs by Dataset Count")
    # plt.xticks(rotation=45, ha="right")
    plt.legend(loc="upper left")
    plt.ylim(0, 1400)
    plt.tight_layout()
    # plt.savefig(figname + ".png", dpi=600, bbox_inches="tight")
    # plt.savefig( figname + ".pdf")


mpl.use("WxAgg")

topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Annika\IQ-X\Feb 2026 UntDDA\MSDial\Marked Low Score\RTP2"
topdir = r"Z:\Group Share\Shivam\IQX Stellar 3\FinalResults\AS"
os.chdir(topdir)

# df = pd.read_csv("Combined_DDA_RT_Corrected_FullNames.csv")
df = pd.read_excel("Combined_DDA_OutlierAnalysis_Results.xlsx", sheet_name="All Data")
df = df[df["Ontology"] != "BRSE"]
df = df[df["Ontology"] != "HBMP"]
# Drop all where Outlier_Result is 1
# df = df[df["Outlier_Prediction"] == 0]
# good_percentage_plot(df)
good_percentage_plot(df, datasets=["OTOT", "OTIT", "ITIT", "Stellar"], type="Total")

df = df[df["Outlier_Prediction"] == 0]
good_percentage_plot(df, datasets=["OTOT", "OTIT", "ITIT", "Stellar"], type="Total", figname="Good_Percentage_by_Dataset_OutlierFiltered")
plt.show()
exit()
print(df.columns)
# Create logSN column
df["logSN"] = df["S/N average"].apply(lambda x: np.log10(x) if x > 0 else 0)
plot_attributes_by_datasetcount(df, attrlist=["S/N average", "logSN", "Total score", "Fill %"])

