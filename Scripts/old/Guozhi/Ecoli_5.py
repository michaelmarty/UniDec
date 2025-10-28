import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import os
from scipy import stats
import unidec.tools as ud
from unidec.LipiDec.AddFeaturestoDF import *

os.chdir("Z:\\Group Share\\Guozhi Zhang\\2D_LC\\230425_maybe_LX_paper")

# Set Paths
path = "230418_EColi_LXMS_Processed.xlsx"

# Load the data
df = pd.read_excel(path)

comparisons = [["Lipid", "ND_minus"], ["Lipid", "ND_plus"],
               ["ND_minus", "ND_plus"]]

labels = [comp[0] + '-' + comp[1] for comp in comparisons]

names = df["Name"]

# Strip off the class from the name
classes = df["Class Name"]
uclasses = np.unique(classes)

print(uclasses)

# Drop classes from the uclasses array
drop_classes = ["LPC", "LPE", "LPE-O", "LPE-O/P", "PS-O/P"]
# drop_classes = []
for drop_class in drop_classes:
    uclasses = uclasses[uclasses != drop_class]

# Move SM in the uclasses array to after HexCer
if "SM" in uclasses and "HexCer" in uclasses:
    insert_index = np.where(uclasses == "HexCer")[0][0] + 1
    uclasses = np.delete(uclasses, np.where(uclasses == "SM"))
    uclasses = np.insert(uclasses, insert_index, "SM")
# Move TG to the front
if "TG" in uclasses:
    uclasses = np.delete(uclasses, np.where(uclasses == "TG"))
    uclasses = np.insert(uclasses, 0, "TG")

# include_classes = ["PE-O"]
# uclasses = include_classes
colors = get_colors(uclasses)

fold_limit = 0.30
ldiff_limt = 0.38
lunsat_limt = 0.07
capsize = 2

pval = 0.05

width = 4
figsize = (2.5, 6)

# Set the font size to 12
plt.rcParams.update({'font.size': 11})

titles = ["", "", "", "", ""]

if True:

    # Volcano Plots of the Data
    for j, label in enumerate(labels):
        figure = plt.figure(figsize=figsize)

        fc = df["Fold Change " + label]
        p = df["P Value " + label]

        axes = figure.subplots(3, 1)

        axes[0].set_title(titles[j])
        # Get average fold change for each class
        for i, class_ in enumerate(uclasses):
            subdf = df[classes == class_]
            fcs = subdf["Fold Change " + label]
            n = len(fcs)
            tval = abs(stats.t.ppf(pval / 2, n))
            class_fc = np.mean(fcs)
            class_std = np.std(fcs) * tval / np.sqrt(n)

            print(class_, class_fc, class_std, n)
            # Make a bar chart
            axes[0].bar(i, class_fc, yerr=class_std, color=colors[i], ecolor="black", capsize=capsize)

            # Add an astrisk if the confidence interval does not include 0
            if class_fc - class_std > 0 or class_fc + class_std < 0:
                if n >= 3 and np.abs(class_fc) > fold_limit:
                    markerypos = np.amax([class_fc + class_std, fold_limit])
                    axes[0].text(i, markerypos, "*", horizontalalignment="center")

            # Remove x-ticks
            axes[0].set_xticks([])
            axes[0].set_ylim(-0.5, 0.8)
            # Remove top and right axes
            axes[0].spines["top"].set_visible(False)
            axes[0].spines["right"].set_visible(False)

            tld = subdf["Tail Length Diff"]
            tud = subdf["Tail Unsat Diff"]
            tl = subdf["Tail Lengths"]
            tu = subdf["Tail Unsaturation"]
            # Remove -1 values
            tld = tld[tl != -1]
            tud = tud[tu != -1]
            fcs = fcs[tl != -1]

            if len(tld) == 0:
                avg_ldiff = 0
                std_ldiff = 0
            else:
                # Calculate the average tail length of the class
                arr = tld * fcs
                avg_ldiff = np.average(arr)
                std_ldiff = np.std(arr) * tval / np.sqrt(len(tld))
            print(class_, avg_ldiff, std_ldiff)

            if len(tud) == 0:
                avg_udiff = 0
                std_udiff = 0
            else:
                # Calculate the average tail unsaturation of the class
                arr = tud * fcs
                avg_udiff = np.average(arr)
                std_udiff = np.std(arr) * tval / np.sqrt(len(tud))
            print(class_, avg_udiff, std_udiff)

            # Plot the tail length and unsaturation
            axes[1].bar(i, avg_ldiff, yerr=std_ldiff, capsize=capsize, color=colors[i])
            axes[1].set_xticks([])
            axes[1].set_ylim(-1, 1)
            # Remove top and right axes
            axes[1].spines["top"].set_visible(False)
            axes[1].spines["right"].set_visible(False)

            # Add an astrisk if the confidence interval does not include 0
            if avg_ldiff - std_ldiff > 0 or avg_ldiff + std_ldiff < 0:
                if n >= 3 and np.abs(avg_ldiff) > ldiff_limt:
                    markerypos = np.amax([avg_ldiff + std_ldiff, ldiff_limt])
                    axes[1].text(i, markerypos, "*", horizontalalignment="center")

            axes[2].bar(i, avg_udiff, yerr=std_udiff, capsize=capsize, color=colors[i])
            axes[2].set_ylim(-1, 1)

            # Add an astrisk if the confidence interval does not include 0
            if avg_udiff - std_udiff > 0 or avg_udiff + std_udiff < 0:
                if n >= 3 and np.abs(avg_udiff) > lunsat_limt:
                    markerypos = np.amax([avg_udiff + std_udiff, lunsat_limt])
                    axes[2].text(i, markerypos, "*", horizontalalignment="center")

        # Label x ticks as uclasses
        # Add Dashed line through 0
        axes[0].axhline(0, color="black", linestyle="-")
        axes[1].axhline(0, color="black", linestyle="-")
        axes[2].axhline(0, color="black", linestyle="-")

        # Add dashed lines at limits
        lstyle = "dotted"
        axes[0].axhline(fold_limit, color="black", linestyle=lstyle)
        axes[0].axhline(-fold_limit, color="black", linestyle=lstyle)
        axes[1].axhline(ldiff_limt, color="black", linestyle=lstyle)
        axes[1].axhline(-ldiff_limt, color="black", linestyle=lstyle)
        axes[2].axhline(lunsat_limt, color="black", linestyle=lstyle)
        axes[2].axhline(-lunsat_limt, color="black", linestyle=lstyle)

        # Set Y labels
        axes[0].set_ylabel("Log2 Fold Change")
        axes[1].set_ylabel("Tail Length Diff.")
        axes[2].set_ylabel("Tail Unsat. Diff.")

        axes[2].set_xticks(np.arange(len(uclasses)))
        axes[2].set_xticklabels(uclasses, rotation=45, horizontalalignment="right", rotation_mode="anchor")
        # Remove top and right axes
        axes[2].spines["top"].set_visible(False)
        axes[2].spines["right"].set_visible(False)
        plt.tight_layout()
        #plt.savefig("Plot" + label + ".pdf")
        #plt.savefig("Plot" + label + ".png")
    plt.show()

