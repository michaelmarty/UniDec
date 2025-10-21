from scipy import stats
import scipy.stats
import matplotlib as mpl
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import pandas as pd
import numpy as np
from Scripts.LipiDec.AddFeaturestoDF import *


def make_fc_barchart(classes, df, groups, column_comparisons, alpha=0.05):
    '''
    Makes a bar chart of the difference in lipid class composition between each comparison group.
    :param names: Names of the lipids
    :param classes: Classes of the lipids to be compared
    :param df: DataFrame containing the lipid data ("Lipid Name" column must be present)
    :param column_comparisons: Two lists of strings (column names) to be compared
    '''

    colors = get_colors(classes)

    class_areas1 = [np.zeros(len(column_comparisons[0])) for i in range(len(classes))]
    class_areas2 = [np.zeros(len(column_comparisons[1])) for i in range(len(classes))]

    for idx, row in df.iterrows():
        name = row["Lipid Name"]
        lipid_class = get_class(name)
        for i in range(len(classes)):
            if lipid_class == classes[i]:
                for j in range(len(column_comparisons[0])):
                    if np.isnan(row[column_comparisons[0][j]]):
                        continue
                    else:
                        class_areas1[i][j] += row[column_comparisons[0][j]]
                for j in range(len(column_comparisons[1])):
                    if np.isnan(row[column_comparisons[1][j]]):
                        continue
                    else:
                        class_areas2[i][j] += row[column_comparisons[1][j]]

    # normalize the areas
    for i in range(len(column_comparisons[0])):
        total = 0
        for j in range(len(classes)):
            total += class_areas1[j][i]
        for j in range(len(classes)):
            class_areas1[j][i] = class_areas1[j][i] / total

    for i in range(len(column_comparisons[1])):
        total = 0
        for j in range(len(classes)):
            total += class_areas2[j][i]
        for j in range(len(classes)):
            class_areas2[j][i] = class_areas2[j][i] / total

    fcs = np.zeros(len(classes))
    stds = np.zeros(len(classes))
    margins = np.zeros(len(classes))
    ddof = 2

    for i in range(len(classes)):
        mean1 = np.mean(class_areas1[i])
        mean2 = np.mean(class_areas2[i])
        std1 = np.std(class_areas1[i], ddof=ddof)
        std2 = np.std(class_areas2[i], ddof=ddof)

        fcs[i] = mean1 / mean2
        stds[i] = fcs[i] * np.sqrt((std1 / mean1) ** 2 + (std2 / mean2) ** 2)
        margins[i] = (stds[i] / np.sqrt(len(column_comparisons[0]) - ddof)) * z_critical

    #z

    # correct the p-values

    rejected, corrected_pvalues, _, _ = multipletests(p_vals, method='fdr_bh')
    for i, lipid_class in enumerate(classes):
        print(
            f"{lipid_class}: raw p={p_vals[i]:.4f}, corrected p={corrected_pvalues[i]:.4f}, significant={rejected[i]}")
    for i in range(len(classes)):
        if corrected_pvalues[i] < 0.05:
            print(
                f"Significant difference in {classes[i]} content between {groups[0]} and {groups[1]} (p = {corrected_pvalues[i]:.4f})")
        else:
            print(
                f"No significant difference in {classes[i]} content between {groups[0]} and {groups[1]} (p = {corrected_pvalues[i]:.4f})")

    sig_labels = []
    for p in corrected_pvalues:
        if p < 0.001:
            sig_labels.append('***')
        elif p < 0.01:
            sig_labels.append('**')
        elif p < 0.05:
            sig_labels.append('*')
        else:
            sig_labels.append('ns')

    # ---- Print results ----
    print("\nZ-test results with Benjamini–Hochberg correction:")
    print("--------------------------------------------------")
    for i in range(len(classes)):
        print(
            f"{classes[i]}: Z = {zscores[i]:.2f}, raw p = {p_vals[i]:.4g}, BH p = {corrected_pvals[i]:.4g}, FC = {fcs[i]:.3f}")

    plt.figure(figsize=(5, 5))

    plt.title(str(groups[0]) + " vs " + str(groups[1]), fontsize=26)
    plt.xlabel("Lipid Class", fontsize=22)
    plt.ylabel("Fold Change\nprotein/membrane", fontsize=22)

    # horizontal line at y=0 (solid)
    plt.hlines(0, -1, len(classes), color='black', linestyle='-')


    n_classes = len(classes)

    x = np.arange(len(classes)) * 0.5
    bar_width = 0.3
    plt.bar(x, fcs, width=bar_width, color=colors, edgecolor='black', linewidth=2, capsize=7, zorder=2)

    # Dotted baseline at y=1, with larger and more spread-out dots
    n_dots = 40  # fewer dots = more spacing
    x_dots = np.linspace(-0.9, n_classes - 0.2, n_dots)  # spreads across the bars
    y_dots = np.ones_like(x_dots) * 1  # y=1 baseline
    plt.plot(x_dots, y_dots, 'o', color='black', markersize=4, zorder=1)

    # Draw bars with thicker error bars and caps
    plt.bar(
        x,
        fcs,
        yerr=margins,
        color=colors,
        capsize=8,  # longer caps
        edgecolor='black',
        linewidth=2,
        width=bar_width,
        zorder=2,
        error_kw={'elinewidth': 2.5, 'capthick': 2.5}  # thicker error lines and caps
    )

    # Replicate points (unfilled circles) above bars
    for i in range(len(classes)):
        g1 = np.array(class_areas1[i])
        fc_points = g1 / np.mean(class_areas2[i])
        jitter = np.random.normal(0, 0.025, len(fc_points))  # small jitter
        plt.scatter(
            x[i] + jitter,
            fc_points,
            facecolors='none',
            edgecolors='black',
            s=25,
            zorder=3  # above bars
        )

    # Labels and ticks
    plt.xticks(x, classes, fontsize=10)
    plt.yticks([0, 0.5, 1.0, 1.5], fontsize=12)
    plt.ylabel("Fold Change\nprotein/membrane", fontsize=14)
    plt.xlabel("Lipid Class", fontsize=14)
    plt.xlim(-0.3, len(classes) - 1 + 0.3)
    plt.title(f"{groups[0]} vs {groups[1]}", fontsize=16)


    # calculate p-values and assign symbols
    #for i in range(len(classes)):
        #g1 = class_areas1[i]
       # g2 = class_areas2[i]
       # t, p = ttest_ind(g1, g2, equal_var=False)
       # if p < 0.001:
       #     sig_labels.append('***')
       # elif p < 0.01:
      #      sig_labels.append('**')
       # elif p < 0.05:
      #      sig_labels.append('*')
     #   else:
      #      sig_labels.append('ns')

    # add the stars/ns above bars
    #for i in range(len(sig_labels)):
        # replicate points
    #    g1 = np.array(class_areas1[i])
      #  g2 = np.array(class_areas2[i])

        # compute the fold change points plotted
     #   fc_points1 = g1 / np.mean(class_areas2[i])
     #   fc_points2 = g2 / np.mean(class_areas2[i])

      #  all_points = np.concatenate([fc_points1, fc_points2])

        # top of error bar
     #  top_error = fcs[i] + margins[i]

        # highest point between replicates and error bar
      #  max_y = max(np.max(all_points), top_error)

        # add a small offset
   #     y = max_y + 0.05  # tweak as needed
    #    plt.text(
       #     x[i], y, sig_labels[i],
           # ha='center', va='bottom',
          #  fontsize=20,  # bigger stars
           # fontname='Arial'
        #    )
    plt.tight_layout()
    plt.show()

z_critical = scipy.stats.norm.ppf(1 - .01 / 2)

# Set backend to Agg
mpl.use('WxAgg')

#file_path = "Z:\\Group Share\\Synapt\\Collaborator Samples (CS)\\lipidomic_analysis_results.xlsx"

df = pd.read_excel("Z:\\Group Share\\JGP\\Christina_letAB_MlaC_pqiAB_data\\Analysis_and_Other\\normalized_abundances.xlsx")

#need to compare samples 2 and 3 in each group.

column_comparisons2 = [("S5_R1", "S5_R2", "S5_R3", "S5_R4", "S5_R5", "S5_R6", "S5_R7", "S5_R8", "S5_R9"),
                       ("S6_R1", "S6_R2", "S6_R3", "S6_R4", "S6_R5", "S6_R6", "S6_R7", "S6_R8", "S6_R9")]

#t_critical = stats.t.ppf(1-0.025, 16)


make_fc_barchart(["CL", "PE", "PG"], df, ["LetAB Protein","LetAB Membrane"],
                 column_comparisons2, z_critical)

