from scipy import stats
import matplotlib as mpl
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
from Scripts.LipiDec.AddFeaturestoDF import *


def make_fc_barchart(classes, df, groups, column_comparisons):
    '''
    Makes a bar chart of the difference in lipid class composition between each comparison group.
    :param names: Names of the lipids
    :param classes: Classes of the lipids to be compared
    :param df: DataFrame containing the lipid data ("Lipid Name" column must be present)
    :param column_comparisons: Two lists of strings (column names) to be compared
    '''

    colors = [ '#fd9121', '#922290', '#0f8f53' ]


    class_areas1 = [np.zeros(len(column_comparisons[0])) for i in range(len(classes))]
    class_areas2 = [np.zeros(len(column_comparisons[1])) for i in range(len(classes))]

    for idx,row in df.iterrows():
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

    #normalize the areas
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
    zscore = np.zeros(len(classes))

    confidence_level = 0.95
    z_critical = stats.norm.ppf(1 - (1 - confidence_level) / 2)


    p_vals = np.zeros(len(classes))
    for i in range(len(classes)):
        n = len(column_comparisons[0])
        mean1 = np.mean(class_areas1[i])
        mean2 = np.mean(class_areas2[i])
        std1 = np.std(class_areas1[i], ddof=1) / np.sqrt(n)
        std2 = np.std(class_areas2[i], ddof=1) / np.sqrt(n)

        fcs[i] = mean1 / mean2
        stds[i] = fcs[i] * np.sqrt((std1/mean1)**2 + (std2/mean2)**2)

        margins[i] = stds[i] * z_critical
        zscore[i] = (fcs[i] - 1.00) / stds[i]
        p_vals[i] = 2.0*(1.0-stats.norm.cdf(abs(zscore[i])))
        if p_vals[i] == 0:
            zs = zscore[i]
            print(zs)
            print(2* stats.norm.cdf(zs))
    print(p_vals)
    #correct the p-values
    rejected,corrected_pvalues,_,_ = multipletests(p_vals, method='fdr_bh')
    for i in range(len(classes)):
        if corrected_pvalues[i] < 0.05:
            print("Significant difference in " + classes[i] + " content between " + groups[0] + " and " + groups[1] + " (p = " + str(corrected_pvalues[i]) + ")")
        else:
            print("No significant difference in " + classes[i] + " content between " + groups[0] + " and " + groups[1] + " (p = " + str(corrected_pvalues[i]) + ")")

    def significance_stars(p_value):
        """Determine number of stars for given p-value."""
        if p_value < 0.0001:
            return '***'
        elif p_value < 0.001:
            return '**'
        elif p_value < 0.01:
            return '*'
        else:
            return 'ns'  # Not significant
    plt.figure(figsize=(6, 7))
    # horizontal line at y=0 (solid)
    plt.hlines(0, -1, len(classes), color='black', linestyle='-')

    n_classes = len(classes)

    x = np.arange(len(classes)) * 0.4
    bar_width = 0.25
    plt.bar(x, fcs, width=bar_width, color=colors, edgecolor='black', linewidth=2, capsize=7, zorder=2)

    # Dotted baseline at y=1, with larger and more spread-out dots
    n_dots = 70  # fewer dots = more spacing
    x_dots = np.linspace(-0.9, n_classes - 0.2, n_dots)  # spreads across the bars
    y_dots = np.ones_like(x_dots) * 1  # y=1 baseline
    plt.plot(x_dots, y_dots, 'o', color='black', markersize=6, zorder=1)

    # Draw bars with thicker error bars and caps
    plt.bar(
        x,
        fcs,
        yerr=margins,
        color=colors,
        capsize=10,  # longer caps
        edgecolor='black',
        linewidth=2,
        width=bar_width,
        zorder=2,
        error_kw={'elinewidth': 2, 'capthick': 2}  # thicker error lines and caps
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
            s=35,
            zorder=3  # above bars
        )
    # Add significance stars above the bars
    for i, p_value in enumerate(corrected_pvalues):  # Use corrected p-values
        stars = significance_stars(p_value)  # Determine significance stars
        if stars:  # Only add stars for significant results
            replicate_points = np.array(class_areas1[i]) / np.mean(class_areas2[i])
            max_height = max(fcs[i] + margins[i], max(replicate_points))  # Bar + error + replicates
            star_y_position = max_height + 0.2  # Add spacing above the maximum point
            plt.text(
                x[i],  # X position
                star_y_position,  # Y position dynamically calculated
                stars,  # Star symbols ('*', '**', '***')
                ha='center', va='bottom',  # Align stars above the bar
                fontsize=18, color='black'  # Customize font size and color
            )
    # Labels and ticks
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    # Make bottom and left spines thicker
    plt.gca().spines['bottom'].set_linewidth(2)  # Set bottom spine thickness
    plt.gca().spines['left'].set_linewidth(2)  # Set left spine thickness
    plt.xticks(x, classes, fontsize=16)
    plt.yticks([0, 0.5, 1.0, 1.5], fontsize=16)
    plt.ylabel("Fold Change\nprotein/membrane", fontsize=20)
    plt.xlabel("Lipid Class", fontsize=20)
    plt.xlim(x[0] - bar_width, x[-1] + bar_width)  # Trim excess space

    plt.ylim(0, 1.6)
    # Make tick marks thicker
    plt.tick_params(axis='both', width=2, length=10)  # Tick mark thickness
    plt.tight_layout()
    plt.show()

#Fig N

from scipy import stats
import matplotlib as mpl
from statsmodels.stats.multitest import multipletests
from matplotlib.patches import Patch
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
from Scripts.LipiDec.AddFeaturestoDF import *

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 18


def make_content_paired_barcharts(classes, df, column_comparisons):
    colors = [ '#fd9121', '#922290', '#0f8f53' ]

    group1_contents = [np.zeros(len(column_comparisons[0])) for i in range(len(classes))]
    group2_contents = [np.zeros(len(column_comparisons[1])) for i in range(len(classes))]

    for idx, row in df.iterrows():
        name = row["Lipid Name"]
        print(name)
        lipid_class = get_class(name)
        for i in range(len(classes)):
            if lipid_class == classes[i]:
                for j in range(len(column_comparisons[0])):
                    if np.isnan(row[column_comparisons[0][j]]):
                        continue
                    else:
                        group1_contents[i][j] += row[column_comparisons[0][j]]
                for j in range(len(column_comparisons[1])):
                    if np.isnan(row[column_comparisons[1][j]]):
                        continue
                    else:
                        group2_contents[i][j] += row[column_comparisons[1][j]]


    group1_means = np.zeros(len(classes))
    group2_means = np.zeros(len(classes))
    group1_stds = np.zeros(len(classes))
    group2_stds = np.zeros(len(classes))
    group1_margins = np.zeros(len(classes))
    group2_margins = np.zeros(len(classes))
    df = len(column_comparisons[0]) - 1
    alpha = 0.05  # for 95% confidence

    # Now perform a t-test to determine if the means are significantly different
    t_stats = np.zeros(len(classes))
    p_vals = np.zeros(len(classes))
    for i in range(len(classes)):
        t_result = stats.ttest_ind(group1_contents[i], group2_contents[i])
        p_vals[i] = t_result.pvalue
        dof = t_result.df
        t_critical = stats.t.ppf(1 - 0.025, dof)
        t_stats[i] = t_result.statistic
        group1_means[i] = np.mean(group1_contents[i]*100)
        group2_means[i] = np.mean(group2_contents[i]*100)
        group1_stds[i] = np.std(group1_contents[i]*100, ddof=1)
        group2_stds[i] = np.std(group2_contents[i]*100, ddof=1)
        group1_margins[i] = (group1_stds[i] / np.sqrt(len(column_comparisons[0])))*t_critical
        group2_margins[i] = (group2_stds[i] / np.sqrt(len(column_comparisons[1])))*t_critical


        #print(f"{classes[i]}: p = {p_vals[i]:.3e}")


        # Apply Benjamini-Hochberg correction
        rejected, corrected_pvals, _, _ = multipletests(p_vals, method='fdr_bh')

        # Print corrected p-values once

    for i in range(len(classes)):
        if rejected[i]:
            print(
                f"Significant difference in {classes[i]} content after BH correction (corrected p = {corrected_pvals[i]:.3e})")
        else:
            print(
                f"No significant difference in {classes[i]} content after BH correction (corrected p = {corrected_pvals[i]:.3e})")


    plt.figure(figsize=(6, 7))
    bar_width = 0.15
    spacing = 0.02
    class_spacing = 0.5
    r1 = np.arange(len(classes)) * class_spacing
    r2 = r1 + bar_width + spacing
    # Center of each class for labels
    # Center x-axis labels between the bars
    r_labels = (r1 + r2) / 2
    print(r_labels)
    plt.xticks(r_labels, classes, fontsize=18)

    plt.bar(r1, group1_means, yerr=group1_margins, color=colors, capsize=5, edgecolor='black',
                width=bar_width, label='LetAB Protein',linewidth=2, error_kw={'elinewidth':2, 'capthick':2})
    plt.bar(r2, group2_means, yerr=group2_margins, color='none', capsize=5, edgecolor=colors,
                hatch='\\', width=bar_width, label='LetAB Membrane',linewidth=2, error_kw={'elinewidth':2, 'capthick':2})
    mpl.rcParams['hatch.linewidth'] = 6  # or whatever thickness you want
    # Plot replicates for all lipids
    for i in range(len(classes)):
        # Protein replicates
        fc_points1 = np.array(group1_contents[i] * 100)
        if len(fc_points1) > 0:
            jitter1 = np.random.normal(0, 0.02, len(fc_points1))
            plt.scatter(
                r1[i] + jitter1,
                fc_points1,
                facecolors='none',
                edgecolors='black',
                s=25,
                zorder=4
            )

        # Membrane replicates
        fc_points2 = np.array(group2_contents[i] * 100)
        if len(fc_points2) > 0:
            jitter2 = np.random.normal(0, 0.02, len(fc_points2))
            plt.scatter(
                r2[i] + jitter2,
                fc_points2,
                facecolors='none',
                edgecolors='black',
                s=35,
                zorder=5
            )

    # Define legend once, outside the loop
    legend_handles = [
        Patch(facecolor='black', edgecolor='black', label='LetAB Protein'),  # solid bars
        Patch(facecolor='white', edgecolor='black', hatch='\\', label='LetAB Membrane')  # striped bars
    ]


    # Plot significance lines and stars
    plt.plot([r1[0], r1[0], r2[0], r2[0]],
             [group1_means[0] + 6, group2_means[0] + 10, group2_means[0] + 10, group2_means[0] + 6],
             color='black', linestyle='-', linewidth=2)
    plt.text((r1[0] + r2[0]) / 2, group2_means[0] + 12, "*", fontsize=18, ha='center')

    plt.plot([r1[1], r1[1], r2[1], r2[1]],
             [group1_means[1] + 6, group1_means[1] + 10, group1_means[1] + 10, group2_means[1] + 10],
             color='black', linestyle='-', linewidth=2)
    plt.text((r1[1] + r2[1]) / 2, group1_means[1] + 12, "*", fontsize=18, ha='center')


    legend_handles = [
        Patch(facecolor='white', edgecolor='black', label='LetAB Protein'),
        Patch(facecolor='white', edgecolor='black', label='LetAB Membrane', hatch='\\')
    ]

    plt.plot([r1[0],r1[0],r2[0],r2[0]],
             [group1_means[0]+6, group2_means[0]+10, group2_means[0]+10, group2_means[0]+6],
             color='black', linestyle='-', linewidth=2)
    plt.text((r1[0]+r2[0])/2, group2_means[0]+12, "*", fontsize=18, ha='center')

    plt.plot([r1[1], r1[1], r2[1], r2[1]],
             [group1_means[1] + 6, group1_means[1] + 10, group1_means[1] + 10, group2_means[1] + 10],
             color='black', linestyle='-', linewidth=2)
    plt.text((r1[1] + r2[1]) / 2, group1_means[1] + 12, "*", fontsize=18, ha='center')

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(2)  # Set bottom spine thickness
    plt.gca().spines['left'].set_linewidth(2)  # Set left spine thickness
    # plt.xticks([r + bar_width/2 for r in range(len(classes))], classes, fontsize=18)
    plt.xlabel('Lipid Class', fontsize=22)
    plt.ylim(0, 100)
    plt.yticks(np.arange(0, 101, 20), fontsize=18)
    plt.ylabel('Summed Intensity %', fontsize=22)
    plt.tight_layout()
    plt.legend(handles=legend_handles, fontsize=14,  loc='upper right', bbox_to_anchor=(1.05, 1.05))
    plt.tick_params(axis='both', width=2, length=10)  # Tick mark thickness
    plt.show()





# Set backend to Agg
mpl.use('WxAgg')

#file_path = "Z:\\Group Share\\Synapt\\Collaborator Samples (CS)\\lipidomic_analysis_results.xlsx"

df = pd.read_excel("Z:\\Group Share\\JGP\\Christina_letAB_MlaC_pqiAB_data\\Analysis_and_Other\\normalized_abundances.xlsx")

#need to compare samples 2 and 3 in each group.

column_comparisons2 = [("S5_R1", "S5_R2", "S5_R3", "S5_R4", "S5_R5", "S5_R6", "S5_R7", "S5_R8", "S5_R9"),
                       ("S6_R1", "S6_R2", "S6_R3", "S6_R4", "S6_R5", "S6_R6", "S6_R7", "S6_R8", "S6_R9")]

#t_critical = stats.t.ppf(1-0.025, 16)


make_fc_barchart(["CL", "PE", "PG"], df, ["LetAB Protein","LetAB Membrane"],
                 column_comparisons2)

df = pd.read_excel("Z:\\Group Share\\JGP\\Christina_letAB_MlaC_pqiAB_data\\Analysis_and_Other\\normalized_abundances.xlsx")

#need to compare samples 2 and 3 in each group.

column_comparisons1 = [("S2_R1", "S2_R2", "S2_R3", "S2_R4", "S2_R5", "S2_R6", "S2_R7", "S2_R8"),
                       ("S3_R1", "S3_R2", "S3_R3", "S3_R4", "S3_R5", "S3_R6", "S3_R7", "S3_R8", "S3_R9")]

column_comparisons2 = [("S5_R1", "S5_R2", "S5_R3", "S5_R4", "S5_R5", "S5_R6", "S5_R7", "S5_R8", "S5_R9"),
                       ("S6_R1", "S6_R2", "S6_R3", "S6_R4", "S6_R5", "S6_R6", "S6_R7", "S6_R8", "S6_R9")]

column_comparisons3 = [("S8_R1", "S8_R2", "S8_R3", "S8_R4", "S8_R5", "S8_R6", "S8_R7", "S8_R8", "S8_R9"),
                       ("S9_R1", "S9_R2", "S9_R3", "S9_R4", "S9_R5", "S9_R6", "S9_R7", "S9_R8", "S9_R9")]

mlac_membrane_columns = ["S3_R1", "S3_R2", "S3_R3", "S3_R4", "S3_R5", "S3_R6", "S3_R7", "S3_R8", "S3_R9"]
letAB_membrane_columns = ["S6_R1", "S6_R2", "S6_R3", "S6_R4", "S6_R5", "S6_R6", "S6_R7", "S6_R8", "S6_R9"]
pqiAB_membrane_columns = ["S9_R1", "S9_R2", "S9_R3", "S9_R4", "S9_R5", "S9_R6", "S9_R7", "S9_R8", "S9_R9"]


#t_crit = 2.306 # I think this is right for n=9, df=8, alpha=0.05, two-tailed


make_content_paired_barcharts(["CL", "PE", "PG"], df, column_comparisons2)