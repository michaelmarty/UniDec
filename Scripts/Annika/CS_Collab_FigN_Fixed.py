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



