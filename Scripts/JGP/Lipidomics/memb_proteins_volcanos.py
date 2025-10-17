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

def make_content_barcharts(classes, df, group, columns, t_critical, plot_std=True):
    rep_contents = [np.zeros(len(columns)) for i in range(len(classes))]
    for idx,row in df.iterrows():
        name = row["Lipid Name"]
        print(name)
        lipid_class = get_class(name)
        for i in range(len(classes)):
            if lipid_class == classes[i]:
                for j in range(len(columns)):
                    if np.isnan(row[columns[j]]):
                        continue
                    else:
                        rep_contents[i][j] += row[columns[j]]

    #now average and plot the absolute contents
    means = np.zeros(len(classes))
    stds = np.zeros(len(classes))

    for i in range(len(classes)):
        means[i] = np.mean(rep_contents[i])
        stds[i] = np.std(rep_contents[i])
    #calculate the margin of error
    margins = np.zeros(len(classes))
    for i in range(len(classes)):
        margins[i] = (stds[i] / np.sqrt(len(columns))) * t_critical

    if plot_std:
        plt.bar(classes, means, yerr=stds, capsize=5, edgecolor='black')
    else:
        plt.bar(classes, means, yerr=margins, capsize=5, edgecolor='black')
    plt.title(group, fontsize=26)
    plt.show()

def make_diff_barchart(classes, df, groups, column_comparisons, t_critical, plot_std=True):
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


    diffs = np.zeros(len(classes))
    stds = np.zeros(len(classes))
    margins = np.zeros(len(classes))

    for i in range(len(classes)):
        mean1 = np.mean(class_areas1[i])
        mean2 = np.mean(class_areas2[i])
        std1 = np.std(class_areas1[i])
        std2 = np.std(class_areas2[i])

        diffs[i] = (mean1 - mean2) * 100
        stds[i] = (np.sqrt(std1**2 + std2**2))*100
        margins[i] = (stds[i] / np.sqrt(len(column_comparisons[0]))) * t_critical

    plt.title(str(groups[0]) + " vs " + str(groups[1]), fontsize=26)
    plt.xlabel("Lipid Class", fontsize=22)
    plt.ylabel("Mol % Difference", fontsize=22)
    plt.subplots_adjust(left=0.20, bottom=0.2)
    if plot_std:
        plt.bar(classes, diffs, yerr=stds, color=colors, capsize=5, edgecolor='black')
        plt.hlines(0, -1, len(classes), color='black', linestyle='-')
    else:
        plt.bar(classes, diffs, yerr=margins, color=colors, capsize=5, edgecolor='black')
        plt.hlines(0, -1, len(classes), color='black', linestyle='-')
    plt.show()

def make_fc_barchart(classes, df, groups, column_comparisons, t_critical, plot_std=True):
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

    for i in range(len(classes)):
        mean1 = np.mean(class_areas1[i])
        mean2 = np.mean(class_areas2[i])
        std1 = np.std(class_areas1[i])
        std2 = np.std(class_areas2[i])

        fcs[i] = mean1 / mean2
        stds[i] = fcs[i] * np.sqrt((std1/mean1)**2 + (std2/mean2)**2)
        margins[i] = (stds[i] / np.sqrt(len(column_comparisons[0]))) * t_critical

    #perform a paired t-test to determine if the fold change is significantly different from 1
    t_stats = np.zeros(len(classes))
    p_vals = np.zeros(len(classes))

    for i in range(len(classes)):
        t_stats[i], p_vals[i] = stats.ttest_ind(class_areas1[i], class_areas2[i])

    #correct the p-values
    rejected,corrected_pvalues,_,_ = multipletests(p_vals, method='fdr_bh')
    for i in range(len(classes)):
        if corrected_pvalues[i] < 0.05:
            print("Significant difference in " + classes[i] + " content between " + groups[0] + " and " + groups[1] + " (p = " + str(corrected_pvalues[i]) + ")")
        else:
            print("No significant difference in " + classes[i] + " content between " + groups[0] + " and " + groups[1] + " (p = " + str(corrected_pvalues[i]) + ")")


    plt.title(str(groups[0]) + " vs " + str(groups[1]), fontsize=26)
    plt.xlabel("Lipid Class", fontsize=22)
    plt.ylabel("Fold Change\nprotein/membrane", fontsize=22)
    plt.hlines(0, -1, len(classes), color='black', linestyle='-')
    plt.hlines(1, -1, len(classes), color='red', linestyle='--')
    plt.subplots_adjust(left=0.20, bottom=0.2)
    if plot_std:
        plt.bar(classes, fcs, yerr=stds, color=colors, capsize=5, edgecolor='black')
    else:
        plt.bar(classes, fcs, yerr=margins, color=colors, capsize=5, edgecolor='black')
    plt.show()

def make_volcano_plot(classes, df, groups, column_comparisons, p_cutoff, fc_cutoff):
    #Just made a random list of these, if there are more classes than this, need to add more markers...
    markers = ['o', '^', 's', 'p', 'h', 'D', 'X', 'v', '<', '>', 'd', 'P', '*', 'H', '8', '1', '2', '3', '4', 'x', '+', '|', '_']


    names = []
    p_values = []
    fcs = []

    for idx,row in df.iterrows():
        names.append(row["Lipid Name"])
        data1 = np.zeros(len(column_comparisons[0]))
        data2 = np.zeros(len(column_comparisons[1]))

        for i in range(len(column_comparisons[0])):
            data1[i] = row[column_comparisons[0][i]]

        for i in range(len(column_comparisons[1])):
            data2[i] = row[column_comparisons[1][i]]

        t_stat, p_val = stats.ttest_ind(data1, data2)
        if(np.mean(data2) == 0 or np.mean(data1) == 0):
            p_values.append(1)
            fcs.append(1)
        else:
            p_values.append(p_val)
            fcs.append(np.mean(data1)/np.mean(data2))

    rejected,corrected_pvalues,_,_ = multipletests(p_values, method='fdr_bh')
    log2_fcs = np.log2(fcs)
    log10_pvals = -np.log10(corrected_pvalues)

    class_separated_names = [[] for i in range(len(classes))]
    class_separated_log2_fcs = [[] for i in range(len(classes))]
    class_separated_log10_pvals = [[] for i in range(len(classes))]

    for i in range(len(names)):
        for j in range(len(classes)):
            if get_class(names[i]) == classes[j]:
                class_separated_names[j].append(names[i])
                class_separated_log2_fcs[j].append(log2_fcs[i])
                class_separated_log10_pvals[j].append(log10_pvals[i])

    #plot each class separately and label each data point with the lipid name
    for i in range(len(classes)):
        plt.scatter(class_separated_log2_fcs[i], class_separated_log10_pvals[i], label=classes[i], marker=markers[i])
        for j in range(len(class_separated_names[i])):
            plt.text(class_separated_log2_fcs[i][j], class_separated_log10_pvals[i][j], class_separated_names[i][j], fontsize=8)

    plt.vlines(-np.log2(fc_cutoff), 0, 20, color='r', linestyle='--')
    plt.vlines(np.log2(fc_cutoff), 0, 20, color='r', linestyle='--')
    plt.hlines(-np.log10(p_cutoff), -8, 8, color='r', linestyle='--')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.xlabel("Log2 Fold Change", fontsize=22)
    plt.ylabel("-Log10 p-value", fontsize=22)
    plt.legend(fontsize=16)
    plt.show()

def make_tail_length_diff_barchart(df, groups, column_comparisons, t_critical, plot_std=True, classes=['CL', 'PE', 'PG']):
    '''
    Makes a bar chart of the difference in lipid class composition between each comparison group.
    :param names: Names of the lipids
    :param classes: Classes of the lipids to be compared
    :param df: DataFrame containing the lipid data ("Lipid Name" column must be present)
    :param column_comparisons: Two lists of strings (column names) to be compared
    '''

    totals1 = [np.zeros(len(column_comparisons[0])) for i in range(len(classes))]
    tls_1 = [np.zeros(len(column_comparisons[0])) for i in range(len(classes))]
    totals2 = [np.zeros(len(column_comparisons[1])) for i in range(len(classes))]
    tls_2 = [np.zeros(len(column_comparisons[1])) for i in range(len(classes))]

    for idx,row in df.iterrows():
        name = row["Lipid Name"]
        lipid_class = get_class(name)
        tname, tl, tu = get_overall_tails(name)

        for c in range(len(classes)):
            if lipid_class == classes[c]:
                for i in range(len(column_comparisons[0])):
                    tls_1[c][i] += row[column_comparisons[0][i]] * tl
                    totals1[c][i] += row[column_comparisons[0][i]]

                for i in range(len(column_comparisons[1])):
                    tls_2[c][i] += row[column_comparisons[1][i]] * tl
                    totals2[c][i] += row[column_comparisons[1][i]]


    for c in range(len(classes)):
        for i in range(len(column_comparisons[0])):
            tls_1[c][i] = tls_1[c][i] / totals1[c][i]
        for i in range(len(column_comparisons[1])):
            tls_2[c][i] = tls_2[c][i] / totals2[c][i]

    tls_1 = np.array(tls_1)
    tls_2 = np.array(tls_2)

    means1 = np.zeros(len(classes))
    means2 = np.zeros(len(classes))
    diffs = np.zeros(len(classes))
    diff_stds = np.zeros(len(classes))
    margins = np.zeros(len(classes))

    for c in range(len(classes)):
        means1[c] = np.mean(tls_1[c,:])
        means2[c] = np.mean(tls_2[c,:])
        diffs[c] = means1[c] - means2[c]
        std1 = np.std(tls_1[c,:])
        std2 = np.std(tls_2[c,:])
        diff_stds[c] = np.sqrt(std1**2 + std2**2)
        margins[c] = (diff_stds[c] / np.sqrt(len(column_comparisons[0]))) * t_critical

    plt.title(str(groups[0]) + " vs " + str(groups[1]), fontsize=26)
    plt.xlabel("Lipid Class", fontsize=22)
    plt.ylabel("Tail Length Difference", fontsize=22)
    plt.hlines(0, -1, len(classes), color='black', linestyle='-')
    plt.subplots_adjust(left=0.20, bottom=0.2)
    if plot_std:
        plt.bar(classes, diffs, yerr=diff_stds, capsize=5, edgecolor='black')
    else:
        plt.bar(classes, diffs, yerr=margins, capsize=5, edgecolor='black')
    plt.show()

def make_content_paired_barcharts(classes, df, group, column_comparisons, t_critical, plot_std=True):
    colors = get_colors(classes)

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

    for i in range(len(classes)):
        group1_means[i] = np.mean(group1_contents[i]*100)
        group2_means[i] = np.mean(group2_contents[i]*100)
        group1_stds[i] = np.std(group1_contents[i]*100)
        group2_stds[i] = np.std(group2_contents[i]*100)
        group1_margins[i] = (group1_stds[i] / np.sqrt(len(column_comparisons[0])))*t_critical
        group2_margins[i] = (group2_stds[i] / np.sqrt(len(column_comparisons[1])))*t_critical


    #Now perform a t-test to determine if the means are significantly different
    t_stats = np.zeros(len(classes))
    p_vals = np.zeros(len(classes))

    for i in range(len(classes)):
        t_stats[i], p_vals[i] = stats.ttest_ind(group1_contents[i], group2_contents[i])


    bar_width = 0.35
    spacing = 0.05

    r1 = np.arange(len(classes))
    r2 = [x + bar_width + spacing for x in r1]

    if plot_std:
        plt.bar(r1, group1_means, yerr=group1_stds, color=colors, capsize=5, edgecolor='black',
                width=bar_width, label='LetAB Protein')
        plt.bar(r2, group2_means, yerr=group2_stds, color=colors, capsize=5, edgecolor='black',
                hatch='\\', width=bar_width, label='LetAB Membrane')
    else:
        plt.bar(r1, group1_means, yerr=group1_margins, color=colors, capsize=5, edgecolor='black',
                width=bar_width, label='LetAB Protein')
        plt.bar(r2, group2_means, yerr=group2_margins, color=colors, capsize=5, edgecolor='black',
                hatch='\\', width=bar_width, label='LetAB Membrane')

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
    plt.xticks([r + bar_width/2 for r in range(len(classes))], classes, fontsize=18)
    plt.xlabel('Lipid Class', fontsize=22)
    plt.ylim(0, 125)
    plt.yticks(np.arange(0, 110, 50), fontsize=18)
    plt.ylabel('Summed Intensity %', fontsize=22)
    plt.tight_layout()
    plt.legend(handles=legend_handles, fontsize=14,  loc='upper right')
    plt.savefig("Z:\\Group Share\\JGP\\CS MS1 Runs\\mzML\\letAB_paired_barcharg.pdf", dpi=300)
    plt.savefig("Z:\\Group Share\\JGP\\CS MS1 Runs\\mzML\\letAB_paired_barcharg.png", dpi=300)
    plt.show()

def perform_pca(classes, df, columns):
    class_total_areas = [np.zeros(len(columns)) for i in range(len(classes))]
    for idx,row in df.iterrows():
        name = row["Lipid Name"]
        lipid_class = get_class(name)
        for i in range(len(classes)):
            if lipid_class == classes[i]:
                for j in range(len(columns)):
                    if np.isnan(row[columns[j]]):
                        continue
                    else:
                        class_total_areas[i][j] += row[columns[j]]

    samples = ["LetAB Protein1", "LetAB Protein2", "LetAB Protein3", "LetAB Protein4", "LetAB Protein5", "LetAB Protein6", "LetAB Protein7", "LetAB Protein8",
               "LetAB Membrane1", "LetAB Membrane2", "LetAB Membrane3", "LetAB Membrane4", "LetAB Membrane5", "LetAB Membrane6", "LetAB Membrane7", "LetAB Membrane8", "LetAB Membrane9"]
    total_areas_df = pd.DataFrame({"CL": class_total_areas[0], "PE": class_total_areas[1], "PG": class_total_areas[2]})

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(total_areas_df)

    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_data)

    pca_df = pd.DataFrame(data=principal_components, columns=["PC1", "PC2"])
    print(pca_df)
    print(f'Explained Variance Ratio: {pca.explained_variance_ratio_}')

    legend_handles = [
        Line2D([0], [0], marker='o', color='w', label='LetAB Protein', markerfacecolor='blue', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='LetAB Membrane', markerfacecolor='red', markersize=10)
    ]

    #Label each data point with a sample name
    for i in range(len(samples)):
        if i >= 8:
            plt.scatter(pca_df["PC1"][i], pca_df["PC2"][i], label=samples[i], color='red')
        else:
            plt.scatter(pca_df["PC1"][i], pca_df["PC2"][i], label=samples[i], color='blue')
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(handles=legend_handles)
    plt.show()

# Set backend to Agg
mpl.use('WxAgg')

file_path = "Z:\\Group Share\\Annika\\CS Collab Data\\Collaborator Samples (CS)\\lipidomic_analysis_results.xlsx"

df = pd.read_excel(r"Z:\Group Share\JGP\Christina_letAB_MlaC_pqiAB_data\Analysis_and_Other\normalized_abundances.xlsx")
df = pd.read_excel(file_path, sheet_name="Aggregated")
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




#make_content_barcharts(["CL", "PE", "PG"], df, "MlaC Protein", letAB_membrane_columns, 2.306, plot_std=True)

t_crit = 2.306
plot_std = False

# perform_pca(["CL", "PE", "PG"], df,
#             ["S2_R1", "S2_R2", "S2_R3", "S2_R4", "S2_R5", "S2_R6", "S2_R7", "S2_R8",
#              "S3_R1", "S3_R2", "S3_R3", "S3_R4", "S3_R5", "S3_R6", "S3_R7", "S3_R8", "S3_R9"])




# make_content_barcharts(["CL", "PE", "PG"], df, "LetAB Membrane", letAB_membrane_columns, t_crit, plot_std=plot_std)
# make_content_barcharts(["CL", "PE", "PG"], df, "MlaC Membrane", mlac_membrane_columns, t_crit, plot_std=plot_std)
# make_content_barcharts(["CL", "PE", "PG"], df, "pqiAB Membrane", pqiAB_membrane_columns, t_crit, plot_std=plot_std)
#
# make_content_paired_barcharts(["CL", "PE", "PG"], df, "MlaC Protein", column_comparisons2, t_crit, plot_std=plot_std)
#
#
# #Plot MlaC
# # make_diff_barchart(["CL", "PE", "PG"], df, ["MlaC Protein","MlaC Membrane"],
# #                    column_comparisons1, t_crit, plot_std=plot_std)
#
# make_fc_barchart(["CL", "PE", "PG"], df, ["MlaC Protein","MlaC Membrane"],
#                  column_comparisons1, t_crit, plot_std=plot_std)

#Plot LetAB
# make_diff_barchart(["CL", "PE", "PG"], df, ["LetAB Protein","LetAB Membrane"],
#                    column_comparisons2, t_crit, plot_std=plot_std)

make_fc_barchart(["CL", "PE", "PG"], df, ["LetAB Protein","LetAB Membrane"],
                 column_comparisons2, t_crit, plot_std=plot_std)

exit()
#plot pqiAB
make_diff_barchart(["CL", "PE", "PG"], df, ["pqiAB Protein","pqiAB Membrane"],
                   column_comparisons3, t_crit, plot_std=plot_std)

make_fc_barchart(["CL", "PE", "PG"], df, ["pqiAB Protein","pqiAB Membrane"],
                 column_comparisons3, t_crit, plot_std=plot_std)


make_volcano_plot(["CL", "PE", "PG"], df, ["MlaC Protein","MlaC Membrane"],
                  column_comparisons1,  0.05, 2)

make_volcano_plot(["CL", "PE", "PG"], df, ["MlaC Protein","MlaC Membrane"],
                  column_comparisons2,  0.05, 2)

make_tail_length_diff_barchart(df, ["MlaC Protein","MlaC Membrane"],
                               column_comparisons1, t_crit, plot_std=plot_std)

make_tail_length_diff_barchart(df, ["LetAB Protein","LetAB Membrane"],
                               column_comparisons2, t_crit, plot_std=plot_std)

exit()

