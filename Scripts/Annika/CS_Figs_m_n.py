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
    zscore = np.zeros(len(classes))
    # df = 7
    # t_critical = stats.t.ppf(1 - 0.025, df)
    confidence_level = 0.95
    z_critical = stats.norm.ppf(1 - (1 - confidence_level) / 2)
    # print(t_critical)

    p_vals = np.zeros(len(classes))
    for i in range(len(classes)):
        n = len(column_comparisons[0])
        mean1 = np.mean(class_areas1[i])
        mean2 = np.mean(class_areas2[i])
        std1 = np.std(class_areas1[i], ddof=1) / np.sqrt(n)
        std2 = np.std(class_areas2[i], ddof=1) / np.sqrt(n)

        fcs[i] = mean1 / mean2
        stds[i] = fcs[i] * np.sqrt((std1/mean1)**2 + (std2/mean2)**2)

        # margins[i] = (stds[i] / np.sqrt(n)) * t_critical
        #margins[i] = stds[i] * t_critical
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


    plt.title(str(groups[0]) + " vs " + str(groups[1]), fontsize=26)
    plt.xlabel("Lipid Class", fontsize=22)
    plt.ylabel("Fold Change\nprotein/membrane", fontsize=22)
    plt.hlines(0, -1, len(classes), color='black', linestyle='-')
    plt.hlines(1, -1, len(classes), color='red', linestyle='--')
    plt.subplots_adjust(left=0.20, bottom=0.2)
    plt.bar(classes, fcs, yerr=margins, color=colors, capsize=5, edgecolor='black')
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

