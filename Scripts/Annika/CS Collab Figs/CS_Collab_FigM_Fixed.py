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
