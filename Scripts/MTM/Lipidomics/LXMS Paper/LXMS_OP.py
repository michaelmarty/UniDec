import matplotlib.pyplot as plt
import os

plt.rcParams.update({'font.size': 11})

sample = "AchR"
os.chdir("C:\Data\Lipidomics\\AchR")

comparisons = [["FT_Ctrl", "Elu_Ctrl"],
               ["E3(+)_Ctrl", "E3(+)_AchR"],
               ["FT_AchR", "Elu_AchR"],
               ]

# Set Paths
path = "LXMS_" + sample + "_Classes.xlsx"

pvaldfs = [
    pd.read_excel(path, sheet_name="Intensity Stats"),
    pd.read_excel(path, sheet_name="Tail Length Stats"),
    pd.read_excel(path, sheet_name="Tail Unsaturation Stats")
]
pvaldfs = pvaldfs[0:1]

ranges = [[0, 80], [35, 45], [1, 5]]
ranges2 = [[0, 3], [-0.5, 0.5], [-0.5, 0.5]]

ylabels = ["Mol % P/O", "Tail Length", "Unsaturation"]
titles = ["Ctrl LX", "Assembly", "AchR LX"]

uclasses = ["PC", "PE", "PI", "PS"]
colors = get_colors(uclasses)


for x, pval_df in enumerate(pvaldfs):

    figure = plt.figure()
    axes = figure.subplots(2, 3)
    for i, comparison in enumerate(comparisons):
        # Get all the things out of the df
        # Find pval_df where A is comparison[0] and B is comparison[1]
        subdf1 = pval_df[(pval_df["A"] == comparison[1]) & (pval_df["B"] == comparison[0])]
        print(comparison)
        percent_po = []
        for j, c in enumerate(uclasses):
            c2 = c + "-O/P"
            # Find where subdf class name is c
            subdf = subdf1[subdf1["Class Name"] == c]
            subdf2 = subdf1[subdf1["Class Name"] == c2]

            m1 = subdf["A_Mean"].values[0]
            s1 = subdf["A_Std"].values[0]
            m2 = subdf["B_Mean"].values[0]
            s2 = subdf["B_Std"].values[0]

            m1_2 = subdf2["A_Mean"].values[0]
            s1_2 = subdf2["A_Std"].values[0]
            m2_2 = subdf2["B_Mean"].values[0]
            s2_2 = subdf2["B_Std"].values[0]

            po_per_1 = 100*m1_2 / (m1 + m1_2)
            po_per_2 = 100*m2_2 / (m2 + m2_2)
            po_per_1sd = po_per_1 * np.sqrt((s1_2 / m1_2) ** 2 + (s1 / m1) ** 2)
            po_per_2sd = po_per_2 * np.sqrt((s2_2 / m2_2) ** 2 + (s2 / m2) ** 2)

            m1 = po_per_1
            s1 = po_per_1sd
            m2 = po_per_2
            s2 = po_per_2sd

            pval = 1

            axis = axes[0][i]
            axis.bar(j, m1, color=colors[j], width=0.45, yerr=s1, ecolor="black")
            axis.bar(j + 0.45, m2, color=colors[j], width=0.45, yerr=s2, ecolor="black")

            max_mean = max(m1 + s1, m2 + s2)
            offset = (ranges[x][1] - ranges[x][0]) * 0.05
            if pval < 0.05:
                axis.plot([j, j + 0.45], [max_mean + offset, max_mean + offset], color="black")

            if j == 0:
                axis.set_ylim(ranges[x][0], ranges[x][1])
                if i == 0:
                    axis.set_ylabel(ylabels[x])
                else:
                    axis.set_yticklabels([])
                axis.set_xticks(np.arange(len(uclasses)) + 0.225, uclasses, rotation=90)
                # REmove the top and right lines
                axis.spines['right'].set_visible(False)
                axis.spines['top'].set_visible(False)
                axis.set_title(titles[i])
            continue
            axis = axes[1][i]
            if x == 0:
                # Plot fold change
                axis.bar(j, fc, color=colors[j], width=0.8, yerr=fcsd, ecolor="black")
                if pval < 0.1:
                    marker = "."
                    if pval < 0.05:
                        marker = "*"
                    # Add an asterisk above the bar
                    ypos = np.max([1.05, fc + fcsd + 0.1])
                    axis.text(j, ypos, marker, horizontalalignment='center',
                              verticalalignment='center')
                if j == 0:
                    if i == 0:
                        axis.set_ylabel("Fold Change")
                    else:
                        axis.set_yticklabels([])

            else:
                axis.bar(j, diff, color=colors[j], width=0.8, yerr=diffsd, ecolor="black")
                if pval < 0.1:
                    marker = "."
                    if pval < 0.05:
                        marker = "*"
                    # Add an asterisk above the bar
                    ypos = np.max([0.05, diff + diffsd + 0.05])
                    axis.text(j, ypos, marker, horizontalalignment='center',
                              verticalalignment='center')
                if j == 0:
                    if i == 0:
                        axis.set_ylabel("Difference")
                    else:
                        axis.set_yticklabels([])

            if j == 0:
                axis.set_ylim(ranges2[x][0], ranges2[x][1])

                axis.set_xticks(np.arange(len(uclasses)), uclasses, rotation=90)
                # REmove the top and right lines
                axis.spines['right'].set_visible(False)
                axis.spines['top'].set_visible(False)

        if x == 0:
            axis.plot([0, len(uclasses) - 1], [1, 1], color="black", linestyle="--")
        else:
            axis.plot([0, len(uclasses) - 1], [0, 0], color="black", linestyle="--")

    plt.tight_layout()
plt.show()
