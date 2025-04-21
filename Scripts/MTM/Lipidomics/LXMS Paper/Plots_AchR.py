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

ranges = [[0, 50], [35, 45], [1, 5]]
ranges2 = [[0, 3], [-0.4, 0.4], [-0.25, 0.25]]

ylabels = ["Mol %", "Tail Length", "Unsaturation"]
titles = ["Ctrl LX", "Assembly", "AchR LX"]

rawdf = pd.read_excel(path, sheet_name="Intensity Stats")
# Strip off the class from the name
classes = rawdf["Class Name"]

# Get unique classes while preserving order
indexes = np.unique(classes, return_index=True)[1]
uclasses = np.array([classes[index] for index in sorted(indexes)])

# include_classes = ["PE-O"]
# uclasses = include_classes
colors = get_colors(uclasses)

for x, pval_df in enumerate(pvaldfs):

    figure = plt.figure()
    axes = figure.subplots(2, 3)
    for i, comparison in enumerate(comparisons):
        print(comparison)
        for j, c in enumerate(uclasses):
            # Get all the things out of the df
            # Find pval_df where A is comparison[0] and B is comparison[1]
            subdf = pval_df[(pval_df["A"] == comparison[1]) & (pval_df["B"] == comparison[0])]
            # Find where subdf class name is c
            subdf = subdf[subdf["Class Name"] == c]

            m1 = subdf["A_Mean"].values[0]
            s1 = subdf["A_Std"].values[0]
            m2 = subdf["B_Mean"].values[0]
            s2 = subdf["B_Std"].values[0]

            fc = subdf["FC"].values[0]
            fcsd = subdf["FC_Std"].values[0]

            diff = subdf["Diff"].values[0]
            diffsd = subdf["Diff_Std"].values[0]

            if len(subdf) == 0:
                pval = 1
            else:
                pval = subdf["Adjusted pval"].values[0]

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
    plt.savefig("LXMS_" + sample + "_" + ylabels[x] + ".png")
    plt.savefig("LXMS_" + sample + "_" + ylabels[x] + ".pdf")
plt.show()
