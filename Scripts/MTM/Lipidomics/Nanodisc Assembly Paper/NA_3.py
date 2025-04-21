import os
from scipy import stats
import matplotlib.pyplot as plt
from Lipidomics_Stats_Tests import *

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

os.chdir("C:\Data\Lipidomics\\NanodiscAssembly")

# Set Paths
sample = "Ecoli"
path = "NA_" + sample + "_Classes.xlsx"

# Load the data
df = pd.read_excel(path, sheet_name="Class Means")

rawdf = pd.read_excel(path, sheet_name="Classes")

# Strip off the class from the name
classes = df["Class Name"]

# Get unique classes while preserving order
indexes = np.unique(classes, return_index=True)[1]
uclasses = np.array([classes[index] for index in sorted(indexes)])

# include_classes = ["PE-O"]
# uclasses = include_classes
colors = get_colors(uclasses)

addendum = "_" + sample

comparisons = [["37_D1", "37_E3"],
               ["25_D1", "25_E3"],
               ["4_D1", "4_E3"],
               ["4_E3", "37_E3"],
               ["25_E3", "37_E3"],
               ["4_E3", "25_E3"],
               ["4_D1", "37_D1"],
               ["25_D1", "37_D1"],
               ["4_D1", "25_D1"],
               ]

pdf = pd.DataFrame(uclasses, columns=["Class Name"])

for sample_types in comparisons:

    figure = plt.figure()

    # Significance test the difference between the two
    comparison = sample_types
    print(comparison)
    pvals = []

    column_names1 = [col for col in rawdf.columns if comparison[0] in col]
    column_names2 = [col for col in rawdf.columns if comparison[1] in col]

    df1 = rawdf[column_names1]
    df2 = rawdf[column_names2]

    m1 = np.mean(df1, axis=1)
    s1 = np.std(df1, axis=1)
    m2 = np.mean(df2, axis=1)
    s2 = np.std(df2, axis=1)
    for j, c in enumerate(uclasses):
        plt.bar(j, m1[j] * 100, color=colors[j], width=0.45, yerr=s1[j] * 100, ecolor="black")
        plt.bar(j + 0.45, m2[j] * 100, color=colors[j], width=0.45, yerr=s2[j] * 100, ecolor="black")

    for j, c in enumerate(uclasses):

        df1 = rawdf[rawdf["Class Name"] == c][column_names1]
        df2 = rawdf[rawdf["Class Name"] == c][column_names2]
        pval = stats.ttest_ind(df1.to_numpy()[0], df2.to_numpy()[0])[1]
        print(c, pval)

        max_mean = max(max(df1.to_numpy()[0] * 100), max(df2.to_numpy()[0] * 100))
        if pval < 0.05:
            plt.plot([j, j + 0.45], [max_mean + 1, max_mean + 1], color="black")
        pvals.append(pval)

    title = comparison[0] + " vs " + comparison[1] + " " + sample

    pdf[title + " pvals"] = pvals
    plt.title(title)
    plt.ylim(0, 50)
    plt.xticks(np.arange(len(uclasses)) + 0.225, uclasses, rotation=90)
    # REmove the top and right lines
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    # Set the font size to 11
    plt.rcParams.update({'font.size': 11})
    plt.ylabel("Mean Relative Intensity (%)")
    plt.tight_layout()
    plt.savefig(title + ".png", dpi=300)
    # plt.show()

writer = pd.ExcelWriter(path, engine='openpyxl', mode='a', if_sheet_exists="replace")
pdf.to_excel(writer, index=False, sheet_name="Class Comparison Pvals")
writer.close()
