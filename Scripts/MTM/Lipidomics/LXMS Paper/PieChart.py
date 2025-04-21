import os
import matplotlib.pyplot as plt

if True:
    os.chdir("C:\Data\Lipidomics\\AchR")

    # Set Paths
    path = "LXMS_AchR_Output.xlsx"

    outfile = "LXMS_AchR_Classes.xlsx"

    drop_classes = ["LPC", "LPE", "LPE-O", "LPE-O/P", "TG", "PG", "Cer", "CL"]
    sample_types = ["FT_Ctrl", "Elu_Ctrl", "FT_AchR", "Elu_AchR", "E3(+)_Ctrl", "E3(+)_AchR"]
    outliers = ["FT_Ctrl_1"]
    outliers = []

    comparisons = [["FT_Ctrl", "Elu_Ctrl"],
                   ["E3(+)_Ctrl", "E3(+)_AchR"],
                   ["FT_AchR", "Elu_AchR"],
                   ]

    n_annulus = 45  # Ellena JF, Blazing MA, McNamee MG. Lipid-protein interactions in reconstituted membranes containing acetylcholine receptor. Biochemistry. 1983;22(24):5523-35. doi: 10.1021/bi00293a012.
    n_disc = 200  # CharmmGUI
    n_empty = 240
    titles = ["Assembly", "AchR LX"]
    figure_save = "AchR_PieChart.pdf"
else:
    os.chdir("C:\Data\Lipidomics\\AmtBAqpZ")
    path = "LXMS_Ecoli_Output.xlsx"

    outfile = "LXMS_Ecoli_Classes.xlsx"

    drop_classes = ["LPE", "LPE-O", "LPE-O/P", ]

    sample_types = ["Elu_Ctrl", "FT_Ctrl", "Elu_AmtB", "FT_AmtB", "E3_Ctrl", "E3_AmtB",
                    "Elu_AqpZ", "FT_AqpZ", "E3_Ctrl", "E3_AqpZ"]
    outliers = []

    comparisons = [["Elu_Ctrl", "FT_Ctrl"],
                   ["E3_Ctrl", "E3_AmtB"],
                   ["Elu_AmtB", "FT_AmtB"],
                   ["E3_Ctrl", "E3_AqpZ"],
                   ["Elu_AqpZ", "FT_AqpZ"],
                   ]

    n_annulus = 75  # Marty, 2016 Ange. Chem.
    n_disc = 150  # Keener, 2019 JACS
    n_empty = 240
    titles = ["Assembly", "AmtB LX", "Assembly", "AqpZ LX"]
    figure_save = "AmtBAqpZ_PieChart.pdf"

df = pd.read_excel(outfile, sheet_name="Class Means")
print(df.columns)

bulkratio = (n_disc - n_annulus) / n_disc
print(bulkratio)

a = 1 / (1 - bulkratio)
b = bulkratio / (1 - bulkratio)

classes = df["Class Name"].values
colors = get_colors(classes)
colors = np.array(colors)
plt.figure()


for i, c in enumerate(comparisons[1:]):
    m1 = df[c[0] + " Mean"].values
    s1 = df[c[0] + " Std"].values
    m2 = df[c[1] + " Mean"].values
    s2 = df[c[1] + " Std"].values

    xa = a * m2 - b * m1

    # sort xa
    sort = np.argsort(xa)
    xa = xa[sort]
    m1 = m1[sort]
    m2 = m2[sort]




    xas = np.sqrt(a ** 2 * s2 ** 2 + b ** 2 * s1 ** 2)

    print(c[0], c[1], xa, xas)

    b2 = xa > 0.01
    b3 = m1 > 0.01
    b4 = m2 > 0.01




    # Make a Pie Chart
    plt.subplot(3, len(comparisons) - 1, i + 1 + (len(comparisons) - 1))
    plt.title(titles[i] + " ND", fontdict={"fontsize": 10})
    plt.pie(m2[b4], labels=classes[sort][b4], autopct='%1.0f%%', shadow=True, startangle=90,
            colors=colors[sort][b4])

    plt.subplot(3, len(comparisons)-1, i+1 + 2* (len(comparisons)-1))
    plt.title(titles[i] + " Annulus", fontdict={"fontsize": 10})
    plt.pie(xa[b2], labels=classes[sort][b2], autopct='%1.0f%%', shadow=True, startangle=90,
            colors=colors[sort][b2])

    plt.subplot(3, len(comparisons)-1, i+1)
    plt.title(titles[i] + " Bulk", fontdict={"fontsize": 10})
    plt.pie(m1[b3], labels=classes[sort][b3], autopct='%1.0f%%', shadow=True, startangle=90,
            colors=colors[sort][b3])
plt.tight_layout()

plt.savefig(figure_save)

plt.show()
