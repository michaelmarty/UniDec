import matplotlib.pyplot as plt
import os
from scipy import stats
# Set default plot font size to 11
plt.rcParams.update({'font.size': 10})
# Set font to ariel
plt.rcParams["font.family"] = "Arial"

path = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\Papers\\LX-MS\\230517_LX_Supplemental_Data.xlsx"

df = pd.read_excel(path, sheet_name="Table S4 - Ecoli Lipid Results", skiprows=1)

os.chdir(os.path.dirname(path))

# Reset class names to remove O/P
df["Class"] = df["Class"].str.replace(" O/P", "")

classes = df["Class"].values
uclasses = np.unique(classes)
print(uclasses)
#uclasses = ["TG"]
#uclasses = ["PI", "PS"]

# Limit DF to only the classes in uclasses
df = df[df["Class"].isin(uclasses)]

colors = get_colors(uclasses)
print(colors)

comps = ["FT_vs_Elu", "AmtB_vs_Ctrl", "AqpZ_vs_Ctrl", "FT_vs_Elu.1", "FT_vs_Elu.2"]
titles = ["Control Exchange", "AmtB Assembly", "AqpZ Assembly", "AmtB LX", "AqpZ LX"]

pcutoff = 0.05
fcutoff = 0.5
markersize = 3
xlim = 1

n_annulus = 45 # Ellena JF, Blazing MA, McNamee MG. Lipid-protein interactions in reconstituted membranes containing acetylcholine receptor. Biochemistry. 1983;22(24):5523-35. doi: 10.1021/bi00293a012.
n_disc = 200
n_empty = 250

bulkratio = (n_disc-n_annulus) / n_disc
print(bulkratio)

#bulkratio = 0.5

R = 8.31446261815324 / 1000.
units = "kJ/mol"
# R = 1.98720425864083 / 1000.
# units = "kcal/mol"
T = 22 + 273.15

dg = False

for i, c in enumerate(comps):
    fcname = "Fold_change_" + c
    pvalname = "P_value_A_" + c

    figure = plt.figure(figsize=(3.25, 3))

    axis = figure.gca()

    comp1, comp2 = c.split("_vs_")

    for j, u in enumerate(uclasses):
        subdf = df[df["Class"] == u]

        mean1 = subdf["Mean_"+comp1].values
        mean2 = subdf["Mean_"+comp2].values

        mean_12 = np.mean([mean1, mean2], axis=0)

        fc = subdf[fcname].values
        pval = subdf[pvalname].values

        logp = -np.log10(pval)

        b1 = pval < pcutoff
        b1fc = np.abs(fc) > fcutoff
        b3 = np.logical_and(b1, b1fc)

        # b2 is not b1
        b2 = np.logical_not(b1)

        results = stats.linregress(mean_12, fc)
        slope = results[0]
        intercept = results[1]
        rvalue = results[2]
        pvalue = results[3]


        if pvalue < 0.05:
            print("Significant")
            print(c, u, slope, intercept, rvalue, pvalue)

        axis.plot(mean_12, fc, "o", label=u, color=colors[j], markersize=markersize)

        #axis.plot(mean_12[b1], fc[b1], "o", label=u, color=colors[j], markersize=markersize)
        #axis.plot(mean_12[b2], fc[b2], "o", color="grey", markersize=markersize)

        # Draw dashed line at pcutoff
        #axis.plot([-xlim, xlim], [-np.log10(pcutoff), -np.log10(pcutoff)], "--", color="k")
        # Draw dashed line a fcutoff
        #axis.plot([fcutoff, fcutoff], [0, 6], "--", color="k")
        #axis.plot([-fcutoff, -fcutoff], [0, 6], "--", color="k")
        # axis.plot([0, 0], [0, 6], "-", color="k")

        #axis.set_ylim([0, 6])
        #axis.set_xlim([-xlim, xlim])

        axis.set_title(titles[i])

        # Turn off top axis
        axis.spines['right'].set_visible(False)
        axis.spines['top'].set_visible(False)

        # Set y axis to cross at x=0
        # axis.spines['left'].set_position(('data', 0))

        axis.set_ylabel("Mean Abundance")

        axis.set_ylabel("Log2 Fold Change")
    plt.legend()
    #plt.tight_layout()

    #plt.savefig(c + ".pdf")
plt.show()
