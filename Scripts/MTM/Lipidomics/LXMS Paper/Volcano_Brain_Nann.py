import matplotlib.pyplot as plt
import os

# Set default plot font size to 11
plt.rcParams.update({'font.size': 10})
# Set font to ariel
plt.rcParams["font.family"] = "Arial"

path = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\Papers\\LX-MS\\LX_Supplemental_Data_Final.xlsx"

df = pd.read_excel(path, sheet_name="Table S5 - Brain Lipid Results", skiprows=1)

os.chdir(os.path.dirname(path))

# Reset class names to remove O/P
df["Class"] = df["Class"].str.replace(" O/P", "")

classes = df["Class"].values
uclasses = np.unique(classes)
print(uclasses)
# uclasses = ["TG"]
# uclasses = ["PI", "PS"]

# Move SM in the uclasses array to after HexCer
if "SM" in uclasses and "HexCer" in uclasses and "Cer" in uclasses:
    uclasses = np.delete(uclasses, np.where(uclasses == "SM"))
    uclasses = np.delete(uclasses, np.where(uclasses == "Cer"))
    insert_index = np.where(uclasses == "HexCer")[0][0] + 1
    uclasses = np.insert(uclasses, insert_index, "Cer")
    uclasses = np.insert(uclasses, insert_index + 1, "SM")

# Move CL to the back
if "CL" in uclasses:
    uclasses = np.delete(uclasses, np.where(uclasses == "CL"))
    uclasses = np.insert(uclasses, len(uclasses), "CL")

# Move TG to the back
if "TG" in uclasses:
    uclasses = np.delete(uclasses, np.where(uclasses == "TG"))
    uclasses = np.insert(uclasses, len(uclasses), "TG")

# Limit DF to only the classes in uclasses
df = df[df["Class"].isin(uclasses)]

colors = get_colors(uclasses)
print(colors)

symbols = get_symbols(uclasses)

comps = ["AChR_vs_Ctrl", "Elu_vs_FT.1"]
titles = ["AChR Assembly", "AChR LX"]

pcutoff = 0.05
fcutoff = 45
markersize = 3
xlim = 200

n_annulus = 45  # Ellena JF, Blazing MA, McNamee MG. Lipid-protein interactions in reconstituted membranes containing acetylcholine receptor. Biochemistry. 1983;22(24):5523-35. doi: 10.1021/bi00293a012.
n_disc = 200  # CharmmGUI
n_empty = 240

bulkratio = (n_disc - n_annulus) / n_disc
fcutoff2 = 45

figure = plt.figure(figsize=(6.5, 3))
axes = figure.subplots(1, 2, sharex=True, sharey=True)

for i, c in enumerate(comps):
    fcname = "Fold_change_" + c
    pvalname = "P_value_A_" + c

    axis = axes[i]

    total_lipids = len(df)
    b1 = df[pvalname] < pcutoff
    b2 = np.abs(df[fcname]) > fcutoff
    b3 = np.logical_and(b1, b2)
    sig_lipids = len(df[b3])
    sdf = df[b3]
    sdf = sdf.sort_values(by=pvalname)
    #print(c, sig_lipids, total_lipids, sig_lipids / total_lipids * 100)

    for j, u in enumerate(uclasses):
        subdf = df[df["Class"] == u]

        fc = subdf[fcname].values
        pval = subdf[pvalname].values

        logp = -np.log10(pval)

        b1 = pval < pcutoff
        b1fc = np.abs(fc) > fcutoff
        b3 = np.logical_and(b1, b1fc)

        # b2 is not b1
        b2 = np.logical_not(b1)

        x = 2 ** fc
        # print(fc2)
        n_ann_corr = (1-x) * n_disc

        data = n_ann_corr

        #print("MINMAX:", np.min(data), np.max(data))

        names = subdf["Name"].values[b3]

        axis.plot(data[b1], logp[b1], marker=symbols[j], linestyle=" ", label=u, color=colors[j], markersize=markersize)
        axis.plot(data[b2], logp[b2], marker=symbols[j], linestyle=" ", color="grey", markersize=markersize)

    # Draw dashed line at pcutoff
    axis.plot([-xlim, xlim], [-np.log10(pcutoff), -np.log10(pcutoff)], "--", color="k")
    # Draw dashed line a fcutoff

    axis.plot([fcutoff2, fcutoff2], [0, 6], "--", color="k")
    axis.plot([-fcutoff2, -fcutoff2], [0, 6], "--", color="k")
    axis.plot([200, 200], [0, 6], "--", color="k")

    axis.set_ylim([0, 6])
    axis.set_xlim([0, xlim])

    axis.set_title(titles[i])

    # Turn off top axis
    axis.spines['right'].set_visible(False)
    axis.spines['top'].set_visible(False)

    # Set y axis to cross at x=0
    # axis.spines['left'].set_position(('data', 0))
    if i == 0:
        axis.set_ylabel("–Log10 P-Value")

    axis.set_xlabel("Effective # of Annular Lipids")

    #if i == 0:
    #    axis.legend()
plt.tight_layout()

plt.savefig("Brain_Nann.pdf")
plt.show()
