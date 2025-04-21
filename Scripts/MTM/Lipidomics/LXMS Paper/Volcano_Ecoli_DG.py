import matplotlib.pyplot as plt
import os

# Set default plot font size to 11
plt.rcParams.update({'font.size': 10})
# Set font to ariel
plt.rcParams["font.family"] = "Arial"

path = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\Papers\\LX-MS\\230517_LX_Supplemental_Data_2.xlsx"

df = pd.read_excel(path, sheet_name="Table S4 - Ecoli Lipid Results", skiprows=1)

os.chdir(os.path.dirname(path))

# Reset class names to remove O/P
df["Class"] = df["Class"].str.replace(" O/P", "")

classes = df["Class"].values
uclasses = np.unique(classes)
print(uclasses)
# uclasses = ["TG"]
# uclasses = ["PI", "PS"]

# Move CL to the back
if "CL" in uclasses:
    uclasses = np.delete(uclasses, np.where(uclasses == "CL"))
    uclasses = np.insert(uclasses, len(uclasses), "CL")

# Limit DF to only the classes in uclasses
df = df[df["Class"].isin(uclasses)]

colors = get_colors(uclasses)
print(colors)
symbols = get_symbols(uclasses)
print(symbols)

comps = ["FT_vs_Elu", "AmtB_vs_Ctrl", "AqpZ_vs_Ctrl", "FT_vs_Elu.1", "FT_vs_Elu.2"]
titles = ["Control Exchange", "AmtB Assembly", "AqpZ Assembly", "AmtB LX", "AqpZ LX"]

pcutoff = 0.05
fcutoff = 0.5
markersize = 3
xlim = 4

n_annulus = 75  # Marty, 2016 Ange. Chem.
n_disc = 150  # Keener, 2019 JACS
n_empty = 240

bulkratio = (n_disc - n_annulus) / n_disc
print(bulkratio)

# bulkratio = 0.5

R = 8.31446261815324 / 1000.
units = "kJ/mol"
# R = 1.98720425864083 / 1000.
# units = "kcal/mol"
T = 22 + 273.15

fcutoff1 = - R * T * np.log(2 ** fcutoff)
print("FCUTOFF:", fcutoff1)
fc2 = (1 / (1 - bulkratio)) * (2 ** fcutoff) - (bulkratio / (1 - bulkratio))
fcutoff2 = - R * T * np.log(fc2)
print("FCUTOFF:", fcutoff2)


figure = plt.figure(figsize=(6.5, 4.5))
axes = figure.subplots(2, 3, sharex=True, sharey=True)

axes[1,0].axis("off")

for i, c in enumerate(comps):
    fcname = "Fold_change_" + c
    pvalname = "P_value_A_" + c

    if i<3:
        axis = axes.flatten()[i]
    else:
        axis = axes.flatten()[i+1]

    total_lipids = len(df)
    b1 = df[pvalname] < pcutoff
    b2 = np.abs(df[fcname]) > fcutoff
    b3 = np.logical_and(b1, b2)
    sig_lipids = len(df[b3])
    sdf = df[b3]
    sdf = sdf.sort_values(by=pvalname)
    print(c, sig_lipids, total_lipids, sig_lipids / total_lipids * 100)

    names = sdf["Name"].to_numpy()

    # Calculate distance from 0, 0
    fc = sdf[fcname].to_numpy()
    pval = sdf[pvalname].to_numpy()
    logp = -np.log10(pval)

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

        fc2 = 2 ** fc
        # print(fc2)
        if i != 0:
            fc2 = (1 / (1 - bulkratio)) * fc2 - (bulkratio / (1 - bulkratio))
            # print(fc2)
            bneg = fc2 < 0
            if np.any(bneg):
                print("Negative Values:", fc2[bneg])
            fc2[bneg] = 1e-12
            fc2 = - R * T * np.log(fc2)
            fc2[bneg] = xlim - 0.5
        else:
            fc2 = - R * T * np.log(fc2)

        data = fc2

        print("MINMAX:", np.amax(data), np.amin(data))

        names = subdf["Name"].values[b3]

        axis.plot(data[b1], logp[b1], symbols[j], label=u, color=colors[j], markersize=markersize)
        axis.plot(data[b2], logp[b2], symbols[j], color="grey", markersize=markersize)

    # Draw dashed line at pcutoff
    axis.plot([-xlim, xlim], [-np.log10(pcutoff), -np.log10(pcutoff)], "--", color="k")
    # Draw dashed line a fcutoff
    if i!=0:
        axis.plot([fcutoff2, fcutoff2], [0, 6], "--", color="k")
        axis.plot([-fcutoff2, -fcutoff2], [0, 6], "--", color="k")
    else:
        axis.plot([fcutoff1, fcutoff1], [0, 6], "--", color="k")
        axis.plot([-fcutoff1, -fcutoff1], [0, 6], "--", color="k")

    axis.set_ylim([0, 6])
    axis.set_xlim([-xlim, xlim])

    axis.set_title(titles[i])

    # Turn off top axis
    axis.spines['right'].set_visible(False)
    axis.spines['top'].set_visible(False)

    # Set y axis to cross at x=0
    # axis.spines['left'].set_position(('data', 0))
    if i == 0 or i == 3:
        axis.set_ylabel("–Log$_{10}$ P-Value")
    if i == 0 or i > 2:
        axis.set_xlabel("ΔG (" + units + ")")

    if i==0:
        axis.set_xticks([-4, -2, 0, 2, 4])
        #axis.set_xticklabels([-4, -2, 0, 2, 4])
        pass

    if i == 0:
        axis.legend()
plt.tight_layout()
plt.savefig("Ecoli_DG.pdf")
plt.show()
