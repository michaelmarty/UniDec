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

comps = ["FT_vs_Elu", "AmtB_vs_Ctrl", "AqpZ_vs_Ctrl", "FT_vs_Elu.1", "FT_vs_Elu.2"]
titles = ["Control Exchange", "AmtB Assembly", "AqpZ Assembly", "AmtB LX", "AqpZ LX"]

pcutoff = 0.05
fcutoff = 0.5
markersize = 3
xlim = 10

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

dg = True

for i, c in enumerate(comps):
    fcname = "Fold_change_" + c
    pvalname = "P_value_A_" + c

    figure = plt.figure(figsize=(3.25, 3))

    axis = figure.gca()

    total_lipids = len(df)
    b1 = df[pvalname] < pcutoff
    b2 = np.abs(df[fcname]) > fcutoff
    b3 = np.logical_and(b1, b2)
    sig_lipids = len(df[b3])
    sdf = df[b3]
    sdf = sdf.sort_values(by=pvalname)
    print(c, sig_lipids, total_lipids, sig_lipids / total_lipids * 100)

    names = sdf["Name"].to_numpy()
    for n in names:
        if "CL" in n:
            print(n)

    # Calculate distance from 0, 0
    fc = sdf[fcname].to_numpy()
    pval = sdf[pvalname].to_numpy()
    logp = -np.log10(pval)

    distance = np.sqrt(fc ** 2 + logp ** 2)

    # Set a cutoff at the 10th most abundant lipid
    dsort = np.sort(distance)[::-1]
    print(dsort)
    if len(dsort) > 10:
        dcutoff = dsort[10]
    else:
        if len(dsort) == 0:
            dcutoff = 0
        else:
            dcutoff = dsort[-1]

    print(dcutoff)

    b1 = distance > dcutoff
    sdf = sdf[b1]

    # sort sdf by pval
    sdf = sdf.sort_values(by=pvalname)
    names = sdf["Name"].to_numpy()
    classes = sdf["Class"].to_numpy()

    # Add text markers for the names
    posindex = 0
    negindex = 0
    for k, n in enumerate(names):
        fc = sdf[fcname].iloc[k]
        x = fc
        y = -np.log10(sdf[pvalname].iloc[k])
        name = n.split("|")[0]

        """
        if fc > 0:
            x = 1 + 0.2
            posindex += 1
            y = 6 - posindex * 0.7
        else:
            x = -0.9
            negindex += 1
            y = 5.8 - negindex * 0.7"""
        cname = classes[k]
        color = colors[list(uclasses).index(cname)]
        axis.text(x, y, name, fontsize=8, color=color, rotation=0)
    print(names)
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

        if dg:
            data = fc2
        else:
            data = fc

        names = subdf["Name"].values[b3]

        axis.plot(data[b1], logp[b1], "o", label=u, color=colors[j], markersize=markersize)
        axis.plot(data[b2], logp[b2], "o", color="grey", markersize=markersize)

        # Draw dashed line at pcutoff
        axis.plot([-xlim, xlim], [-np.log10(pcutoff), -np.log10(pcutoff)], "--", color="k")
        # Draw dashed line a fcutoff
        axis.plot([fcutoff, fcutoff], [0, 6], "--", color="k")
        axis.plot([-fcutoff, -fcutoff], [0, 6], "--", color="k")
        # axis.plot([0, 0], [0, 6], "-", color="k")

        axis.set_ylim([0, 6])
        axis.set_xlim([-xlim, xlim])

        axis.set_title(titles[i])

        # Turn off top axis
        axis.spines['right'].set_visible(False)
        axis.spines['top'].set_visible(False)

        # Set y axis to cross at x=0
        # axis.spines['left'].set_position(('data', 0))

        axis.set_ylabel("Log10 P-Value")
        if True:
            if dg:
                axis.set_xlabel("Delta G (" + units + ")")
            else:
                axis.set_xlabel("Log2 Fold Change")
    plt.legend()
    plt.tight_layout()

    # plt.savefig(c + ".pdf")
plt.show()
