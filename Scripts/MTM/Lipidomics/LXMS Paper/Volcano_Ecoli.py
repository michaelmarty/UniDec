import matplotlib.pyplot as plt

# Set default plot font size to 11
plt.rcParams.update({'font.size': 11})
# Set font to ariel
plt.rcParams["font.family"] = "Arial"

path = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\Papers\\LX-MS\\230517_LX_Supplemental_Data.xlsx"

df = pd.read_excel(path, sheet_name="Table S4 - Ecoli Lipid Results", skiprows=1)

print(df.columns)

classes = df["Class"].values
uclasses = np.unique(classes)
print(uclasses)
colors = get_colors(uclasses)
print(colors)

comps = ["FT_vs_Elu", "AmtB_vs_Ctrl", "AqpZ_vs_Ctrl", "FT_vs_Elu.1", "FT_vs_Elu.2"]
titles = ["Control Exchange", "AmtB Assembly", "AqpZ Assembly", "AmtB LX", "AqpZ LX"]

figure = plt.figure(figsize=(6.25, 5))
axes = figure.subplots(2, 3)

pcutoff = 0.05
fcutoff = 0.5
markersize = 3

for i, c in enumerate(comps):
    fcname = "Fold_change_" + c
    pvalname = "P_value_A_" + c

    axis = axes.flatten()[comps.index(c)]

    total_lipids = len(df)
    b1 = df[pvalname] < pcutoff
    b2 = np.abs(df[fcname]) > fcutoff
    b3 = np.logical_and(b1, b2)
    sig_lipids = len(df[b3])
    sdf = df[b3]

    print(c, sig_lipids, total_lipids, sig_lipids / total_lipids * 100)


    # sort sdf by pval
    sdf = sdf.sort_values(by=pvalname)
    names = sdf["Name"].to_numpy()
    classes = sdf["Class"].to_numpy()

    # Add text markers for the names
    posindex = 0
    negindex = 0
    for k, n in enumerate(names):
        fc = sdf[fcname].iloc[k]
        if fc > 0:
            x = 1 + 0.2
            posindex += 1
            y = 6-posindex*0.7
        else:
            x = -0.9
            negindex += 1
            y = 5.8 - negindex * 0.7
        cname = classes[k]
        color = colors[list(uclasses).index(cname)]
        axis.text(x, y, n[3:], fontsize=9, color=color, rotation=0)

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

        names = subdf["Name"].values[b3]

        axis.plot(fc[b1], logp[b1], "o", label=u, color=colors[j], markersize=markersize)
        axis.plot(fc[b2], logp[b2], "o", label=u, color="grey", markersize=markersize)



        # Draw dashed line at pcutoff
        axis.plot([-1, 1], [-np.log10(pcutoff), -np.log10(pcutoff)], "--", color="k")
        # Draw dashed line a fcutoff
        axis.plot([fcutoff, fcutoff], [0, 6], "--", color="k")
        axis.plot([-fcutoff, -fcutoff], [0, 6], "--", color="k")
        #axis.plot([0, 0], [0, 6], "-", color="k")

        axis.set_ylim([0, 6])
        axis.set_xlim([-1, 1])


        axis.set_title(titles[i])

        # Turn off top axis
        axis.spines['right'].set_visible(False)
        axis.spines['top'].set_visible(False)

        # Set y axis to cross at x=0
        #axis.spines['left'].set_position(('data', 0))
        axis.set_xlabel("Log2 Fold Change")
        if i%3 == 0:
            axis.set_ylabel("Log10 P-Value")

# Remove plot at bottom right
axes.flatten()[-1].axis("off")

plt.tight_layout()
plt.show()
