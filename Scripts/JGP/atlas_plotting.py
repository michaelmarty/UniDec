import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set backend to Agg
mpl.use('WxAgg')

level1_id_counts_tdp = [206,253,219,194,191,208,200,221,228]
level1_id_counts_tdp_sd = np.std(level1_id_counts_tdp)
level2A_id_counts_tdp = [39,66,44,34,51,45,41,59,40]
level2A_id_counts_tdp_sd = np.std(level2A_id_counts_tdp)
level3_id_counts_tdp = [75,111,92,82,90,107,76,94,91]
level3_id_counts_tdp_sd = np.std(level3_id_counts_tdp)

level1_id_counts_imp_lg = [105,114,117,87,112,104,95,124,125]
level1_id_counts_imp_lg_sd = np.std(level1_id_counts_imp_lg)
level2A_id_counts_imp_lg = [185,219,206,217,210,194,196,231,207]
level2A_id_counts_imp_lg_sd = np.std(level2A_id_counts_imp_lg)
level3_id_counts_imp_lg = [150,175,169,181,187,157,164,188,146]
level3_id_counts_imp_lg_sd = np.std(level3_id_counts_imp_lg)

level1_id_counts_imp_sg = [95,101,101,85,92,91,76,99,111]
level1_id_counts_imp_sg_sd = np.std(level1_id_counts_imp_sg)
level2A_id_counts_imp_sg = [161,186,185,157,188,163,145,143,166]
level2A_id_counts_imp_sg_sd = np.std(level2A_id_counts_imp_sg)
level3_id_counts_imp_sg = [103,118,111,103,117,123,150,53,130]
level3_id_counts_imp_sg_sd = np.std(level3_id_counts_imp_sg)

imp_sg_values = np.array([np.mean(level1_id_counts_imp_sg), np.mean(level2A_id_counts_imp_sg), np.mean(level3_id_counts_imp_sg)])
print("IMP-SG", imp_sg_values)
print("IMP-SG Total:", np.sum(imp_sg_values))
imp_sg_stds = np.array([level1_id_counts_imp_sg_sd, level2A_id_counts_imp_sg_sd, level3_id_counts_imp_sg_sd])
imp_lg_values = np.array([np.mean(level1_id_counts_imp_lg), np.mean(level2A_id_counts_imp_lg), np.mean(level3_id_counts_imp_lg)])
print("IMP-LG", imp_lg_values)
print("IMP-LG Total:", np.sum(imp_lg_values))
imp_lg_stds = np.array([level1_id_counts_imp_lg_sd, level2A_id_counts_imp_lg_sd, level3_id_counts_imp_lg_sd])
tdp_values = np.array([np.mean(level1_id_counts_tdp), np.mean(level2A_id_counts_tdp), np.mean(level3_id_counts_tdp)])
print("TDP", tdp_values)
print("TDP Total:", np.sum(tdp_values))
tdp_stds = np.array([level1_id_counts_tdp_sd, level2A_id_counts_tdp_sd, level3_id_counts_tdp_sd])

values = np.array([imp_lg_values, imp_sg_values, tdp_values])
stds = np.array([imp_lg_stds, imp_sg_stds, tdp_stds])

categories = ["IMP-LG", "IMP-SG", "TDP"]
if True:
    fig, ax = plt.subplots()
    bottom = np.zeros(len(categories))
    for i in range(values.shape[1]):
        ax.bar(categories, values[:, i], yerr=stds[:, i], bottom=bottom, capsize=5, label="Ambiguity Level {}".format(i+1))
        bottom += values[:, i]

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    ax.plot([0.1,0.1,0.9,0.9], [np.sum(imp_lg_values), np.sum(imp_lg_values) + 20, np.sum(imp_lg_values) + 20, np.sum(imp_sg_values)], color="black", lw=2)
    ax.text(0.5, np.sum(imp_lg_values) + 25, "*", ha="center", va="center", fontsize=16)

    ax.plot([0.15,0.15,1.9,1.9],
            [np.sum(imp_lg_values), np.sum(imp_lg_values) + 50, np.sum(imp_lg_values) + 50, np.sum(tdp_values)],
            color="black", lw=2)
    ax.text(1, np.sum(imp_lg_values) + 55, "*", ha="center", va="center", fontsize=16)

    plt.yticks(np.arange(0, 101, 25), fontsize=14)
    plt.xticks(fontsize=20)
    plt.ylabel("Number of IDs", fontsize=20)
    plt.legend(loc="upper right", bbox_to_anchor=(1,1.25), fontsize=12)
    ax.set_yticks(np.arange(0,601,100))
    plt.tight_layout()
    plt.savefig("D:\\Johnny\\SmithLab\\Atlas paper\\Figures\\id_counts_stacked.pdf", dpi=300)
    plt.savefig("D:\\Johnny\\SmithLab\\Atlas paper\\Figures\\id_counts_stacked.png", dpi=300)
    plt.show()

if True:
    tdp_totals = np.zeros(9)
    for i in range(9):
        tdp_totals[i] = np.sum([level1_id_counts_tdp[i], level2A_id_counts_tdp[i], level3_id_counts_tdp[i]])
    tdp_totals_sd = np.std(tdp_totals)

    imp_lg_totals = np.zeros(9)
    for i in range(9):
        imp_lg_totals[i] = np.sum([level1_id_counts_imp_lg[i], level2A_id_counts_imp_lg[i], level3_id_counts_imp_lg[i]])
    imp_lg_totals_sd = np.std(imp_lg_totals)

    imp_sg_totals = np.zeros(9)
    for i in range(9):
        imp_sg_totals[i] = np.sum([level1_id_counts_imp_sg[i], level2A_id_counts_imp_sg[i], level3_id_counts_imp_sg[i]])
    imp_sg_totals_sd = np.std(imp_sg_totals)

    tdp_efficiencies = tdp_totals/200
    tdp_efficiencies_sd = np.std(tdp_efficiencies)
    imp_lg_efficiencies = imp_lg_totals/200
    imp_lg_efficiencies_sd = np.std(imp_lg_efficiencies)
    imp_sg_efficiencies = imp_sg_totals/120
    imp_sg_efficiencies_sd = np.std(imp_sg_efficiencies)

    efficiencies = np.array([np.mean(imp_lg_efficiencies), np.mean(imp_sg_efficiencies), np.mean(tdp_efficiencies)])
    stds = np.array([imp_lg_efficiencies_sd, imp_sg_efficiencies_sd, tdp_efficiencies_sd])

    x = np.arange(len(categories))
    fig, ax = plt.subplots()

    dots = ax.errorbar(x, efficiencies, yerr=stds, fmt='o', capsize=5, color="black", markersize=8)

    ax.plot([0,0,0.975,0.975],
            [efficiencies[1]+stds[1]+0.1,efficiencies[1]+stds[1]+0.3,efficiencies[1]+stds[1]+0.3, efficiencies[1]+stds[1]+0.1],
            color="black", lw=2)

    ax.text(0.5, efficiencies[1]+stds[1]+0.35, "*", ha="center", va="center", fontsize=20, fontweight="bold")

    ax.plot([0,0,2,2],
            [efficiencies[1] + stds[1] + 0.40, efficiencies[1] + stds[1] + 0.60, efficiencies[1] + stds[1] + 0.60, efficiencies[1] + stds[1] + 0.40],
            color="black", lw=2)

    ax.text(1, efficiencies[1]+stds[1]+0.65, "*", ha="center", va="center", fontsize=20, fontweight="bold")

    ax.plot([1.025,1.025,2,2],
            [efficiencies[1] + stds[1] + 0.1, efficiencies[1] + stds[1] + 0.3, efficiencies[1] + stds[1] + 0.3, efficiencies[1] + stds[1] + 0.1],
            color="black", lw=2)
    ax.text(1.5, efficiencies[1]+stds[1]+0.35, "*", ha="center", va="center", fontsize=20, fontweight="bold")

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    ax.set_ylabel("Proteoform Identification\nRate (IDs/minute)", fontsize=20)
    ax.set_yticks(np.arange(0, 4.1, 1))
    ax.set_yticklabels(np.arange(0, 4.1, 1), fontsize=14)
    ax.set_xticks(np.arange(len(categories)))
    ax.set_xticklabels(categories, fontsize=20)
    ax.set_xlim(-0.25, 2.25)
    plt.tight_layout()
    plt.savefig("D:\\Johnny\\SmithLab\\Atlas paper\\Figures\\id_efficiencies.pdf", dpi=300)
    plt.savefig("D:\\Johnny\\SmithLab\\Atlas paper\\Figures\\id_efficiencies.png", dpi=300)
    plt.show()



