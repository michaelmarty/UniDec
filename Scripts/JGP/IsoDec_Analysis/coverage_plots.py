import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as stats
from matplotlib import rcParams


# Set backend to Agg
mpl.use('WxAgg')

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

colors = ["skyblue", "darkblue", "yellow", "red", "green"]
algorithms = ["ID-4", "ID-8", "FLD", "TFD", "eTH"]

uvpd_coverages = [71.04,69.88,54.05,69.50, 72.59]
uvpd_fragsmatched = [59.91,59.29,60.62,34.84, 29.39]

etd_coverages = [76.06,76.45,55.6,76.06, 76.06]
etd_fragsmatched = [63.47,64.60,48.32,24.62, 24.60]

ecd_coverages = [70.59,68.63,40.51,66.67]
ecd_fragsmatched = [43.75,42.20,41.82,40.81]

show_plots = True
save_plots = False

xlabel_size = 28
ylabel_size = 28
ytick_size = 18

#Make a 3 category paired barplot of the coverages on each dataset
if True:
    frag_types = ["UVPD", "ETD", "ECD"]
    x = np.arange(len(frag_types))
    width = 0.15
    spacing = 0.025
    fig, ax = plt.subplots(figsize=(8,5))

    id8_coverages = [uvpd_coverages[0], etd_coverages[0], ecd_coverages[0]]
    id4_coverages = [uvpd_coverages[1], etd_coverages[1], ecd_coverages[1]]
    fld_coverages = [uvpd_coverages[2], etd_coverages[2], ecd_coverages[2]]
    tfd_coverages = [uvpd_coverages[3], etd_coverages[3], ecd_coverages[3]]
    ethrash_coverages = [uvpd_coverages[4], etd_coverages[4], 0]

    bars1 = ax.bar(x - 2 * width - (2 * spacing), id4_coverages, width, label="IsoDec-4bin", color=colors[0], edgecolor='black')
    bars2 = ax.bar(x - 1 * width - (1 * spacing), id8_coverages, width, label="IsoDec-8bin", color=colors[1], edgecolor='black')
    bars3 = ax.bar(x, fld_coverages, width, label="FLASHDeconv", color=colors[2], edgecolor='black')
    bars4 = ax.bar(x + 1 * width + (1 * spacing), tfd_coverages, width, label="TopFD", color=colors[3], edgecolor='black')
    bars5 = ax.bar(x + 2 * width + (2 * spacing), ethrash_coverages, width, label="eTHRASH", color=colors[4], edgecolor='black')

    plt.ylim(0,100)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    ax.set_xticks(x)
    ax.set_xticklabels(frag_types, fontsize=xlabel_size)
    ax.set_ylabel("Sequence Coverage (%)", fontsize=ylabel_size)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_yticklabels(np.arange(0, 101, 20), fontsize=ytick_size)
    plt.tight_layout()
    plt.savefig("C:\\users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\coverage_pairedbar.png",
                dpi=300, transparent=True)
    plt.savefig("C:\\users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\coverage_pairedbar.pdf",
                dpi=300, transparent=True)
    plt.show()

if True:
    id8_frags = [uvpd_fragsmatched[0], etd_fragsmatched[0], ecd_fragsmatched[0]]
    id4_frags = [uvpd_fragsmatched[1], etd_fragsmatched[1], ecd_fragsmatched[1]]
    fld_frags = [uvpd_fragsmatched[2], etd_fragsmatched[2], ecd_fragsmatched[2]]
    tfd_frags = [uvpd_fragsmatched[3], etd_fragsmatched[3], ecd_fragsmatched[3]]
    ethrash_frags = [uvpd_fragsmatched[4], etd_fragsmatched[4], 0]

    frag_types = ["UVPD", "ETD", "ECD"]
    x = np.arange(len(frag_types))
    width = 0.15
    spacing = 0.025

    fig, ax = plt.subplots(figsize=(8,5))
    bars1 = ax.bar(x-2*width - (2*spacing), id4_frags, width, label="IsoDec-4bin", color=colors[0], edgecolor='black')
    bars2 = ax.bar(x-1*width - (1*spacing), id8_frags, width, label="IsoDec-8bin", color=colors[1], edgecolor='black')
    bars3 = ax.bar(x, fld_frags, width, label="FLASHDeconv", color=colors[2], edgecolor='black')
    bars4 = ax.bar(x+1*width + (1*spacing), tfd_frags, width, label="TopFD", color=colors[3], edgecolor='black')
    bars5 = ax.bar(x+2*width + (2*spacing), ethrash_frags, width, label="eTHRASH", color=colors[4], edgecolor='black')

    plt.ylim(0,100)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    ax.set_xticks(x)
    ax.set_xticklabels(frag_types, fontsize=xlabel_size)
    ax.set_ylabel("Fragments Matched (%)", fontsize=ylabel_size)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_yticklabels(np.arange(0, 101, 20), fontsize=ytick_size)
    plt.legend(loc='upper right', fontsize=18)
    plt.tight_layout()
    plt.savefig("C:\\users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\fragsmatched_pairedbar.png",
                dpi=300, transparent=True)
    plt.savefig("C:\\users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\fragsmatched_pairedbar.pdf",
                dpi=300, transparent=True)
    plt.show()

if False:
    prsm_coverages = [63.00,62.70,49.38,51.06]
    prsm_coverages_sd =  [20.23,20.44,21.82,17.23]

    #Coverage significant differences:
    #ID-4 vs FLASHDeconv
    #ID-4 vs TopFD
    #ID-8 vs TopFD
    #ID-8 vs FLASHDeconv

    prsm_fragsmatched = [42.50,42.00,45.21,35.85]
    prsm_fragsmatched_sd = [11.02,11.17,14.62,14.18]

    #Plot the prsm coverages
    if False:
        plt.bar(algorithms, prsm_coverages, yerr=prsm_coverages_sd, color=colors, capsize=5, edgecolor='black')

        plt.plot([0,0,2,2],
                 [prsm_coverages[0] + prsm_coverages_sd[0] + 9, prsm_coverages[0] + prsm_coverages_sd[0] + 11,
                  prsm_coverages[0] + prsm_coverages_sd[0] + 11, prsm_coverages[0]+ prsm_coverages_sd[0] + 9],
                 color='black', linestyle='-', linewidth=2)
        plt.text(1, prsm_coverages[0] + prsm_coverages_sd[0] + 10, "*", fontsize=14)

        plt.plot([0, 0, 3, 3],
                 [prsm_coverages[0] + prsm_coverages_sd[0] + 13, prsm_coverages[0] + prsm_coverages_sd[0] + 15,
                  prsm_coverages[0] + prsm_coverages_sd[0] + 15, prsm_coverages[0] + prsm_coverages_sd[0] + 13],
                 color='black', linestyle='-', linewidth=2)
        plt.text(1.5, prsm_coverages[0] + prsm_coverages_sd[0] + 14, "*", fontsize=14)

        plt.plot([1, 1, 2, 2],
                 [prsm_coverages[0] + prsm_coverages_sd[0] + 1, prsm_coverages[0] + prsm_coverages_sd[0] + 3,
                  prsm_coverages[0] + prsm_coverages_sd[0] + 3, prsm_coverages[0] + prsm_coverages_sd[0] + 1],
                 color='black', linestyle='-', linewidth=2)
        plt.text(1.5, prsm_coverages[0] + prsm_coverages_sd[0] + 2, "*", fontsize=14)

        plt.plot([1, 1, 3, 3],
                 [prsm_coverages[0] + prsm_coverages_sd[0] + 5, prsm_coverages[0] + prsm_coverages_sd[0] + 7,
                  prsm_coverages[0] + prsm_coverages_sd[0] + 7, prsm_coverages[0] + prsm_coverages_sd[0] + 5],
                 color='black', linestyle='-', linewidth=2)
        plt.text(2.25, prsm_coverages[0] + prsm_coverages_sd[0] + 6, "*", fontsize=14)

        plt.ylabel("Sequence Coverage (%)", fontsize=ylabel_size)
        plt.ylim(0,100)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.yticks(np.arange(0, 101, 20), fontsize=ytick_size)
        plt.xticks(fontsize=xlabel_size)
        plt.tight_layout()
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\prsm_coverage.png", dpi=300)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\prsm_coverage.pdf", dpi=300)
        plt.show()

    if False:
        #Plot the prsm frags matched
        #Fragsmatched significant differences:
        #ID4 vs FLD
        #ID4 vs TFD
        #ID8 vs TFD
        #ID8 vs FLD
        #FLD vs TFD

        plt.bar(algorithms, prsm_fragsmatched, yerr=prsm_fragsmatched_sd, color=colors, capsize=5, edgecolor='black')

        plt.plot([2.05,2.05,3,3],
                 [prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 2, prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 4,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 4, prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 2],
                 color='black', linestyle='-', linewidth=2)
        plt.text(2.5, prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 3, "*", fontsize=14)

        plt.plot([1,1,1.95,1.95],
                 [prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 2,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 4,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 4,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 2],
                 color='black', linestyle='-', linewidth=2)
        plt.text(1.5, prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 3, "*", fontsize=14)

        plt.plot([1, 1, 3,3],
                 [prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 6,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 8,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 8,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 6],
                 color='black', linestyle='-', linewidth=2)
        plt.text(2, prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 7, "*", fontsize=14)


        plt.plot([0,0,1.95,1.95],
                 [prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 10,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 12,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 12,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 10],
                 color='black', linestyle='-', linewidth=2)
        plt.text(1, prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 11, "*", fontsize=14)

        plt.plot([0, 0, 3,3],
                 [prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 14,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 16,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 16,
                  prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 14],
                 color='black', linestyle='-', linewidth=2)
        plt.text(1, prsm_fragsmatched[2] + prsm_fragsmatched_sd[2] + 15, "*", fontsize=14)

        plt.ylabel("Fragments Matched (%)", fontsize=ylabel_size)
        plt.ylim(0,100)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.yticks(np.arange(0, 101, 20), fontsize=ytick_size)
        plt.xticks(fontsize=xlabel_size)
        plt.tight_layout()
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\prsm_fragsmatched.png", dpi=300)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\prsm_fragsmatched.pdf", dpi=300)
        plt.show()

if False:
    algorithms = ["ID-4", "ID-8", "ID-THRASH"]
    colors = ["skyblue", "darkblue", "orange"]

    uvpd_coverages = [71.04, 69.88, 60.62]
    uvpd_fragsmatched = [59.91, 59.29, 57.14]

    etd_coverages = [76.06, 76.45, 74.9]
    etd_fragsmatched = [63.47, 64.60, 63.03]

    ecd_coverages = [70.59, 68.63, 72.55]
    ecd_fragsmatched = [43.75, 42.20, 37.0]


    frag_types = ["UVPD", "ETD", "ECD"]
    x = np.arange(len(frag_types))
    width = 0.3
    spacing = 0.025


    id8_coverages = [uvpd_coverages[0], etd_coverages[0], ecd_coverages[0]]
    id8_frags = [uvpd_fragsmatched[0], etd_fragsmatched[0], ecd_fragsmatched[0]]
    id4_coverages = [uvpd_coverages[1], etd_coverages[1], ecd_coverages[1]]
    id4_frags = [uvpd_fragsmatched[1], etd_fragsmatched[1], ecd_fragsmatched[1]]
    idthrash_coverages = [uvpd_coverages[2], etd_coverages[2], ecd_coverages[2]]
    idthrash_frags = [uvpd_fragsmatched[2], etd_fragsmatched[2], ecd_fragsmatched[2]]


    if True:
        fig, ax = plt.subplots(figsize=(8, 5))
        bars1 = ax.bar(x-width - spacing, id4_coverages, width, label="IsoDec-4bin", color=colors[0], edgecolor='black')
        bars2 = ax.bar(x, id8_coverages, width, label="IsoDec-8bin", color=colors[1], edgecolor='black')
        bars3 = ax.bar(x+width + spacing, idthrash_coverages, width, label="IsoDec-THRASH", color=colors[2], edgecolor='black')

        plt.ylim(0,100)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        ax.set_xticks(x)
        ax.set_xticklabels(frag_types, fontsize=xlabel_size)
        ax.set_ylabel("Sequence Coverage (%)", fontsize=ylabel_size)
        ax.set_yticks(np.arange(0, 101, 20))
        ax.set_yticklabels(np.arange(0, 101, 20), fontsize=ytick_size)
        plt.tight_layout()
        if save_plots:
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\coverage_pairedbar_thrash.png", dpi=300)
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\coverage_pairedbar_thrash.pdf", dpi=300)
        plt.show()

    if True:
        fig, ax = plt.subplots(figsize=(8, 5))
        bars1 = ax.bar(x - width - spacing, id4_frags, width, label="IsoDec-4bin", color=colors[0],
                       edgecolor='black')
        bars2 = ax.bar(x, id8_frags, width, label="IsoDec-8bin", color=colors[1], edgecolor='black')
        bars3 = ax.bar(x + width + spacing, idthrash_frags, width, label="IsoDec-THRASH", color=colors[2],
                       edgecolor='black')

        plt.legend(fontsize=18)
        plt.ylim(0, 100)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        ax.set_xticks(x)
        ax.set_xticklabels(frag_types, fontsize=xlabel_size)
        ax.set_ylabel("Fragments Matched (%)", fontsize=ylabel_size)
        ax.set_yticks(np.arange(0, 101, 20))
        ax.set_yticklabels(np.arange(0, 101, 20), fontsize=ytick_size)
        plt.tight_layout()
        if save_plots:
            plt.savefig(
                "D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\fragsmatched_pairedbar_thrash.png",
                dpi=300)
            plt.savefig(
                "D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\fragsmatched_pairedbar_thrash.pdf",
                dpi=300)
        plt.show()

if False:
    z = stats.norm.ppf(0.975)

    categories = ["ID-FLD", "ID-TFD"]
    colors = ["green", "purple"]
    delta_coverages = [5.67,12.38]
    delta_coverages_std = [6.85,9.75]
    coverage_margins = z * (delta_coverages_std/np.sqrt(1754))

    delta_fragsmatched = [0.25,8.97]
    delta_fragsmatched_std = [6.71,10.10]
    fragsmatched_margins = z * (delta_fragsmatched_std/np.sqrt(1754))


    if True:
        plt.bar(categories, delta_coverages, yerr=coverage_margins, color=colors, capsize=5, edgecolor='black')
        plt.plot([-1, 2], [0, 0], color='black', linewidth=2)
        plt.xlim(-0.5,1.5)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xticks(fontsize=xlabel_size)
        plt.ylabel("Delta Sequence\nCoverage (%)", fontsize=ylabel_size)
        plt.yticks(np.arange(-20, 21, 5), fontsize=ytick_size)
        plt.tight_layout()
        plt.ylim(0, 15)
        plt.ylabel("Delta Sequence\nCoverage (%)", fontsize=ylabel_size)
        plt.tight_layout()
        if save_plots:
            plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\delta_coverage_bar.png",
                        transparent=True, dpi=300)
            plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\delta_coverage_bar.pdf",
                        transparent=True, dpi=300)
        plt.show()

    if True:
        plt.bar(categories, delta_fragsmatched, yerr=fragsmatched_margins, color=colors, capsize=5, edgecolor='black')
        plt.plot([-1, 2], [0, 0], color='black', linewidth=2)
        plt.xlim(-0.5,1.5)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xticks(fontsize=xlabel_size)
        plt.ylabel("Delta Fragments\nMatched (%)", fontsize=ylabel_size)
        plt.yticks(np.arange(-20, 21, 5), fontsize=ytick_size)
        plt.tight_layout()
        plt.ylim(-5, 10)
        plt.tight_layout()
        if save_plots:
            plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\delta_fragsmatched_bar.png",
                        transparent=True, dpi=300)
            plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\delta_fragsmatched_bar.pdf",
                        transparent=True, dpi=300)
        plt.show()

if False:
    algs = ["ID-8", "FLD", "TFD"]
    colors = ["darkblue", "yellow", "red"]

    id8_pvalid = [97.22,94.44,92.96]
    id8_numvalid = [210,204,132]

    fld_pvalid = [96.09,88.49,80.41]
    fld_numvalid = [123,123,78]

    tfd_pvalid = [97.14,94.20,87.43]
    tfd_numvalid = [102,130,167]

    if True:
        z = 9.925

        id8_fld_pvalid_diffs = np.array(id8_pvalid) - np.array(fld_pvalid)
        id8_fld_pvalid_diffs_mean = np.mean(id8_fld_pvalid_diffs)
        id8_fld_pvalid_diffs_std = np.std(id8_fld_pvalid_diffs)
        id8_fld_pvalid_diffs_margin = z * (id8_fld_pvalid_diffs_std/np.sqrt(len(id8_fld_pvalid_diffs)))

        id8_tfd_pvalid_diffs = np.array(id8_pvalid) - np.array(tfd_pvalid)
        id8_tfd_pvalid_diffs_mean = np.mean(id8_tfd_pvalid_diffs)
        id8_tfd_pvalid_diffs_std = np.std(id8_tfd_pvalid_diffs)
        id8_tfd_pvalid_diffs_margin = z * (id8_tfd_pvalid_diffs_std/np.sqrt(len(id8_tfd_pvalid_diffs)))

        id8_fld_pvalid_t_stat, id8_fld_pvalid_pval = stats.ttest_rel(id8_pvalid, fld_pvalid)
        id8_tfd_pvalid_t_stat, id8_tfd_pvalid_pval = stats.ttest_rel(id8_pvalid, tfd_pvalid)

        print("Pvalid ID-8 vs FLD pval", id8_fld_pvalid_pval)
        print("Pvalid ID-8 vs TFD pval", id8_tfd_pvalid_pval)

        comparisons = ["ID-8 - FLD", "ID-8 - TFD"]
        plt.bar(comparisons, [id8_fld_pvalid_diffs_mean, id8_tfd_pvalid_diffs_mean],
                yerr=[id8_fld_pvalid_diffs_margin, id8_tfd_pvalid_diffs_margin], color=["green", "purple"],
                capsize=5, edgecolor='black')
        plt.axhline(0, color="k")
        plt.yticks(np.arange(-10, 21, 5), fontsize=ytick_size)
        plt.ylim(-10,20)
        plt.xticks(fontsize=xlabel_size)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.ylabel("Difference in\nValidated Clusters (%)", fontsize=ylabel_size)
        plt.tight_layout()
        if False:
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\pvalid_diffs_bar.png", transparent=true, dpi=300)
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\pvalid_diffs_bar.pdf", transparent=true, dpi=300)
        plt.show()

    if True:
        id8_fld_numvalid_diffs = np.array(id8_numvalid) - np.array(fld_numvalid)
        id8_fld_numvalid_diffs_mean = np.mean(id8_fld_numvalid_diffs)
        id8_fld_numvalid_diffs_std = np.std(id8_fld_numvalid_diffs)
        id8_fld_numvalid_diffs_margin = z * (id8_fld_numvalid_diffs_std/np.sqrt(len(id8_fld_numvalid_diffs)))

        id8_tfd_numvalid_diffs = np.array(id8_numvalid) - np.array(tfd_numvalid)
        id8_tfd_numvalid_diffs_mean = np.mean(id8_tfd_numvalid_diffs)
        id8_tfd_numvalid_diffs_std = np.std(id8_tfd_numvalid_diffs)
        id8_tfd_numvalid_diffs_margin = z * (id8_tfd_numvalid_diffs_std/np.sqrt(len(id8_tfd_numvalid_diffs)))

        #Perform a paired sample t-test on the numvalids
        id8_fld_numvalid_t_stat, id8_fld_numvalid_pval = stats.ttest_rel(id8_numvalid, fld_numvalid)
        id8_tfd_numvalid_t_stat, id8_tfd_numvalid_pval = stats.ttest_rel(id8_numvalid, tfd_numvalid)
        print("Numvalid ID-8 vs FLD pval", id8_fld_numvalid_pval)
        print("Numvalid ID-8 vs TFD pval", id8_tfd_numvalid_pval)

        comparisons = ["ID-8 - FLD", "ID-8 - TFD"]
        plt.bar(comparisons, [id8_fld_numvalid_diffs_mean, id8_tfd_numvalid_diffs_mean],
                yerr=[id8_fld_numvalid_diffs_margin, id8_tfd_numvalid_diffs_margin], color=["green", "purple"],
                capsize=5, edgecolor='black')
        plt.axhline(0, color="k")
        plt.xticks(fontsize=xlabel_size)
        plt.yticks(np.arange(-150, 251, 50), fontsize=ytick_size)
        plt.ylim(-150,250)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.ylabel("Difference in\nValidated Clusters", fontsize=ylabel_size)
        plt.tight_layout()
        if False:
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\numvalid_diffs_bar.png", dpi=300)
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\numvalid_diffs_bar.pdf", dpi=300)
        plt.show()

