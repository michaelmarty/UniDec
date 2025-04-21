import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
import scipy.stats as stats
from matplotlib import rcParams

# Set backend to Agg
mpl.use('WxAgg')

rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['font.family'] = 'Arial'

def perform_ttest(mean1, std1, mean2, std2, n1, n2):
    tstat = (mean1 - mean2) / np.sqrt((std1**2/n1) + (std2**2/n2))

    df = n1 + n2 - 2

    p = 2 * (1 - stats.t.cdf(abs(tstat), df))

    return tstat, p



algorithms = ["ID-4", "ID-8", "FLD", "TFD"]

colors = ["skyblue", "darkblue", "yellow", "red"]
show_plots = True
save_plots = True

xlabel_size = 24
ylabel_size = 28
ytick_size = 18

if False:
    times = [33.51, 34.73, 29.28, 161.05]

    plt.bar(algorithms, times, color=colors, capsize=3, edgecolor='black')
    plt.ylabel("Time (min)", fontsize=ylabel_size)
    plt.ylim(0,200)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.yticks(np.arange(0, 201, 50), fontsize=ytick_size)
    plt.xticks(fontsize=xlabel_size)
    plt.tight_layout()
    if save_plots:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\tdp_decontimes.pdf",
                    transparent=True, dpi=300)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\tdp_decontimes.png",
                    transparent=True, dpi=300)
    if show_plots:
        plt.show()

if True:
    algorithms = ["ID-8", "FLD", "TFD"]
    colors = ["darkblue", "yellow", "red"]

    avg_unamb_prsms = [645.33,783.17,545.5]
    unamb_prsms_std = [48.10,37.81,38.82]
    avg_amb_prsms = [34.3,54.0,33.17]
    amb_prsms_std = [3.04,9.81,6.57]

    if True:
        id8_fld_tstat, id8_fld_p = perform_ttest(avg_unamb_prsms[0], unamb_prsms_std[0], avg_unamb_prsms[1], unamb_prsms_std[1], 6,6)
        print("ID-8 vs FLD Unambiguous PrSMs: t-statistic:", id8_fld_tstat, "p-value:", id8_fld_p)
        id8_tfd_tstat, id8_tfd_p = perform_ttest(avg_unamb_prsms[0], unamb_prsms_std[0], avg_unamb_prsms[2], unamb_prsms_std[2], 6,6)
        print("ID-8 vs TFD Unambiguous PrSMs: t-statistic:", id8_tfd_tstat, "p-value:", id8_tfd_p)
        fld_tfd_tstat, fld_tfd_p = perform_ttest(avg_unamb_prsms[1], unamb_prsms_std[1], avg_unamb_prsms[2], unamb_prsms_std[2], 6,6)
        print("FLD vs TFD Unambiguous PrSMs: t-statistic:", fld_tfd_tstat, "p-value:", fld_tfd_p)

        id8_fld_tstat, id8_fld_p = perform_ttest(avg_amb_prsms[0], amb_prsms_std[0], avg_amb_prsms[1], amb_prsms_std[1], 6,6)
        print("ID-8 vs FLD Ambiguous PrSMs: t-statistic:", id8_fld_tstat, "p-value:", id8_fld_p)
        id8_tfd_tstat, id8_tfd_p = perform_ttest(avg_amb_prsms[0], amb_prsms_std[0], avg_amb_prsms[2], amb_prsms_std[2], 6,6)
        print("ID-8 vs TFD Ambiguous PrSMs: t-statistic:", id8_tfd_tstat, "p-value:", id8_tfd_p)
        fld_tfd_tstat, fld_tfd_p = perform_ttest(avg_amb_prsms[1], amb_prsms_std[1], avg_amb_prsms[2], amb_prsms_std[2], 6,6)
        print("FLD vs TFD Ambiguous PrSMs: t-statistic:", fld_tfd_tstat, "p-value:", fld_tfd_p)

    legend_handles = [
        Patch(facecolor='white', edgecolor='black', label='Unambiguous'),
        Patch(facecolor='white', edgecolor='black', label='Ambiguous', hatch='\\')]

    fig, ax = plt.subplots(clip_on=False)

    bar_width = 0.40
    bar_spacing = 0.05

    index = np.arange(len(algorithms))

    bar1 = ax.bar(index - bar_width/2 - bar_spacing/2, avg_unamb_prsms, bar_width, yerr=unamb_prsms_std, capsize=3, color=colors, edgecolor='black', label="Unambiguous")
    bar2 = ax.bar(index + bar_width/2 + bar_spacing/2, avg_amb_prsms, bar_width, yerr=amb_prsms_std, capsize=3, color=colors, hatch='\\', edgecolor='black', label="Ambiguous")

    ##region Significance Brackets
    #Add brackets denoting statistically significant differences at alpha = 0.01
    plt.plot([1.25,1.25,2.225,2.225],
             [avg_unamb_prsms[1]+ 60, avg_unamb_prsms[1] + 75, avg_unamb_prsms[1] + 75, avg_unamb_prsms[1] + 60],
             linewidth=2, color="black")
    plt.text(1.725, avg_unamb_prsms[1] + 50, "*", fontsize=18)

    plt.plot([0.225, 0.225, 1.20, 1.20],
             [avg_unamb_prsms[1] + 60, avg_unamb_prsms[1] + 75, avg_unamb_prsms[1] + 75, avg_unamb_prsms[1] + 60],
             linewidth=2, color="black")
    plt.text(0.725, avg_unamb_prsms[1] + 50, "*", fontsize=18)


    plt.plot([-0.225, -0.225, 0.75, 0.75],
             [avg_unamb_prsms[1]+ 140, avg_unamb_prsms[1] + 155, avg_unamb_prsms[1] + 155, avg_unamb_prsms[1] + 140],
             linewidth=2, color="black")
    plt.text(0.225, avg_unamb_prsms[1] + 130, "*", fontsize=18)

    plt.plot([-0.225, -0.225, 1.775, 1.775],
             [avg_unamb_prsms[1] + 190, avg_unamb_prsms[1] + 205, avg_unamb_prsms[1] + 205, avg_unamb_prsms[1] + 190],
             linewidth=2, color="black")
    plt.text(0.775, avg_unamb_prsms[1] + 180, "*", fontsize=18)

    plt.plot([0.8, 0.8, 1.775, 1.775],
             [avg_unamb_prsms[1] + 140, avg_unamb_prsms[1] + 155, avg_unamb_prsms[1] + 155, avg_unamb_prsms[1] + 140],
             linewidth=2, color="black")
    plt.text(1.225, avg_unamb_prsms[1] + 130, "*", fontsize=18)
    ##endregion

    plt.ylabel("PrSM Count", fontsize=ylabel_size)
    plt.legend(handles=legend_handles, fontsize=ytick_size, loc=[0.15,1.15], handlelength=4, handleheight=1)
    plt.ylim(0, 1000)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.yticks(np.arange(0, 1001, 200), fontsize=ytick_size)
    plt.xticks(fontsize=xlabel_size)
    plt.xticks(index, algorithms)
    plt.tight_layout()
    if save_plots:
        plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\Papers\\IsoDec\\Figures\\tdp_prsms_count_new.pdf",
                    dpi=300, transparent=True)
        plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\tdp_prsms_count_new.png",
                    dpi=300, transparent=True)
    if show_plots:
        plt.show()

if True:
    times = [33.51, 34.73, 29.28, 161.05]
    total_features = [3753364, 3745614, 2297734, 5764411]
    time_per_cluster = ((np.array(times)/60) / np.array(total_features)) * 1000000000

    plt.bar(algorithms, time_per_cluster, color=colors, capsize=3, edgecolor='black', clip_on=False)
    plt.yticks(np.arange(0, 501, 100), fontsize=ytick_size)
    plt.ylabel("Time per Cluster (ns)", fontsize=ylabel_size)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xticks(fontsize=xlabel_size)
    plt.tight_layout()
    if save_plots:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\tdp_time_per_cluster.pdf",
                    transparent=True, dpi=300)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\tdp_time_per_cluster.png",
                    transparent=True, dpi=300)
    if show_plots:
        plt.show()