import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set backend to Agg
mpl.use('WxAgg')

colors = ["skyblue", "darkblue", "yellow", "red"]
algorithms = ["ID-4", "ID-8", "FLD", "TFD"]

total_features = [1631.17, 1673.83, 1148.17, 1646]
total_features_std = [56.12, 65.84, 47.71, 56.59]
#significant differences (p<0.01):
#ID-4 vs FLD
#ID-8 vs FLD
#FLD vs TFD

features_matched = [76.12, 75.08, 64.17, 65.88]
features_matched_std = [6.13, 6.17, 6.11, 5.27]
#significant differences (p<0.01):
#ID-4 vs FLD
#ID-4 vs TFD
#ID-8 vs FLD
#ID-8 vs TFD

intensity_matched = [95.45, 95.54, 92.88, 95.02]
intensity_matched_std = [2.01, 2.06, 2.84, 1.50]
#significant differences (p<0.01):
#ID-4 vs FLD
#ID-8 vs FLD

ms1_clusters = [159669.67, 158966.17, 151274.83, 301496.0]
ms1_clusters_std = [6448.23, 6394.54, 2308.78, 11233.93]
#significant differences (p<0.01):
#ID-4 vs TFD
#ID-8 vs TFD
#FLD vs TFD


ms2_clusters = [465891, 465302.83, 231680.83, 659239.17]
ms2_clusters_std = [82163.30, 85566.84, 23673.39, 87871.23]
#significant differences (p<0.01):
#ID-4 vs FLD
#ID-4 vs TFD
#ID-8 vs FLD
#ID-8 vs TFD
#FLD vs TFD

pks_in_masses = [70.61,70.22,60.19,33.04]
pks_in_masses_std = [0.89,0.87,0.98,1.07]
#significant differences (p<0.01):
#ID-4 vs FLD
#ID-4 vs TFD
#ID-8 vs FLD
#ID-8 vs TFD


xlabel_size = 20
ylabel_size = 20
ytick_size = 16

show_plots = True
save_plots = False


if True:
    #Make a six-panel figure with all of the metrics!
    fig, axs = plt.subplots(3, 2, figsize=(10, 20))

    ## region Total MS1 Clusters
    axs[0,0].bar(algorithms, ms1_clusters, yerr=ms1_clusters_std, color=colors,
                 capsize=5, edgecolor='black', clip_on=False)

    axs[0,0].plot([2, 2, 3, 3],
             [ms1_clusters[3] + ms1_clusters_std[3] + 10000, ms1_clusters[3] + ms1_clusters_std[3] + 15000,
              ms1_clusters[3] + ms1_clusters_std[3] + 15000, ms1_clusters[3] + ms1_clusters_std[3] + 10000],
             color='black', linewidth=1)
    axs[0,0].text(2.5, ms1_clusters[3] + ms1_clusters_std[3] + 10000, "*", fontsize=14, ha='center')

    axs[0,0].plot([1, 1, 3, 3],
             [ms1_clusters[3] + ms1_clusters_std[3] + 30000, ms1_clusters[3] + ms1_clusters_std[3] + 35000,
              ms1_clusters[3] + ms1_clusters_std[3] + 35000, ms1_clusters[3] + ms1_clusters_std[3] + 30000],
             color='black', linewidth=1)
    axs[0,0].text(2, ms1_clusters[3] + ms1_clusters_std[3] + 30000, "*", fontsize=14, ha='center')

    axs[0,0].plot([0, 0, 3, 3],
             [ms1_clusters[3] + ms1_clusters_std[3] + 52500, ms1_clusters[3] + ms1_clusters_std[3] + 57500,
              ms1_clusters[3] + ms1_clusters_std[3] + 57500, ms1_clusters[3] + ms1_clusters_std[3] + 52500],
             color='black', linewidth=1)
    axs[0,0].text(1.5, ms1_clusters[3] + ms1_clusters_std[3] + 52500, "*", fontsize=14, ha='center')

    axs[0,0].spines['top'].set_visible(False)
    axs[0,0].spines['right'].set_visible(False)
    axs[0,0].set_ylim(0, 400000)
    axs[0,0].set_yticks(np.arange(0, 400001, 100000), labels=["0", "1e5", "2e5", "3e5", "4e5"], fontsize=ytick_size)
    axs[0,0].set_ylabel("MS1 Clusters", fontsize=ylabel_size)
    axs[0, 0].set_xticks(np.arange(0, len(algorithms)), labels=algorithms, fontsize=xlabel_size)

    ## endregion Total MS1 Clusters

    ## region Total MS2 Clusters
    axs[0,1].bar(algorithms, ms2_clusters, yerr=ms2_clusters_std, color=colors, capsize=5, edgecolor='black', clip_on=False)

    axs[0,1].plot([1, 1, 2, 2],
             [ms2_clusters[0] + ms2_clusters_std[0] + 15000, ms2_clusters[0] + ms2_clusters_std[0] + 25000,
              ms2_clusters[0] + ms2_clusters_std[0] + 25000, ms2_clusters[0] + ms2_clusters_std[0] + 15000],
             color='black', linewidth=1)
    axs[0,1].text(1.5, ms2_clusters[0] + ms2_clusters_std[0] + 15000, "*", fontsize=14, ha='center')

    axs[0,1].plot([0, 0, 2, 2],
             [ms2_clusters[0] + ms2_clusters_std[0] + 65000, ms2_clusters[0] + ms2_clusters_std[0] + 75000,
              ms2_clusters[0] + ms2_clusters_std[0] + 75000, ms2_clusters[0] + ms2_clusters_std[0] + 65000],
             color='black', linewidth=1)
    axs[0,1].text(1, ms2_clusters[0] + ms2_clusters_std[0] + 65000, "*", fontsize=14, ha='center')

    axs[0,1].plot([2, 2, 3, 3],
             [ms2_clusters[3] + ms2_clusters_std[3] + 15000, ms2_clusters[3] + ms2_clusters_std[3] + 25000,
              ms2_clusters[3] + ms2_clusters_std[3] + 25000, ms2_clusters[3] + ms2_clusters_std[3] + 15000],
             color='black', linewidth=1)
    axs[0,1].text(2.5, ms2_clusters[3] + ms2_clusters_std[3] + 15000, "*", fontsize=14, ha='center')

    axs[0,1].plot([1, 1, 3, 3],
             [ms2_clusters[3] + ms2_clusters_std[3] + 65000, ms2_clusters[3] + ms2_clusters_std[3] + 75000,
              ms2_clusters[3] + ms2_clusters_std[3] + 75000, ms2_clusters[3] + ms2_clusters_std[3] + 65000],
             color='black', linewidth=1)
    axs[0,1].text(2, ms2_clusters[3] + ms2_clusters_std[3] + 65000, "*", fontsize=14, ha='center')

    axs[0,1].plot([0, 0, 3, 3],
             [ms2_clusters[3] + ms2_clusters_std[3] + 115000, ms2_clusters[3] + ms2_clusters_std[3] + 125000,
              ms2_clusters[3] + ms2_clusters_std[3] + 125000, ms2_clusters[3] + ms2_clusters_std[3] + 115000],
             color='black', linewidth=1)
    axs[0,1].text(1.5, ms2_clusters[3] + ms2_clusters_std[3] + 115000, "*", fontsize=14, ha='center')

    axs[0,1].spines['top'].set_visible(False)
    axs[0,1].spines['right'].set_visible(False)
    axs[0,1].set_ylabel("MS2 Clusters", fontsize=ylabel_size)
    axs[0,1].set_yticks(np.arange(0, 800001, 200000), labels=["0", "2e5", "4e5", "6e5", "8e5"], fontsize=ytick_size)
    axs[0,1].set_xticks(np.arange(0, len(algorithms)), labels=algorithms, fontsize=xlabel_size)

    ## endregion Total MS2 Clusters

    ## region MS1 Features
    axs[1,0].bar(algorithms, total_features, yerr=total_features_std, color=colors, capsize=5, edgecolor='black',
            clip_on=False)

    ## region Add significance brackets and stars
    axs[1,0].plot([0, 0, 2, 2],
             [total_features[0] + total_features_std[0] + 350, total_features[0] + total_features_std[0] + 400,
              total_features[0] + total_features_std[0] + 400, total_features[0] + total_features_std[0] + 350],
             color='black', linewidth=1, clip_on=False)
    axs[1,0].text(1, total_features[0] + total_features_std[0] + 375, "*", fontsize=14, ha='center')

    axs[1,0].plot([1, 1, 1.98, 1.98],
             [total_features[0] + total_features_std[0] + 150, total_features[0] + total_features_std[0] + 200,
              total_features[0] + total_features_std[0] + 200, total_features[0] + total_features_std[0] + 150],
             color='black', linewidth=1, clip_on=False)
    axs[1,0].text(1.5, total_features[0] + total_features_std[0] + 175, "*", fontsize=14, ha='center')

    axs[1,0].plot([2.02, 2.02, 3, 3],
             [total_features[0] + total_features_std[0] + 150, total_features[0] + total_features_std[0] + 200,
              total_features[0] + total_features_std[0] + 200, total_features[0] + total_features_std[0] + 150],
             color='black', linewidth=1, clip_on=False)
    axs[1,0].text(2.5, total_features[0] + total_features_std[0] + 175, "*", fontsize=14, ha='center')
    ## endregion Add significance brackets and stars

    axs[1,0].spines['top'].set_visible(False)
    axs[1,0].spines['right'].set_visible(False)
    axs[1,0].set_ylabel("Total Features", fontsize=ylabel_size)
    axs[1,0].set_yticks(np.arange(0, 2001, 1000), labels=np.arange(0,2001,1000), fontsize=ytick_size)
    axs[1,0].set_ylim(0, 2000)
    axs[1,0].set_xticks(np.arange(0, len(algorithms)), labels=algorithms, fontsize=xlabel_size)

    ## endregion MS1 Features

    ## region MS1 Clusters in Features
    axs[1,1].bar(algorithms, pks_in_masses, yerr=pks_in_masses_std, color=colors, capsize=5, edgecolor='black',
            clip_on=False)

    ## region Add significance brackets and stars
    axs[1,1].plot([0, 0, 2, 2],
             [pks_in_masses[0] + pks_in_masses_std[0] + 3, pks_in_masses[0] + pks_in_masses_std[0] + 5,
              pks_in_masses[0] + pks_in_masses_std[0] + 5, pks_in_masses[0] + pks_in_masses_std[0] + 3],
             color='black', linewidth=1, clip_on=False)
    axs[1,1].text(1, pks_in_masses[0] + pks_in_masses_std[0] + 4, "*", fontsize=14, ha='center')

    axs[1,1].plot([0, 0, 3, 3],
             [pks_in_masses[0] + pks_in_masses_std[0] + 9, pks_in_masses[0] + pks_in_masses_std[0] + 11,
              pks_in_masses[0] + pks_in_masses_std[0] + 11, pks_in_masses[0] + pks_in_masses_std[0] + 9],
             color='black', linewidth=1, clip_on=False)
    axs[1,1].text(1.5, pks_in_masses[0] + pks_in_masses_std[0] + 10, "*", fontsize=14, ha='center')

    axs[1,1].plot([1, 1, 2, 2],
             [pks_in_masses[0] + pks_in_masses_std[0] + 15, pks_in_masses[0] + pks_in_masses_std[0] + 17,
              pks_in_masses[0] + pks_in_masses_std[0] + 17, pks_in_masses[0] + pks_in_masses_std[0] + 15],
             color='black', linewidth=1, clip_on=False)
    axs[1,1].text(1.5, pks_in_masses[0] + pks_in_masses_std[0] + 16, "*", fontsize=14, ha='center')

    axs[1,1].plot([1, 1, 3, 3],
             [pks_in_masses[0] + pks_in_masses_std[0] + 21, pks_in_masses[0] + pks_in_masses_std[0] + 23,
              pks_in_masses[0] + pks_in_masses_std[0] + 23, pks_in_masses[0] + pks_in_masses_std[0] + 21],
             color='black', linewidth=1, clip_on=False)
    axs[1,1].text(2, pks_in_masses[0] + pks_in_masses_std[0] + 22, "*", fontsize=14, ha='center')

    ## endregion Add significance brackets and stars

    axs[1,1].spines['top'].set_visible(False)
    axs[1,1].spines['right'].set_visible(False)
    axs[1,1].set_ylabel("MS1 Clusters\nin Features (%)", fontsize=ylabel_size)
    axs[1,1].set_yticks(np.arange(0, 101, 20), labels=np.arange(0,101,20), fontsize=ytick_size)
    axs[1,1].set_xticks(np.arange(0, len(algorithms)), labels=algorithms, fontsize=xlabel_size)

    ## endregion MS1 Clusters in Features

    ##region Features Matched
    axs[2,0].bar(algorithms, features_matched, yerr=features_matched_std, color=colors, capsize=5, edgecolor='black',
            clip_on=False)

    ## region Add significance brackets and stars
    axs[2,0].plot([0, 0, 2, 2],
             [features_matched[0] + features_matched_std[0] + 15, features_matched[0] + features_matched_std[0] + 17,
              features_matched[0] + features_matched_std[0] + 17, features_matched[0] + features_matched_std[0] + 15],
             color='black', linewidth=1, clip_on=False)
    axs[2,0].text(1, features_matched[0] + features_matched_std[0] + 16, "*", fontsize=14, ha='center')

    axs[2,0].plot([0, 0, 3, 3],
             [features_matched[0] + features_matched_std[0] + 21, features_matched[0] + features_matched_std[0] + 23,
              features_matched[0] + features_matched_std[0] + 23, features_matched[0] + features_matched_std[0] + 21],
             color='black', linewidth=1, clip_on=False)
    axs[2,0].text(1.5, features_matched[0] + features_matched_std[0] + 22, "*", fontsize=14, ha='center')

    axs[2,0].plot([1, 1, 2, 2],
             [features_matched[0] + features_matched_std[0] + 1, features_matched[0] + features_matched_std[0] + 3,
              features_matched[0] + features_matched_std[0] + 3, features_matched[0] + features_matched_std[0] + 1],
             color='black', linewidth=1, clip_on=False)
    axs[2,0].text(1.5, features_matched[0] + features_matched_std[0] + 1.5, "*", fontsize=14, ha='center')

    axs[2,0].plot([1, 1, 3, 3],
             [features_matched[0] + features_matched_std[0] + 9, features_matched[0] + features_matched_std[0] + 11,
              features_matched[0] + features_matched_std[0] + 11, features_matched[0] + features_matched_std[0] + 9],
             color='black', linewidth=1, clip_on=False)
    axs[2,0].text(1.9, features_matched[0] + features_matched_std[0] + 10, "*", fontsize=14, ha='center')

    ## endregion Add significance brackets and stars

    axs[2,0].spines['top'].set_visible(False)
    axs[2,0].spines['right'].set_visible(False)
    axs[2,0].set_ylabel("Features Matched (%)", fontsize=ylabel_size)
    axs[2,0].set_yticks(np.arange(0, 101, 20), labels=np.arange(0,101,20), fontsize=ytick_size)
    axs[2,0].set_xticks(np.arange(0, len(algorithms)), labels=algorithms, fontsize=xlabel_size)
    axs[2,0].set_ylim(0, 100)

    ## endregion Features Matched

    ##region Feature Intensity Matched
    axs[2,1].bar(algorithms, intensity_matched, yerr=intensity_matched_std,
            color=colors, capsize=5, edgecolor='black', clip_on=False)

    ## region Add significance brackets and stars
    axs[2,1].plot([0, 0, 2, 2],
             [intensity_matched[0] + intensity_matched_std[0] + 7, intensity_matched[0] + intensity_matched_std[0] + 8,
              intensity_matched[0] + intensity_matched_std[0] + 8, intensity_matched[0] + intensity_matched_std[0] + 7],
             color='black', linewidth=1, clip_on=False)
    axs[2,1].text(1, intensity_matched[0] + intensity_matched_std[0] + 7, "*", fontsize=14, ha='center')

    axs[2,1].plot([1, 1, 2, 2],
             [intensity_matched[0] + intensity_matched_std[0] + 1, intensity_matched[0] + intensity_matched_std[0] + 2,
              intensity_matched[0] + intensity_matched_std[0] + 2, intensity_matched[0] + intensity_matched_std[0] + 1],
             color='black', linewidth=1, clip_on=False)
    axs[2,1].text(1.5, intensity_matched[0] + intensity_matched_std[0] + 0.75, "*", fontsize=14, ha='center')

    ## endregion Add significance brackets and stars

    axs[2,1].spines['top'].set_visible(False)
    axs[2,1].spines['right'].set_visible(False)
    axs[2,1].set_ylabel("Feature Intensity\nMatched (%)", fontsize=ylabel_size, ha='center')
    axs[2,1].set_yticks(np.arange(0, 101, 20), labels=np.arange(0,101,20), fontsize=ytick_size)
    axs[2,1].set_xticks(np.arange(0, len(algorithms)), labels=algorithms, fontsize=xlabel_size)
    axs[2,1].set_ylim(0, 100)
    ## endregion Feature Intensity Matched

    #Add letter labels to each subplot in the upper right corner
    xpos = 0.9
    ypos = 1

    axs[0,0].text(xpos,ypos, "A", transform=axs[0,0].transAxes, fontsize=30)
    axs[0,1].text(xpos,ypos, "B", transform=axs[0,1].transAxes, fontsize=30)
    axs[1,0].text(xpos,ypos, "C", transform=axs[1,0].transAxes, fontsize=30)
    axs[1,1].text(xpos,ypos, "D", transform=axs[1,1].transAxes, fontsize=30)
    axs[2,0].text(xpos,ypos, "E", transform=axs[2,0].transAxes, fontsize=30)
    axs[2,1].text(xpos,ypos, "F", transform=axs[2,1].transAxes, fontsize=30)

    plt.tight_layout()
    plt.savefig(
        "D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\isodec_td_decon.png",
        dpi=300)
    plt.savefig(
        "D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\isodec_td_decon.pdf",
        dpi=300)
    plt.show()