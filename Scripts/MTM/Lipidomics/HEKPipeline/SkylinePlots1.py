import matplotlib.pyplot as plt

from PublicScripts.Lipids.SkylineResultsPlotter import *

file = r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs\Molecule Transition Results_stats_working.xlsx"
file = r"Z:\Group Share\Annika\Stellar\FAM\CB2 Discs\Molecule Transition Results_stats.xlsx"

classdf = pd.read_excel(file, sheet_name="Class Stats")
normdf = pd.read_excel(file, sheet_name="Molecule Stats")
print(normdf)

set1 = ["FT 1", "FT 2", "FT 3"]
set2 = ["E 1", "E 2", "E 3"]
set3 = ["Memb Extr 1", "Memb Extr 2", "Memb Extr 3"]

# Figure with pie 1 volcano pie 2 in a row
# plt.figure(figsize=(18, 6))
# ax1 = plt.subplot(1, 3, 1)
# pie_chart(classdf, "Set1_Mean", ax=ax1, title="Depleted Nanodiscs", otherthresh=0.015)
# ax2 = plt.subplot(1, 3, 2)
# make_volcano_plot(normdf, title="Depleted ND vs. CB2 ND", ax=ax2, fold_range=np.log2(1.75), xlims=[-6, 6])
# ax3 = plt.subplot(1, 3, 3)
# pie_chart(classdf, "Set2_Mean", ax=ax3, title="CB2 Nanodiscs", otherthresh=0.015)
# plt.tight_layout()
# plt.savefig(os.path.splitext(file)[0] + "_volcano_pies.png")
# plt.savefig(os.path.splitext(file)[0] + "_volcano_pies.pdf")
# plt.show()

# Class plots bar chart
# plt.figure(figsize=(10, 6))
# ax1 = plt.subplot(1, 1, 1)
# make_class_plots(classdf, ax=ax1)
# plt.tight_layout()
# plt.savefig(os.path.splitext(file)[0] + "_class_bar.png")
# plt.savefig(os.path.splitext(file)[0] + "_class_bar.pdf")
# plt.show()

# Heatmap of the fold changes for each molecule in the volcano plot that passes thresholds
# plt.figure(figsize=(5, 10))
# ax1 = plt.subplot(1, 1, 1)
# normdfsorted = filter_normdf(normdf, fold_range=np.log2(1.75))
# heatmap_plot(normdfsorted, ["Log2FC", "StdDevFC"], title="Fold Change Heatmap", ax=ax1, lindex=0, colormap="RdBu_r",
#              normtype="zerocentered", fontsize=9, color_classes=True, colorbar="Log2 Fold Change")
# plt.tight_layout()
# plt.savefig(os.path.splitext(file)[0] + "_heatmap.png", dpi=600)
# plt.savefig(os.path.splitext(file)[0] + "_heatmap.pdf")
# plt.show()


# uclasses = normdf["Molecule List Name"].unique()
# for c in uclasses:
#     plt.figure(figsize=(10, 6))
#     ax = plt.subplot(1, 1, 1)
#     lipid_bar_chart_compare(normdf, classes=[c], ax=ax, fontsize=10, s1name="Depleted Nanodiscs", s2name="CB2 Nanodiscs", log2fclimit=np.log2(1))
#     plt.tight_layout()
#     plt.show()
# exit()
# plt.figure(figsize=(10, 6))
# ax = plt.subplot(1, 1, 1)
# lipid_bar_chart_compare(normdf, classes=["SM"], ax=ax, fontsize=10, s1name="Depleted Nanodiscs", s2name="CB2 Nanodiscs", log2fclimit=np.log2(1))
# plt.tight_layout()
# # plt.savefig(os.path.splitext(file)[0] + "_PE_bar.png", dpi=600)
# # plt.savefig(os.path.splitext(file)[0] + "_PE_bar.pdf")
# plt.show()

# Stacked pie chart to show the membrane
plt.figure(figsize=(6, 8))
ax1 = plt.subplot(2, 1, 1)
pie_chart(classdf, "Set1_Mean", ax=ax1, title="HEK Membrane Extract", otherthresh=0.015)
ax2 = plt.subplot(2, 1, 2)
make_volcano_plot(normdf, title="Membrane vs. CB2 ND", ax=ax2, fold_range=np.log2(1.75), xlims=[-6, 6])
# ax3 = plt.subplot(1, 3, 3)
# pie_chart(classdf, "Set2_Mean", ax=ax3, title="CB2 Nanodiscs", otherthresh=0.015)
plt.tight_layout()
plt.savefig(os.path.splitext(file)[0] + "_volcano_pies_2.png")
plt.savefig(os.path.splitext(file)[0] + "_volcano_pies_2.pdf")
plt.show()