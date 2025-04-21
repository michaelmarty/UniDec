import matplotlib.pyplot as plt
import tools as ud
from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 12

path = "C:\\Data\\Lipidomics\\AmtB_LXMS.xlsx"

df = pd.read_excel(path)

df = set_tails(df, "Name")
df = set_class_name(df, "Name")

df = calc_tail_diff(df)

unique_classes = np.unique(df["Class Name"])
exclude = ["LPE"]
colordict = {"CL": "orange", "PG": "green", "PE": "purple"}
figure = plt.figure(figsize=(4, 8))
axes = figure.subplots(3,1)
for c in unique_classes:
    if c not in exclude:
        subdf = df[df["Class Name"] == c]
        # subdf = subdf[subdf["Flag_A_AmtB_vs_Ctrl"] != 0]

        n = len(subdf)
        tval = ud.get_tvalue(n - 1)
        fc = subdf["Fold_change_AmtB_vs_Ctrl"]

        # plt.plot(subdf["Tail Unsaturation"], subdf["Fold_change_AqpZ_vs_Ctrl"], "o")

        avg_fold_change = np.mean(fc)
        std_fold_change = np.std(fc)
        axes[0].bar(c, avg_fold_change, yerr=std_fold_change * tval / np.sqrt(n), capsize=5, color=colordict[c])

        avg_ldiff = np.average(subdf["Tail Length Diff"], weights=fc)
        std_ldiff = ud.weighted_std(subdf["Tail Length Diff"], fc)
        axes[1].bar(c, avg_ldiff, yerr=std_ldiff * tval / np.sqrt(n), capsize=5, color=colordict[c])

        avg_udiff = np.average(subdf["Tail Unsat Diff"], weights=fc)
        std_udiff = ud.weighted_std(subdf["Tail Unsat Diff"], fc)
        axes[2].bar(c, avg_udiff, yerr=std_udiff * tval / np.sqrt(n), capsize=5, color=colordict[c])

axes[0].set_ylabel("Avg. Log2Fold Change")
axes[1].set_ylabel("Avg. Lipid Length Diff.")
axes[2].set_ylabel("Avg. Lipid Unsat. Diff.")
# Remove top and right axis lines
for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # Remove tick marker but leave label
    ax.tick_params(axis='x', which='both', length=0)
    # Add horizontal line at y=0
    ax.axhline(y=0, color='k', linestyle='-')

# Save Fig
outpath = path[:-5] + "_out.png"
plt.tight_layout()
plt.savefig(outpath, dpi=300, bbox_inches='tight')


plt.show()

# exit()
outpath = path[:-5] + "_out.xlsx"
df.to_excel(outpath)
