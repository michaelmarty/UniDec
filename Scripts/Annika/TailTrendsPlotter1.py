import os
import PublicScripts.Lipids.LipidFunctions as LipidFunctions
import pandas as pd
# python
import numpy as np
import matplotlib.pyplot as plt

# matplotlib.use("WxAgg")

os.chdir(r"Z:\Group Share\Annika\Stellar\Untargeted DDA")

file = r"PeakMaster_2_2025_11_10_14_39_41.csv"
column_name = "Metabolite name"
rt_name = "Average Rt(min)"
class_name = "Ontology"
refcol = "Reference RT"

# os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")
#
# file = r"C1_Ref_Testing_FilteredHighLimited10.csv"
# column_name = "Name"
# rt_name = "RT"
# class_name = "CompoundClass"
# refcol = "RTP"
# refcol= "Reference RT"

df = pd.read_csv(file)
# Remove any with no MS2 in molecule name
df = df[~df[column_name].str.contains("no MS2", na=False)]
# Remove any with lowscore in name
df = df[~df[column_name].str.contains("low score", na=False)]
# Take only ones where MS/MS assigned is TRUE
df = df[df["MS/MS assigned"] == True]
df = df[df["MS/MS matched"] == True]
# Take only with total score above 1.5
df = df[df["Total score"] >= 1.5]


df = LipidFunctions.clean_df(df, column=column_name)
df = LipidFunctions.set_tails(df, columnname=column_name)

print(np.average(df[rt_name] - df[refcol])),
# plt.hist(df[rt_name] - df[refcol], bins=50)
# plt.show()

# plt.plot(df[refcol], df[rt_name], 'o')
# plt.xlabel('Reference RT (min)')
# plt.ylabel('Average RT (min)')
# plt.plot([0, max(df[refcol])], [0, max(df[refcol])], 'r--')
# plt.show()
# exit()

# df = LipidFunctions.encode_headgroups(df, classname=class_name, rtcol_name=rt_name)
df["Residuals"] = np.abs(df[rt_name] - df[refcol])

# Remove residuals above 0.25 min
filtereddf = df[df["Residuals"] <= 2]
# Filter to remove tail length -1
filtereddf = filtereddf[filtereddf["Tail Length"] >= 0]
(print(len(filtereddf)))
# Write to file
if "Filtered" not in file:
    filtereddf.to_csv("Filtered_" + file, index=False)
df = filtereddf

# Stats on residuals
print(f'Mean Residual: {np.mean(df["Residuals"]):.3f} min')
print(f'Median Residual: {np.median(df["Residuals"]):.3f} min')
print(f'STD of Residuals: {np.std(df["Residuals"]):.3f} min')

classes = df[class_name].unique()
print(classes)

# classes = ["PC", "PE", "TG"]
# classes = ["LPC"]

for c in classes:
    # Isolate FA ontology
    fa_df = df[df[class_name] == c]
    if len(fa_df) < 20:
        continue
    #
    # print(f'Mean Residual: {np.mean(fa_df["Residuals"]):.3f} min')
    # print(f'STD of Residuals: {np.std(fa_df["Residuals"]):.3f} min')
    plt.figure()
    classstd = np.std(fa_df["Residuals"])

    unsats = np.unique(fa_df["Tail Unsaturation"].values)
    colors = plt.cm.viridis(np.linspace(0, 1, len(unsats)))

    for i, u in enumerate(unsats):
        subdf = fa_df[fa_df["Tail Unsaturation"] == u]
        # Sort subdf by tail lengths
        subdf = subdf.sort_values(by=["Tail Length"])
        tail_lengths = subdf["Tail Length"].values
        rts = subdf[rt_name].values
        refrts = subdf[refcol].values

        # plt.subplot(122)
        plt.plot(tail_lengths, refrts, marker='x', color=colors[i])

        #
        # plt.subplot(121)
        plt.plot(tail_lengths, rts, marker='o', linestyle="--", label=f'{c} Unsaturation {u}', color=colors[i])

        # Mark outliers
        for j in range(len(rts)):
            if np.abs(rts[j] - refrts[j]) > 0.5:
                plt.plot(tail_lengths[j], rts[j], marker='s', markersize=10, markerfacecolor='none', markeredgecolor='r', markeredgewidth=2)


    plt.xlabel('Tail Length')
    plt.ylabel('Retention Time (min)')
    plt.legend()
plt.show()