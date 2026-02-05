import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from sklearn.neighbors import LocalOutlierFactor
# from Scripts.MTM.Lipidomics.LipidFunctions import encode_headgroups
import matplotlib as mpl
import mplcursors
from PublicScripts.Lipids.LipidFunctions import detect_outliers_class

mpl.use("WxAgg")

def get_labels(df):
    labels = []
    for i, row in df.iterrows():
        string = (f"{row['Metabolite name']}\n"
                  f"dRT: {row['RTdiff']:.2f} min\nScore: {row['Total score']}\nS/N: {row['S/N average']}")

        # Add which contrib is highest
        contrib_cols = [col for col in df.columns if col.endswith('_contrib')]
        max_contrib = 0
        max_feature = ""
        for col in contrib_cols:
            if row[col] > max_contrib:
                max_contrib = row[col]
                max_feature = col.replace('_contrib', '')
        outlier = row['outlier']
        if outlier:
            string += f"\nOutlier Factor: {max_feature}"

        labels.append(string)
    return labels


os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")
posfile = "PosIDsRTP2_HeadgroupTailCheck_Network.csv"
negfile = "NegIDsRTP2_HeadgroupTailCheck_Network.csv"
posdf = pd.read_csv(posfile)
negdf = pd.read_csv(negfile)

# print(posdf.columns)
rtdiff = posdf['Average Rt(min)'] - posdf['Reference RT']
posdf['RTdiff'] = rtdiff
nrtdiff = negdf['Average Rt(min)'] - negdf['Reference RT']
negdf['RTdiff'] = nrtdiff
posdf["LogSN"] = np.log10(posdf["S/N average"] + 1)
negdf["LogSN"] = np.log10(negdf["S/N average"] + 1)
# Setup headgroup encoding
#'HG_Fragment_Match', 'HG_NL_Match', 'Tail_Fragment_Match', 'Tail_NL_Match'
features = ['RTdiff', 'Total score', 'LogSN', 'HG_Fragment_Match', 'HG_NL_Match', 'Tail_Fragment_Match', 'Tail_NL_Match', "z_Anomaly_Score"]

posdf = detect_outliers_class(posdf, features)
negdf = detect_outliers_class(negdf, features)



posdf.to_csv("PosIDsRTP2_OutlierAnalysis.csv", index=False)
negdf.to_csv("NegIDsRTP2_OutlierAnalysis.csv", index=False)

x = 'Average Rt(min)'
y = 'Average Mz'

plt.figure(figsize=(12, 6))
plt.subplot(1, 4, 1)

plt.plot(posdf["mahalanobis"].values, marker='o', linestyle='', label='Positive Mode')
plt.plot(negdf["mahalanobis"].values, marker='x', linestyle='', label='Negative Mode')
# plt.axhline(y=3, color='r', linestyle='--', label='Outlier Threshold')
plt.xlabel('Index')
plt.ylabel('Negative Outlier Factor Score')
plt.title('Outlier Scores for Lipid Identifications')
plt.legend()

plt.subplot(142)


sc1 = plt.scatter(posdf[x], posdf[y], c=posdf['outlier'], cmap='coolwarm', label='Positive Mode')

cursor = mplcursors.cursor(sc1, hover=True)
poslabels = get_labels(posdf)

neglabels = get_labels(negdf)
@cursor.connect("add")
def on_add(sel):
    # sel.index is the point index for scatter
    sel.annotation.set_text(poslabels[sel.index])
    sel.annotation.get_bbox_patch().set_alpha(0.9)

plt.xlabel(x)
plt.ylabel(y)
plt.title('Positive Mode Outlier Detection')
plt.subplot(1, 4, 3)
sc2 = plt.scatter(negdf[x], negdf[y], c=negdf['outlier'], cmap='coolwarm', label='Negative Mode')
cursor2 = mplcursors.cursor(sc2, hover=True)
@cursor2.connect("add")
def on_add2(sel):
    # sel.index is the point index for scatter
    sel.annotation.set_text(neglabels[sel.index])
    sel.annotation.get_bbox_patch().set_alpha(0.9)
plt.xlabel(x)
plt.ylabel(y)
plt.title('Negative Mode Outlier Detection')

# Create a histogram of top contrib feature for only outliers
plt.subplot(1, 4, 4)
all = pd.concat([posdf, negdf])
outliers = all[all['outlier'] == True]
all_top_features = outliers['Top contrib feature']
top_feature_counts = all_top_features.value_counts()
top_feature_counts.plot(kind='bar', color='skyblue')
plt.xlabel('Feature')
plt.ylabel('Count')
plt.title('Top Contributing Features for Outliers')


plt.tight_layout()
plt.show()

