import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import LocalOutlierFactor
from Scripts.MTM.Lipidomics.LipidFunctions import encode_headgroups

def detect_outliers(df, features, n_neighbors=20, contamination=0.05):
    lof = LocalOutlierFactor(n_neighbors=n_neighbors, contamination=contamination)
    df['outlier'] = lof.fit_predict(df[features])
    df['outlier_score'] = lof.negative_outlier_factor_

    print("Outlier detection complete. Number of outliers detected:", sum(df['outlier'] == -1))
    print(df[df['outlier'] == -1][["Metabolite name"] + features + ['outlier_score']])
    return df

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")
posfile = "PosIDsRTP2_HeadgroupTailCheck.csv"
negfile = "NegIDsRTP2_HeadgroupTailCheck.csv"
posdf = pd.read_csv(posfile)
negdf = pd.read_csv(negfile)

print(posdf.columns)
rtdiff = posdf['Average Rt(min)'] - posdf['Reference RT']
posdf['RTdiff'] = rtdiff
nrtdiff = negdf['Average Rt(min)'] - negdf['Reference RT']
negdf['RTdiff'] = nrtdiff
posdf["LogSN"] = np.log10(posdf["S/N average"] + 1)
negdf["LogSN"] = np.log10(negdf["S/N average"] + 1)
# Setup headgroup encoding
posdf = encode_headgroups(posdf)
negdf = encode_headgroups(negdf)

features = ['Headgroup Encoding', 'RTdiff', 'Total score', 'LogSN', 'HG_Fragment_Match', 'HG_NL_Match', 'Tail_Fragment_Match', 'Tail_NL_Match']

posdf = detect_outliers(posdf, features)
negdf = detect_outliers(negdf, features)

posdf.to_csv("PosIDsRTP2_OutlierAnalysis.csv", index=False)
negdf.to_csv("NegIDsRTP2_OutlierAnalysis.csv", index=False)

x = 'Average Rt(min)'
y = 'Average Mz'

plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.scatter(posdf[x], posdf[y], c=posdf['outlier'], cmap='coolwarm', label='Positive Mode')
plt.xlabel(x)
plt.ylabel(y)
plt.title('Positive Mode Outlier Detection')
plt.subplot(1, 2, 2)
plt.scatter(negdf[x], negdf[y], c=negdf['outlier'], cmap='coolwarm', label='Negative Mode')
plt.xlabel(x)
plt.ylabel(y)
plt.title('Negative Mode Outlier Detection')
plt.tight_layout()
plt.show()

