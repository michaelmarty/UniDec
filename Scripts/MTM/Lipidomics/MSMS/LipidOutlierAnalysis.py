import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from sklearn.neighbors import LocalOutlierFactor
# from Scripts.MTM.Lipidomics.LipidFunctions import encode_headgroups
from sklearn.covariance import EmpiricalCovariance, MinCovDet
import matplotlib as mpl
from scipy.spatial.distance import mahalanobis
import mplcursors

mpl.use("WxAgg")

# def detect_outliers(df, features, n_neighbors=20, contamination="auto"):
#     lof = LocalOutlierFactor(n_neighbors=n_neighbors, contamination=contamination)
#     df['outlier'] = lof.fit_predict(df[features])
#     df['outlier_score'] = lof.negative_outlier_factor_
#
#     # Sort by outlier score
#     df = df.sort_values(by='outlier_score').reset_index(drop=True)
#
#     print("Outlier detection complete. Number of outliers detected:", sum(df['outlier'] == -1))
#     print(df[df['outlier'] == -1][["Metabolite name"] + features + ['outlier_score']])
#     return df

def detect_outliers_class(df, features):
    newdf = pd.DataFrame()
    for c in df['Ontology'].unique():
        # Get subset for class
        subdf = df[df['Ontology'] == c]

        # Clean up things that will cause issues with covariance matrix
        if len(subdf) < 3:
            newdf = pd.concat([newdf, subdf])
            continue

        # Drop any features with zero variance
        variances = subdf[features].var()
        valid_features = variances[variances > 0].index.tolist()
        # Sort features by variance
        valid_features = sorted(valid_features, key=lambda x: variances[x], reverse=True)

        # Check for two columns of valid features that are identical
        to_drop = set()
        for i in range(len(valid_features)):
            for j in range(i+1, len(valid_features)):
                if np.all(subdf[valid_features[i]] == subdf[valid_features[j]]):
                    to_drop.add(valid_features[j])
        valid_features = [f for f in valid_features if f not in to_drop]

        # Remove class if not enough valid features
        if len(subdf) < len(valid_features)+1:
            valid_features = valid_features[:len(subdf)-1]

        outputdf = detect_outliers_mahalanobis(subdf, valid_features)
        newdf = pd.concat([newdf, outputdf])
    df = newdf.reset_index(drop=True)
    # Sort by outlier score
    df = df.sort_values(by='mahalanobis').reset_index(drop=True)

    # Print number of outliers detected
    print("Outlier detection complete. Number of outliers detected:", np.nansum(df['outlier']), "out of", len(df))
    return df


def detect_outliers_mahalanobis(df, features, threshold=10.0):
    # Calculate Mahalanobis distance for each point
    data = df[features].values
    mean = np.mean(data, axis=0)
    cov = np.cov(data, rowvar=False)
    try:
        inv_cov = np.linalg.inv(cov)
    except:
        inv_cov = np.linalg.pinv(cov)

    outdf = df.copy()
    outdf['mahalanobis'] = [mahalanobis(row, mean, inv_cov) for row in data]
    outdf['outlier'] = outdf['mahalanobis'] > threshold

    #Robust Mahalanobis
    mcd = MinCovDet().fit(data)
    robust_mean = mcd.location_
    robust_cov = mcd.covariance_
    try:
        robust_inv_cov = np.linalg.inv(robust_cov)
    except:
        robust_inv_cov = np.linalg.pinv(robust_cov)

    robust_mahalanobis = [mahalanobis(row, robust_mean, robust_inv_cov) for row in data]
    outdf['mahalanobis'] = robust_mahalanobis
    outdf['outlier'] = outdf['mahalanobis'] > threshold

    eigvals, V = np.linalg.eigh(cov)
    eigvals = np.maximum(eigvals, 1e-12)  # Prevent division by zero
    z = V.T @ (data - mean).T
    y = z.T / np.sqrt(eigvals)
    contrib = np.square(y)
    # Add contribution columns to outdf with feature_contrib name
    for i, feature in enumerate(features):
        outdf[f'{feature}_contrib'] = contrib[:, i]
    # Sort by outlier score
    outdf = outdf.sort_values(by='mahalanobis').reset_index(drop=True)

    # if sum(outdf['outlier']) != 0:
    #     print("Mahalanobis outlier detection complete. Number of outliers detected:", sum(outdf['outlier']))
    #     print(outdf[outdf['outlier']][["Metabolite name"] + features + ['mahalanobis']])
    return outdf

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
posfile = "PosIDsRTP2_HeadgroupTailCheck.csv"
negfile = "NegIDsRTP2_HeadgroupTailCheck.csv"
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
features = ['RTdiff', 'Total score', 'LogSN', 'HG_Fragment_Match', 'HG_NL_Match']

posdf = detect_outliers_class(posdf, features)
negdf = detect_outliers_class(negdf, features)



posdf.to_csv("PosIDsRTP2_OutlierAnalysis.csv", index=False)
negdf.to_csv("NegIDsRTP2_OutlierAnalysis.csv", index=False)

x = 'Average Rt(min)'
y = 'Average Mz'

plt.figure(figsize=(12, 6))
plt.subplot(1, 3, 1)

plt.plot(posdf["mahalanobis"].values, marker='o', linestyle='', label='Positive Mode')
plt.plot(negdf["mahalanobis"].values, marker='x', linestyle='', label='Negative Mode')
plt.axhline(y=-1.5, color='r', linestyle='--', label='Outlier Threshold')
plt.xlabel('Index')
plt.ylabel('Negative Outlier Factor Score')
plt.title('Outlier Scores for Lipid Identifications')
plt.legend()

plt.subplot(132)


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
plt.subplot(1, 3, 3)
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
plt.tight_layout()
plt.show()

