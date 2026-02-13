import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from PublicScripts.Lipids.LipidOutlierAnalysis import dataset_split_rf_data, random_forest, plot_rf_results, \
    prep_for_rf, plsr, random_split_rf_data, apply_rf_prediction
import joblib

mpl.use("WxAgg")

topdir = r"Z:\Group Share\Shivam\DDA"
os.chdir(topdir)

filepath = r"Combined_DDA_OutlierAnalysis.xlsx"
df = pd.read_excel(filepath)
# Drop any with "AqX" in the dataset name, since those are all not annotated and will just add noise to the model
df = df[~df["Dataset"].str.contains("AqX", case=False, na=False)]

drop_columns = ["Post curation result", "MS/MS assigned", "Annotation tag (VS1.0)", "INCHIKEY", "m/z matched",
                "MS/MS matched", "Manually modified for quantification", "Manually modified for annotation",
                "Isotope tracking parent ID", "Isotope tracking weight number", "Spectrum reference file name"]

# Drop columns
for col in drop_columns:
    if col in df.columns:
        df = df.drop(columns=[col])

features = ['RTdiff', 'LogSN', 'Headgroup_Match', 'Tail_Match_Fraction',
                "Class_z_Anomaly_Score", "Tail_z_Anomaly_Score","mahalanobis", "Simple dot product", "Fill %", "RT matched"
                "Weighted dot product", "Reverse dot product", "Match_Grade_Encoded", "Headgroup Encoding","Adduct Encoding"]

not_features = ["Metabolite name", "Formula", "Ontology", "SMILES", "Comment", "MS1 isotopic spectrum",
                "MS/MS spectrum", "Ref Spec", "Polarity", "Dataset", "Dataset Count", "Tails", "T1", "T2", "T3",
                "T4", "Low_Quality_Flag", 'Adduct type', "Match_Grade"]

drop_features = ['Alignment ID', 'Average Rt(min)', 'Average Mz', 'Reference RT', 'Reference m/z', "RT matched",
                 "Adduct_Match", "outlier", "Fill %", "S/N average", "Matched peaks percentage",
                 "m/z similarity", "RT similarity", "RTdiff", "Class_Anomaly_Score", "Tail_Anomaly_Score",
                 "Tail_z_Anomaly_Score", "Class_z_Anomaly_Score", "Tail_Match_Fraction", 'Reverse dot product',
                 "LogSN", "mahalanobis", "Total score", "Headgroup_Match", "Weighted dot product", "Class_degree"]
drop_features = []
columns = df.columns.tolist()
features = []
for col in columns:
    if col not in not_features and col not in drop_features:
        if ("Class_z_" not in col):
            features.append(col)

features = ["Simple dot product", "Matched peaks count", "Class_strength"]

df, features = prep_for_rf(df, features, one_hot_encode=False)

print("Features:", features)

fimps = []
accuracies = []
fimps2 = []
accuracies2 = []
savemodel = None
for dataset in df["Dataset"].unique():
    # print(f"Processing dataset: {dataset}")
    x_train, x_test, y_train, y_test = dataset_split_rf_data(df, features, "Low_Quality_Flag",
                                                             test_dataset=None, training_dataset=dataset)
    # x_train, x_test, y_train, y_test = dataset_split_rf_data(df, features, "Low_Quality_Flag",
    #                                                          test_dataset=dataset, training_dataset=None)
    # x_train, x_test, y_train, y_test = random_split_rf_data(df, features, "Low_Quality_Flag", random_state=None)
    plsrmodel, feature_importance2, y_pred2, report2 = plsr(x_train, x_test, y_train, y_test, features)
    clf, feature_importance, y_pred, report = random_forest(x_train, x_test, y_train, y_test, features, n_estimators=25)

    if dataset == "Stellar":
        savemodel = clf
    accuracies.append(report["accuracy"])
    print(f"Accuracy for {dataset}: {report['accuracy']:.3f}")
    # print("Feature Importances:")
    # print(feature_importance)
    # plot_rf_results(feature_importance)
    # Save model for future use
    fimps.append(feature_importance)
    fimps2.append(feature_importance2)
    accuracies2.append(report2["accuracy"])
print("Average Accuracy: {:.3f}".format(np.mean(accuracies)))
print("Average Accuracy PLSR: {:.3f}".format(np.mean(accuracies2)))
# Plot average feature importance across all datasets
avg_fimp = pd.concat(fimps, axis=1).mean(axis=1).sort_values(ascending=False)
avg_fimp2 = pd.concat(fimps2, axis=1).mean(axis=1).sort_values(ascending=False)

df = apply_rf_prediction(df, savemodel)
# Write all to sheet and then split by "Outlier_Prediction into two additional sheets
with pd.ExcelWriter("Combined_DDA_OutlierAnalysis_Results.xlsx") as writer:
    df.to_excel(writer, sheet_name="All Data", index=False)
    df[df["Outlier_Prediction"] == 0].to_excel(writer, sheet_name="Predicted Good", index=False)
    df[df["Outlier_Prediction"] == 1].to_excel(writer, sheet_name="Predicted Bad", index=False)
print(df)


# Feature importance plots
plt.figure(figsize=(10, 6))
ax = plt.subplot(1, 2, 1)
plot_rf_results(avg_fimp, ax=ax)
ax = plt.subplot(1, 2, 2)
plot_rf_results(avg_fimp2, ax=ax)
plt.tight_layout()

# Plots of distribution of each feature for good vs bad IDs
plt.figure(figsize=(10, 6))
for i, feature in enumerate(features):
    ax = plt.subplot(1, len(features), i+1)
    for label in df["Low_Quality_Flag"].unique():
        subset = df[df["Low_Quality_Flag"] == label]
        ax.hist(subset[feature], bins=20, alpha=0.5, label=str(label))
    ax.set_title(feature)
    ax.legend()
plt.tight_layout()



plt.show()

joblib.dump(clf, "RF_Model.pkl")