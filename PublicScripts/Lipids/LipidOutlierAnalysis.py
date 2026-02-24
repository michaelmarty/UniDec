import os
import pandas as pd
from PublicScripts.Lipids.LipidFunctions import (add_ref_spec_xml, class_network_analysis,
                                                 tail_network_analysis, merge_pos_neg, detect_outliers_class,
                                                 get_unique_names, encode_headgroups)
from PublicScripts.Lipids.LipidFragPrediction import assign_df_fragments
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import PLSRegression
import joblib

drop_columns = ["MS/MS assigned", "m/z matched", "Manually modified for annotation", "RT matched", "MS/MS matched",
                "Post curation result", "Manually modified for quantification",
                "Isotope tracking parent ID", "Isotope tracking weight number", "Spectrum reference file name",
                "INCHIKEY", "Annotation tag (VS1.0)",
                "Matched peaks percentage", "Matched peaks count","Fill %",
                "RT similarity", "m/z similarity", "Formula"]

def comment_parser(df):
    if "Comment" not in df.columns:
        return df
    if df["Comment"].isnull().all():
        return df
    # See if any of the comments contain "Low quality"
    ucomments = df["Comment"].unique()
    if not any("Low quality" in str(c) for c in ucomments):
        return df
    # Create a new column "Low_Quality_Flag"
    df["Low_Quality_Flag"] = df["Comment"].apply(lambda x: True if "Low quality" in str(x) else False)
    return df

def basic_class_network_analysis(df, tol=0.3):
    posdf = df[df["Polarity"] == "Positive"]
    negdf = df[df["Polarity"] == "Negative"]

    posdf, posgraphs = class_network_analysis(posdf, tol=tol, add_all_cols=True)
    negdf, neggraphs = class_network_analysis(negdf, tol=tol, add_all_cols=True)

    # Recombine DFs
    combined_df = pd.concat([posdf, negdf], ignore_index=True)

    return combined_df


def outlier_analysis(posdf, negdf, tol=0.05, rttol=0.1, contribs=False, add_all_cols_heads=False,
                     add_all_cols_tails=False, drop_cols=True):
    if drop_cols:
        # Drop unnecessary columns if present
        for col in drop_columns:
            if col in posdf.columns:
                posdf = posdf.drop(columns=[col])
            if col in negdf.columns:
                negdf = negdf.drop(columns=[col])

    print("Headgroup and Tail Fragment Checks...")
    # Perform Tail and Headgroup Fragment Checks
    posdf = assign_df_fragments(posdf, tol=tol)
    negdf = assign_df_fragments(negdf, tol=tol)

    print("MSMS Networking Analysis...")
    posdf, posgraphs = class_network_analysis(posdf, tol=tol, add_all_cols=add_all_cols_heads)
    negdf, neggraphs = class_network_analysis(negdf, tol=tol, add_all_cols=add_all_cols_heads)
    # headgraphs = posgraphs + neggraphs
    #
    posdf, posgraphs = tail_network_analysis(posdf, min_count=5, mz_tol=tol, add_all_cols=add_all_cols_tails)
    negdf, neggraphs = tail_network_analysis(negdf, min_count=5, mz_tol=tol, add_all_cols=add_all_cols_tails)
    # tailgraphs = posgraphs + neggraphs

    # Merge DFs
    posdf, negdf = merge_pos_neg(posdf, negdf, rttol)

    # Combined Pos and Neg
    combined_df = pd.concat([posdf, negdf], ignore_index=True)

    # Outlier Detection
    rtdiff = combined_df['Average Rt(min)'] - combined_df['Reference RT']
    combined_df['RTdiff'] = rtdiff
    combined_df["LogSN"] = np.log10(combined_df["S/N average"] + 1)

    features = ['RTdiff', 'Total score', 'LogSN', 'Headgroup_Match', 'Tail_Match_Fraction',
                "Class_z_Anomaly_Score", "Tail_z_Anomaly_Score"]
    combined_df = detect_outliers_class(combined_df, features, separate_adducts=True, contribs=contribs)

    if not add_all_cols_tails:
        # Drop all columns that start with "T1" or so on
        cols_to_drop = [col for col in combined_df.columns if col.startswith(('T1', 'T2', 'T3', 'T4', 'T5'))]
        combined_df = combined_df.drop(columns=cols_to_drop)

    # Sort by Class and then Metabolite Name
    combined_df = combined_df.sort_values(by=["Ontology", "Metabolite name"])

    # Parse Comments if present
    combined_df = comment_parser(combined_df)
    return combined_df

def outlier_pipeline(posfile, negfile, posxml="", negxml="", tol=0.05):
    print("Importing Files:", posfile, negfile)
    posdf = pd.read_csv(posfile)
    negdf = pd.read_csv(negfile)

    print("Adding Reference Spectra from XML if needed...")
    # Add in Reference Spectra if not already present
    if posxml != "" and "Ref Spec" not in posdf.columns:
        posdf = add_ref_spec_xml(posdf, posxml)
    if negxml != "" and "Ref Spec" not in negdf.columns:
        negdf = add_ref_spec_xml(negdf, negxml)

    combined_df = outlier_analysis(posdf, negdf, tol=tol)

    # Write Outputs
    combined_df.to_excel("Combined_OutlierAnalysis.xlsx", index=False, sheet_name="Outlier_Analysis")

    namedf = get_unique_names(combined_df)
    # Write to separate sheet in the same Excel file
    with pd.ExcelWriter("Combined_OutlierAnalysis.xlsx", engine='openpyxl', mode='a') as writer:
        namedf.to_excel(writer, sheet_name='Unique_Names', index=False)

def prep_for_rf(df, features, one_hot_encode=False):
    # if "Anomaly_Score" in column, clip it between -5 and 5 for better visualization
    for feature in features:
        if "Anomaly_Score" in feature:
            df[feature] = df[feature].clip(-5, 5)
        if "mahalanobis" in feature:
            df[feature] = df[feature].clip(-10, 10)

    # Encode "Match Grade" based on first letter to integer
    if "Match_Grade" in df.columns and "Match_Grade_Encoded" in features:
        grade_mapping = {'A': 5, 'B': 4, 'C': 3, 'D': 2, 'E': 1}
        df["Match_Grade_Encoded"] = df["Match_Grade"].str[0].map(grade_mapping)
        if one_hot_encode:
            features.remove("Match_Grade_Encoded")
            for grade in grade_mapping.keys():
                colname = f"Match_Grade_{grade}"
                df[colname] = (df["Match_Grade"].str[0] == grade).astype(int)
                features.append(colname)

    # Encode headgroups if not already encoded
    if "Headgroup Encoding" not in df.columns and "Headgroup Encoding" in features:
        df = encode_headgroups(df)
        if one_hot_encode:
            # Remove "Headgroup Encoding" from features and add one hot encoded columns instead
            features.remove("Headgroup Encoding")
            onehotfeatures = df.columns[df.columns.str.startswith("HG_")].tolist()
            features.extend(onehotfeatures)

    # Encode Adducts if not already encoded
    if "Adduct Encoding" not in df.columns and "Adduct Encoding" in features:
        adduct_mapping = {adduct: idx for idx, adduct in enumerate(df["Adduct type"].unique())}
        df["Adduct Encoding"] = df["Adduct type"].map(adduct_mapping)
        if one_hot_encode:
            features.remove("Adduct Encoding")
            for adduct in df["Adduct type"].unique():
                colname = f"Adduct_{adduct}"
                df[colname] = (df["Adduct type"] == adduct).astype(int)
                features.append(colname)

    return df, features

def random_split_rf_data(df, features, target="Low_Quality_Flag", test_size=0.1, random_state=42):
    X = df[features].fillna(0)
    y = df[target]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    return X_train, X_test, y_train, y_test

def dataset_split_rf_data(df, features, target="Low_Quality_Flag", dataset_col="Dataset", test_dataset="Stellar",
                          training_dataset=None, use_percent=1.0, random_state=42):
    if training_dataset is None:
        train_df = df[df[dataset_col] != test_dataset]
        test_df = df[df[dataset_col] == test_dataset]
    elif test_dataset is None:
        train_df = df[df[dataset_col] == training_dataset]
        test_df = df[df[dataset_col] != training_dataset]
    else:
        train_df = df[df[dataset_col] == training_dataset]
        test_df = df[df[dataset_col] == test_dataset]
    if use_percent < 1.0:
        train_df = train_df.sample(frac=use_percent, random_state=random_state)
    X_train = train_df[features].fillna(0)
    y_train = train_df[target]
    X_test = test_df[features].fillna(0)
    y_test = test_df[target]
    return X_train, X_test, y_train, y_test

def random_forest(X_train, X_test, y_train, y_test, features, n_estimators=50, random_state=42, verbose=False):
    clf = RandomForestClassifier(n_estimators=n_estimators, random_state=random_state)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    cr = classification_report(y_test, y_pred, output_dict=True)
    importances = clf.feature_importances_
    feature_importance = pd.Series(importances, index=features).sort_values(ascending=False)
    if verbose:
        print(classification_report(y_test, y_pred))
        print("Feature Importances:")
        print(feature_importance)
    return clf, feature_importance, y_pred, cr

def apply_rf_prediction(df, model=None):
    if model is None:
        modelpath = "RF_Model.pkl"
        if not os.path.exists(modelpath):
            raise ValueError("Model is None. Please provide a trained model.")
        model = joblib.load(modelpath)
    elif type(model) == str:
        if not os.path.exists(model):
            raise ValueError(f"Model file {model} does not exist.")
        model = joblib.load(model)
    features = model.feature_names_in_
    X = df[features].fillna(0)
    y_pred = model.predict(X)
    df["Outlier_Prediction"] = y_pred
    return df

def vip(model):
    t = model.x_scores_
    w = model.x_weights_
    q = model.y_loadings_
    p, h = w.shape
    vips = np.zeros((p,))
    s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
    total_s = np.sum(s)
    for i in range(p):
        weight = np.array([(w[i,j] / np.linalg.norm(w[:,j]))**2 for j in range(h)])
        vips[i] = np.sqrt(p*(s.T @ weight)/total_s)
    return vips

def plsr(X_train, X_test, y_train, y_test, features, n_components=2, verbose=False):
    pls = PLSRegression(n_components=n_components)
    pls.fit(X_train, y_train)
    y_pred = pls.predict(X_test)

    # Feature importance can be approximated by the absolute value of the coefficients
    vip_scores = vip(pls)
    feature_importance = pd.Series(vip_scores, index=features).sort_values(ascending=False)
    # Report
    y_pred_class = (y_pred > 0.5).astype(int)
    cr = classification_report(y_test, y_pred_class, output_dict=True)
    if verbose:
        print(classification_report(y_test, y_pred_class))
        print("Feature Importances:")
        print(feature_importance)

    return pls, feature_importance, y_pred_class, cr

def plot_rf_results(feature_importance, ax=None):
    # Plot the random forest results
    if ax is None:
        plt.figure(figsize=(10, 6))
        ax = plt.gca()
    feature_importance.plot(kind='bar', ax=ax)
    ax.set_ylabel('Importance Score')
    ax.set_xlabel('Features')
    if ax is None:
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":

    topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\HEK\New Libs Tests"
    posfile = "PosIDsRTP2.csv"
    negfile = "NegIDsRTP2.csv"
    posfile = "PosIDsC1_Ref.csv"
    negfile = "NegIDsC1_Ref.csv"
    posxml = ""
    negxml = ""

    os.chdir(topdir)
    outdf = outlier_pipeline(posfile, negfile, posxml, negxml)

