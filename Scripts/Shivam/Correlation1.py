import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from PublicScripts.Lipids.LipidFunctions import encode_headgroups
# from sklearn.naive_bayes import GaussianNB # Not as good
# Import neural net classifier from sklearn.neural_network
# from sklearn.neural_network import MLPClassifier # Not as good as random forest
# from sklearn.inspection import permutation_importance
# from sklearn.svm import SVC
# from sklearn.inspection import permutation_importance

if __name__ == "__main__":
    os.chdir(r"Z:\Group Share\Shivam\Reducing_FP\Tagged_Peaks")
    filepath = r"Combined_OutlierAnalysis.xlsx"

    features = ['RTdiff', 'LogSN', 'Headgroup_Match', 'Tail_Match_Fraction',
                "Class_z_Anomaly_Score", "Tail_z_Anomaly_Score","mahalanobis", "Simple dot product",
                "Weighted dot product", "Reverse dot product"]
    # 'Total score'

    df = pd.read_excel(filepath)
    # df = encode_headgroups(df)
    # features = features + ["Headgroup Encoding"]

    # if "Anomaly_Score" in column, clip it between -5 and 5 for better visualization
    for feature in features:
        if "Anomaly_Score" in feature:
            df[feature] = df[feature].clip(-5, 5)
        if "mahalanobis" in feature:
            df[feature] = df[feature].clip(-10, 10)

    # Encode "Match Grade" based on first letter to integer
    if "Match_Grade" in df.columns:
        grade_mapping = {'A': 5, 'B': 4, 'C': 3, 'D': 2, 'E': 1}
        df["Match_Grade_Encoded"] = df["Match_Grade"].str[0].map(grade_mapping)
        features.append("Match_Grade_Encoded")

    # Set up a random forest to predict Low_Quality_Flag

    X = df[features].fillna(0)
    y = df["Low_Quality_Flag"]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)
    clf = RandomForestClassifier(n_estimators=50, random_state=42)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    print(classification_report(y_test, y_pred))
    importances = clf.feature_importances_
    feature_importance = pd.Series(importances, index=features).sort_values(ascending=False)
    print("Feature Importances:")
    print(feature_importance)

    # Plot the random forest results
    plt.figure(figsize=(10, 6))
    feature_importance.plot(kind='bar')
    plt.title('Feature Importances from Random Forest Classifier')
    plt.ylabel('Importance Score')
    plt.xlabel('Features')
    plt.tight_layout()
    plt.show()


    # Plot features with whether Low_Quality_Flag is true or false
    for feature in features:
        plt.figure(figsize=(8, 6))
        for flag in [True, False]:
            subset = df[df["Low_Quality_Flag"] == flag]
            data = subset[feature].dropna()

            plt.hist(data, bins=50, alpha=0.5, label=f'Low_Quality_Flag={flag}')
        plt.title(f'Distribution of {feature} by Low_Quality_Flag')
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
    plt.show()