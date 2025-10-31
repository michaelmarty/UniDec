import os
import LipidFunctions
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_predict
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
# python
from sklearn.preprocessing import SplineTransformer
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
import numpy as np
import matplotlib.pyplot as plt

def plsr(df, xcols, ycol):

    X = df[xcols].values
    Y = df[ycol].values

    pls = PLSRegression(n_components=2)
    pls.fit(X, Y)
    # Y_pred = cross_val_predict(pls, X, Y, cv=5)
    Y_pred = pls.predict(X)

    r2 = r2_score(Y, Y_pred)
    rmse = np.sqrt(mean_squared_error(Y, Y_pred))

    print(f'R^2: {r2:.3f}, RMSE: {rmse:.3f}')

    plt.scatter(Y, Y_pred)
    plt.xlabel('Measured')
    plt.ylabel('Predicted')
    plt.title('PLSR Prediction')
    plt.plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--')
    plt.show()

def mars(df, xcols, ycol, test_size=0.2, random_state=42,
                     n_knots=10, degree=3, plot=True):
    """
    Spline-based alternative to MARS using scikit-learn's SplineTransformer.
    Returns: (model, r2, rmse)
    """
    X = df[xcols].values
    y = df[ycol].values

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
    X_train, X_test, y_train, y_test = X, X, y, y

    pipeline = make_pipeline(
        SplineTransformer(n_knots=n_knots, degree=degree, include_bias=False),
        LinearRegression()
    )

    pipeline.fit(X_train, y_train)
    y_pred = pipeline.predict(X_test)

    r2 = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    # Print outliers where y_test - y_pred > 0.5
    # outliers = np.where(np.abs(y_test - y_pred) > 0.5)[0]
    # if len(outliers) > 0:
    #     print("Outliers:")
    #     for idx in outliers:
    #         print(f"Measured: {y_test[idx]}, Predicted: {y_pred[idx]}, Difference: {y_test[idx] - y_pred[idx]}")

    if plot:
        plt.figure()
        plt.scatter(y_test, y_pred, alpha=0.7)
        mn, mx = np.min(y_test), np.max(y_test)
        plt.plot([mn, mx], [mn, mx], "r--")
        plt.xlabel("Measured")
        plt.ylabel("Predicted")
        plt.title(f"Spline Regression (R^2={r2:.3f}, RMSE={rmse:.3f})")
        plt.show()

    return pipeline, r2, rmse

os.chdir(r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults")

file = r"PosIDsC1.csv"
column_name = "Metabolite name"
rt_name = "Average Rt(min)"
class_name = "Ontology"
df = pd.read_csv(file)
df = LipidFunctions.clean_df(df, column=column_name)
df = LipidFunctions.set_tails(df, columnname=column_name)
df = LipidFunctions.encode_headgroups(df, classname=class_name, rtcol_name=rt_name)
print(df.columns)

print(df[class_name].unique())

# Find columns with HG_ in them
hg_columns = [col for col in df.columns if col.startswith("HG_")]
print("Headgroup encoding columns:", hg_columns)

xc1 = ["Tail Length", "Tail Unsaturation"]
xcols = ["Tail Length", "Tail Unsaturation"] + hg_columns
print("X columns:", xcols)

# Isolate FA ontology

fa_df = df[df[class_name] == "TG"]

# plsr(df, xcols=["Tail Length", "Tail Unsaturation", "Headgroup Encoding"], ycol=rt_name)
# plsr(fa_df, xcols=["Tail Length", "Tail Unsaturation"], ycol=rt_name)
# plsr(df , xcols=xcols, ycol=rt_name)
# mars(fa_df, xcols = xc1, ycol=rt_name)
# mars(df, xcols=xcols, ycol=rt_name)

# random_forest(fa_df, xcols=["Tail Length", "Tail Unsaturation"], ycol=rt_name)
#random_forest(df, xcols=["Tail Length", "Tail Unsaturation", "Headgroup Encoding"], ycol=rt_name)



unsats = np.unique(fa_df["Tail Unsaturation"].values)

for u in unsats:
    subdf = fa_df[fa_df["Tail Unsaturation"] == u]
    tail_lengths = subdf["Tail Length"].values
    rts = subdf[rt_name].values
    areas = subdf["Total score"].values

    # Perform cubic interpolation for smoother curve
    sorted_indices = np.argsort(tail_lengths)
    tail_lengths = tail_lengths[sorted_indices]
    rts = rts[sorted_indices]

    if len(tail_lengths) < 4:
        plt.plot(tail_lengths, rts, marker='o', label=f'Unsaturation {u}')
        continue

    # # If any tail lengths are duplicated, average their rts
    # unique_lengths, indices = np.unique(tail_lengths, return_index=True)
    # averaged_rts = []
    # for length in unique_lengths:
    #     averaged_rts.append(np.average(rts[tail_lengths == length], weights = areas[tail_lengths == length]))
    # tail_lengths = unique_lengths
    # rts = np.array(averaged_rts)
    #
    # if len(tail_lengths) < 4:
    #     p=interp.interp1d(tail_lengths, rts, kind='linear', fill_value="extrapolate")
    # else:
    #     p = interp.interp1d(tail_lengths, rts, kind='cubic', fill_value="extrapolate")
    p = np.polyfit(tail_lengths, rts, 3, w=areas)
    p = np.poly1d(p)
    xnew = np.linspace(np.min(tail_lengths), np.max(tail_lengths), num=100)
    ynew = p(xnew)

    outliers = np.where(np.abs(rts - p(tail_lengths)) > 0.25)[0]
    if len(outliers) > 0:
        print(f"Outliers for unsaturation {u}:")
        for idx in outliers:
            print(f"Tail Length: {tail_lengths[idx]}, Rt: {rts[idx]}, Predicted Rt: {p(tail_lengths[idx])}")
            print("Name:", subdf.iloc[idx][column_name])

    plt.plot(tail_lengths, rts, marker='o', label=f'Unsaturation {u}')
    plt.plot(xnew, ynew, linestyle='--')
plt.xlabel('Tail Length')
plt.ylabel('Retention Time (min)')
plt.show()