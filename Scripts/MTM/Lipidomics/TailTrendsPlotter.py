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

def random_forest(df, xcols, ycol):
    X = df[xcols].values
    Y = df[ycol].values

    rf = RandomForestRegressor(n_estimators=100, random_state=42)
    rf.fit(X, Y)
    Y_pred = rf.predict(X)

    # X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
    #
    # rf = RandomForestRegressor(n_estimators=100, random_state=42)
    # rf.fit(X_train, Y_train)
    #
    # Y_pred = rf.predict(X_test)
    #
    # r2 = r2_score(Y_test, Y_pred)
    # rmse = np.sqrt(mean_squared_error(Y_test, Y_pred))
    #
    # print(f'Random Forest R^2: {r2:.3f}, RMSE: {rmse:.3f}')
    #

    diffs = Y - Y_pred
    # Count how many are within 0.5 minutes
    within_05 = np.sum(np.abs(diffs) <= 0.5)
    total = len(Y)
    print(f'Random Forest: {within_05}/{total} predictions within 0.5 minutes ({(within_05/total)*100:.2f}%)')
    # Count how many are outside 0.5 minutes
    outside_05 = np.sum(np.abs(diffs) > 0.5)
    print(f'Random Forest: {outside_05}/{total} predictions outside 0.5 minutes ({(outside_05/total)*100:.2f}%)')


    plt.scatter(Y, Y_pred)
    plt.xlabel('Measured')
    plt.ylabel('Predicted')
    plt.title('Random Forest Prediction')
    plt.plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--')
    plt.show()

file = r"C:\Data\IQ-X\MartyLab_AX\PosIDs.csv"
column_name = "Metabolite name"
rt_name = "Average Rt(min)"
class_name = "Ontology"
df = pd.read_csv(file)
df = LipidFunctions.clean_df(df, column=column_name)
df = LipidFunctions.set_tails(df, columnname=column_name)
df = LipidFunctions.encode_headgroups(df, classname=class_name, rtcol_name=rt_name)
print(df.columns)

print(df[class_name].unique())


# Isolate FA ontology

fa_df = df[df[class_name] == "DG"]

# plsr(df, xcols=["Tail Length", "Tail Unsaturation", "Headgroup Encoding"], ycol=rt_name)
#plsr(fa_df, xcols=["Tail Length", "Tail Unsaturation"], ycol=rt_name)

random_forest(fa_df, xcols=["Tail Length", "Tail Unsaturation"], ycol=rt_name)
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

    outliers = np.where(np.abs(rts - p(tail_lengths)) > 0.5)[0]
    if len(outliers) > 0:
        print(f"Outliers for unsaturation {u}:")
        for idx in outliers:
            print(f"Tail Length: {tail_lengths[idx]}, Rt: {rts[idx]}, Predicted Rt: {p(tail_lengths[idx])}")

    plt.plot(tail_lengths, rts, marker='o', label=f'Unsaturation {u}')
    plt.plot(xnew, ynew, linestyle='--')
plt.xlabel('Tail Length')
plt.ylabel('Retention Time (min)')
plt.show()