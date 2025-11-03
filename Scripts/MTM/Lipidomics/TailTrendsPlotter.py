import os
import LipidFunctions
import pandas as pd
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
import matplotlib

matplotlib.use("WxAgg")

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
                     n_knots=4, degree=2, plot=True):
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
refcol = "Reference RT"

file = r"C1_Ref_Testing.csv"
column_name = "Name"
rt_name = "RT"
class_name = "CompoundClass"
refcol = "RTP"
df = pd.read_csv(file)
df = LipidFunctions.clean_df(df, column=column_name)
df = LipidFunctions.set_tails(df, columnname=column_name)
df = LipidFunctions.encode_headgroups(df, classname=class_name, rtcol_name=rt_name)
df["Residuals"] = np.abs(df[rt_name] - df[refcol])

# Stats on residuals
print(f'Mean Residual: {np.mean(df["Residuals"]):.3f} min')
print(f'Median Residual: {np.median(df["Residuals"]):.3f} min')
print(f'STD of Residuals: {np.std(df["Residuals"]):.3f} min')

print(df[class_name].unique())

# Find columns with HG_ in them
hg_columns = [col for col in df.columns if col.startswith("HG_")]


xc1 = ["Tail Length", "Tail Unsaturation"]
xcols = ["Tail Length", "Tail Unsaturation"] + hg_columns


# Isolate FA ontology

fa_df = df[df[class_name] == "PC"]

print(f'Mean Residual: {np.mean(fa_df["Residuals"]):.3f} min')
print(f'Median Residual: {np.median(fa_df["Residuals"]):.3f} min')
print(f'STD of Residuals: {np.std(fa_df["Residuals"]):.3f} min')

classstd = np.std(fa_df["Residuals"])

# plsr(df, xcols=["Tail Length", "Tail Unsaturation", "Headgroup Encoding"], ycol=rt_name)
# plsr(fa_df, xcols=["Tail Length", "Tail Unsaturation"], ycol=rt_name)
# plsr(df , xcols=xcols, ycol=rt_name)
# mars(fa_df, xcols = xc1, ycol=rt_name)
# mars(df, xcols=xcols, ycol=rt_name)

# random_forest(fa_df, xcols=["Tail Length", "Tail Unsaturation"], ycol=rt_name)
# random_forest(df, xcols=["Tail Length", "Tail Unsaturation", "Headgroup Encoding"], ycol=rt_name)

# # Create 3D Plot of Tail Length, Tail Unsaturation vs RT
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(121, projection='3d')
# ax.set_title('3D Plot of Tail Length, Tail Unsaturation vs RT')
# ax.set_xlabel('Tail Length')
# ax.set_ylabel('Tail Unsaturation')
# ax.set_zlabel('Retention Time (min)')
#
#
#
# # Fit 3D data using polynomial surface fitting
# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.linear_model import LinearRegression
# X = fa_df[["Tail Length", "Tail Unsaturation"]].values
# Y = fa_df[rt_name].values
# poly = PolynomialFeatures(degree=3)
# X_poly = poly.fit_transform(X)
# model = LinearRegression()
# model.fit(X_poly, Y)
# # Create grid to plot surface
# xlen = np.linspace(fa_df["Tail Length"].min(), fa_df["Tail Length"].max(), 10)
# xunsat = np.linspace(fa_df["Tail Unsaturation"].min(), fa_df["Tail Unsaturation"].max(), 10)
# xlen, xunsat = np.meshgrid(xlen, xunsat)
# xgrid = np.column_stack((xlen.ravel(), xunsat.ravel()))
# xgrid_poly = poly.transform(xgrid)
# ygrid = model.predict(xgrid_poly)
# ygrid = ygrid.reshape(xlen.shape)
# ax.plot_surface(xlen, xunsat, ygrid, alpha=0.5, color='r')
#
#
#
# residuals = np.abs(Y - model.predict(X_poly))
#
# # Color by residuals
# colors = plt.cm.viridis((residuals - residuals.min()) / (residuals.max() - residuals.min()))
# ax.scatter(fa_df["Tail Length"], fa_df["Tail Unsaturation"], fa_df[rt_name], marker='o', c=colors)
#
# # Plot residuals
# ax = fig.add_subplot(122)
# ax.set_title('Residuals of 3D Fit')
# ax.set_xlabel('Index')
# ax.set_ylabel('Residual (min)')
# ax.scatter(range(len(residuals)), residuals, marker='o', c=colors)
#
#
# plt.figure()

unsats = np.unique(fa_df["Tail Unsaturation"].values)
colors = plt.cm.viridis(np.linspace(0, 1, len(unsats)))

for i, u in enumerate(unsats):
    subdf = fa_df[fa_df["Tail Unsaturation"] == u]
    # Sort subdf by tail lengths
    subdf = subdf.sort_values(by=["Tail Length"])
    tail_lengths = subdf["Tail Length"].values
    rts = subdf[rt_name].values
    refrts = subdf[refcol].values


    # print(subdf[[column_name, "Tail Length", rt_name, refcol]])
    # areas = subdf["Total score"].values
    # residualcolors = plt.cm.viridis(subdf["Residuals"].values / np.max(subdf["Residuals"].values))


    # plt.subplot(122)
    plt.plot(tail_lengths, refrts, marker='x', color=colors[i])

    #
    # plt.subplot(121)
    plt.plot(tail_lengths, rts, marker='o', linestyle="--", label=f'Unsaturation {u}', color=colors[i])

    # Mark outliers
    for j in range(len(rts)):
        if np.abs(rts[j] - refrts[j]) > 0.25:
            plt.plot(tail_lengths[j], rts[j], marker='s', markersize=10, markerfacecolor='none', markeredgecolor='r', markeredgewidth=2)

    # if len(tail_lengths) < 6:
    #     continue
    # p = np.polyfit(tail_lengths, rts, 3)#, w=areas)
    # p = np.poly1d(p)
    # xnew = np.linspace(np.min(tail_lengths), np.max(tail_lengths), num=100)
    # ynew = p(xnew)
    # resid = rts - p(tail_lengths)
    # plt.plot(tail_lengths, rts, marker='o', label=f'Unsaturation {u}', color=colors[i])
    #
    # plt.plot(xnew, ynew, linestyle='--', color=colors[int(u)])
plt.xlabel('Tail Length')
plt.ylabel('Retention Time (min)')
plt.legend()
plt.show()