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

# matplotlib.use("WxAgg")

os.chdir(r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults")

file = r"PosIDsC1.csv"
column_name = "Metabolite name"
rt_name = "Average Rt(min)"
class_name = "Ontology"
refcol = "Reference RT"

os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

file = r"C1_Ref_Testing_HighW10.csv"
column_name = "Name"
rt_name = "RT"
class_name = "CompoundClass"
refcol = "RTP"
# refcol= "Reference RT"

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

# Isolate FA ontology
fa_df = df[df[class_name] == "PC"]

print(f'Mean Residual: {np.mean(fa_df["Residuals"]):.3f} min')
print(f'Median Residual: {np.median(fa_df["Residuals"]):.3f} min')
print(f'STD of Residuals: {np.std(fa_df["Residuals"]):.3f} min')

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
    plt.plot(tail_lengths, rts, marker='o', linestyle="--", label=f'Unsaturation {u}', color=colors[i])

    # Mark outliers
    for j in range(len(rts)):
        if np.abs(rts[j] - refrts[j]) > 0.25:
            plt.plot(tail_lengths[j], rts[j], marker='s', markersize=10, markerfacecolor='none', markeredgecolor='r', markeredgewidth=2)


plt.xlabel('Tail Length')
plt.ylabel('Retention Time (min)')
plt.legend()
plt.show()