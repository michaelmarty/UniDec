import os
import pandas as pd
from PublicScripts.Lipids.LipidOutlierAnalysis import apply_rf_prediction, basic_class_network_analysis
from PublicScripts.Lipids.LipidFunctions import msdial_to_transitionlist, tl_to_il, derplicate_il

topdir = r"Z:\Group Share\Shivam\DDA"
os.chdir(topdir)

filepath = r"Combined_DDA_Full.xlsx"
df = pd.read_excel(filepath, sheet_name="AX Total")
print(len(df))

df = basic_class_network_analysis(df)
print(df.columns)
df = apply_rf_prediction(df)

# Print number with Outerlier Prediction of "Good"
good_count = len(df[df["Outlier_Prediction"] == 0])
print(f"Number of good predictions: {good_count} out of {len(df)}")

dfgood = df[df["Outlier_Prediction"] == 0]
tl_df = msdial_to_transitionlist(dfgood)
ildf = tl_to_il(tl_df)
ildf = derplicate_il(ildf)

# Write out transition list and isolation list to separate sheets in the same Excel file
with pd.ExcelWriter("AqX_Predictions.xlsx") as writer:
    df.to_excel(writer, sheet_name="Predictions", index=False)
    tl_df.to_excel(writer, sheet_name="Transition List", index=False)
    ildf.to_excel(writer, sheet_name="Isolation List", index=False)
