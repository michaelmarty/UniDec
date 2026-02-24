from PublicScripts.Lipids.LipidFunctions import *
import os

topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Shivam\IQX_Stellar_fixed"
os.chdir(topdir)

# df = pd.read_csv("Combined_DDA_RT_Corrected_FullNames.csv")
df = pd.read_excel("Combined_DDA_OutlierAnalysis_Results.xlsx", sheet_name="All Data")
# Drop all where Outlier_Result is 1
df = df[df["Outlier_Prediction"] == 0]

# Select only Dataset of "OTOT"
df = df[df["Dataset"] == "OTIT"]
print(len(df))
mzrt_plot(df)