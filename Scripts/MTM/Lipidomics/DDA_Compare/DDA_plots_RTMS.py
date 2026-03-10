from PublicScripts.Lipids.LipidFunctions import *
import os

topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Annika\IQ-X\Feb 2026 UntDDA\MSDial\Marked Low Score\RTP2"
topdir = r"Z:\Group Share\Annika\IQX vs Stellar DDA\RTP2 Peak Masters and xml"
topdir = r"Z:\Group Share\Shivam\IQX Stellar 3\FinalResults"
os.chdir(topdir)

# df = pd.read_csv("Combined_DDA_RT_Corrected_FullNames.csv")
df = pd.read_excel("Combined_DDA_OutlierAnalysis_Results.xlsx", sheet_name="All Data")
# Drop all where Outlier_Result is 1
df = df[df["Outlier_Prediction"] == 0]

# Select only Dataset of "OTOT"
df = df[df["Dataset"] == "OTIT"]
print(len(df))
mzrt_plot(df, write_file=True)