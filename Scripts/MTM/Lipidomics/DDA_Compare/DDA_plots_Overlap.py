import os
from PublicScripts.Lipids.LipidFunctions import *
from venn import venn
import matplotlib.pyplot as plt
import matplotlib as mpl
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *

mpl.use("WxAgg")

topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Annika\IQ-X\Feb 2026 UntDDA\MSDial\Marked Low Score\RTP2"
topdir = r"Z:\Group Share\Annika\IQX vs Stellar DDA\RTP2 Peak Masters and xml"
topdir = r"Z:\Group Share\Shivam\Final_results"
topdir = r"Z:\Group Share\Shivam\IQX Stellar 3\FinalResults\AS"
os.chdir(topdir)

df = pd.read_excel("Combined_DDA_Full.xlsx", sheet_name="Full Data")
df = pd.read_excel("Combined_DDA_OutlierAnalysis_Results.xlsx", sheet_name="All Data")
# df = df[df["Outlier_Prediction"] == 0]

# # Drop any with class "BRSE"
# df = df[df["Ontology"] != "BRSE"]
# df = df[df["Ontology"] != "HBMP"]

datasets = ["OTOT", "OTIT", "ITIT", "Stellar"]
datasets = ["OTOT", "OTIT",  "ITIT","Stellar", "Stellar_both"]
compare_classes_plot(df, datasets=datasets, use_simple_names=False, drop_low_quality=True)


