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
os.chdir(topdir)

df = pd.read_excel("Combined_DDA_Full.xlsx", sheet_name="Full Data")

compare_classes_plot(df)


