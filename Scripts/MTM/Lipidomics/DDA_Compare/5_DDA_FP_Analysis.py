import os
import pandas as pd
import PublicScripts.Lipids.LipidOutlierAnalysis as ot
from PublicScripts.Lipids.LipidFunctions import smiles_analysis


topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

filepath = r"Combined_DDA_Full.xlsx"
df = pd.read_excel(filepath)

# Keep only rows labeled True
df = df[df["Overall Match Status"] == True].copy()

df = smiles_analysis(df)

posdf = df[df["Polarity"] == "Positive"]
negdf = df[df["Polarity"] == "Negative"]

outdf = ot.outlier_analysis(
    posdf,
    negdf,
    tol=0.3,
    contribs=False,
    add_all_cols_heads=True,
    add_all_cols_tails=False,
    drop_cols=False,
    do_tails=False
)

outdf.to_excel("Combined_DDA_OutlierAnalysis.xlsx", index=False)