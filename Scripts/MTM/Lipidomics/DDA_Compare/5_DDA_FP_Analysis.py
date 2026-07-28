import os
import pandas as pd
import PublicScripts.Lipids.LipidOutlierAnalysis as ot
from PublicScripts.Lipids.LipidFunctions import (smiles_analysis, class_network_analysis)


topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

filepath = r"Combined_Fragment_Assignment_Colored.xlsx"
df = pd.read_excel(filepath)

# Keep only rows labeled True
df = df[df["Overall Match Status"] == True].copy()

df = smiles_analysis(df)

posdf = df[df["Polarity"] == "Positive"]
negdf = df[df["Polarity"] == "Negative"]

if negdf.empty:
    negdf = posdf.iloc[0:0].copy()

if negdf.empty:
    print("Only positive-mode data found. Running positive-mode analysis only.")

    outdf, _ = class_network_analysis(
        posdf,
        tol=0.3,
        add_all_cols=True
    )

elif posdf.empty:
    print("Only negative-mode data found. Running negative-mode analysis only.")

    outdf, _ = class_network_analysis(
        negdf,
        tol=0.3,
        add_all_cols=True
    )

else:
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