from pathlib import Path
import sys

# Make the package available when running this script directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import os
import pandas as pd
from PublicScripts.Lipids.RT_Correction import correct_dataset_rts


topdir = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\Final Results\DDA Pipeline"
os.chdir(topdir)

df = pd.read_csv("Combined_DDA.csv")

# Drop all but Stellar datasets
# df = df[df["Dataset"].isin(["ITIT", "Stellar"])]

df = correct_dataset_rts(df, "Average Rt(min)", ref_dataset="Stellar")

# Write out corrected dataframe for reference
df.to_csv("Combined_DDA_RT_Corrected.csv", index=False)