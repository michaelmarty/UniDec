import os
import pandas as pd
from PublicScripts.Lipids.RT_Correction import correct_dataset_rts

topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Shivam\IQX_Stellar_fixed"
os.chdir(topdir)

df = pd.read_csv("Combined_DDA.csv")

# Drop all but Stellar datasets
# df = df[df["Dataset"].isin(["ITIT", "Stellar"])]

correct_dataset_rts(df, "Average Rt(min)", ref_dataset="Stellar")

# Write out corrected dataframe for reference
df.to_csv("Combined_DDA_RT_Corrected.csv", index=False)