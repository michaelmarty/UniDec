import os
import pandas as pd
from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *


topdir = r"Z:\Group Share\Shivam\DDA"
topdir = r"Z:\Group Share\Shivam\IQX Stellar 3\FinalResults"
os.chdir(topdir)

df = pd.read_csv("Combined_DDA_RT_Corrected.csv")

df = full_name_assignment(df, rttol=0.2, rtcol="Average Rt(min)")
df = get_number_of_datasets(df, namecol="Metabolite name")
write_dfs_to_excel(df)
