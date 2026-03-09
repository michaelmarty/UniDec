from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *
import pandas as pd
from PublicScripts.Lipids.RT_Correction import correct_dataset_rts

topdir = r"Z:\Group Share\Shivam\IQX Stellar 3\FinalResults"
os.chdir(topdir)

datasets = ["ITIT", "OTIT", "OTOT", "Stellar"]

df = prep_datasets(topdir, datasets) # 1
df = correct_dataset_rts(df, "Average Rt(min)", ref_dataset="Stellar", plot=False) # 2
df = dataset_split(df) # 3
df = outlier_setup_dda(df) # 5


