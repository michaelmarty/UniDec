from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *
import pandas as pd
from PublicScripts.Lipids.RT_Correction import correct_dataset_rts

topdir = r"Z:\Group Share\Annika\DDA_Compare"
os.chdir(topdir)

datasets = ["Stellar", "Stellar_Both"]

if False:
    df = prep_datasets(topdir, datasets, dtype="mzmine") # 1
    df = correct_dataset_rts(df, "Average Rt(min)", ref_dataset="Stellar", plot=False) # 2
    df = dataset_split(df) # 3
    df = outlier_setup_dda(df) # 5




mpl.use("WxAgg")


df = pd.read_excel("Combined_DDA_OutlierAnalysis.xlsx")
compare_classes_plot(df, datasets=datasets, use_simple_names=False, drop_low_quality=True)

