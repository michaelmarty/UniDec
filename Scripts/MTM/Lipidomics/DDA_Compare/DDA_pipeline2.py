from PublicScripts.Lipids.Lipid_DDA_Dataset_Compare import *
import pandas as pd
from PublicScripts.Lipids.RT_Correction import correct_dataset_rts

mpl.use("WxAgg")

topdir = r"Z:\Group Share\Annika\DDA_Compare\mzmine"
os.chdir(topdir)

datasets = ["Stellar", "OTOT", "OTIT", "ITIT", "IQX", "Apex"]

if True:
    df = prep_datasets(topdir, datasets, dtype="mzmine", rttol=0.3) # 1
    df = correct_dataset_rts(df, "Average Rt(min)", ref_dataset="Stellar", plot=False) # 2
    df = correct_dataset_rts(df, "Reference RT", ref_dataset="Stellar", plot=False)  # 2
    df = dataset_split(df, rttol=0.4) # 3
    df = outlier_setup_dda(df, do_heads=False) # 5

df = pd.read_excel("Combined_DDA_OutlierAnalysis.xlsx")
compare_classes_plot(df, datasets=datasets, use_simple_names=False, drop_low_quality=True)

