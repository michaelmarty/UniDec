import pandas as pd
import numpy as np
import os
from Scripts.LipiDec.AddFeaturestoDF import get_overall_tails, get_class

def fix_name(name):
    name = name.replace("("," ")
    name = name.replace(")","")
    return name


mergeResult_path = r"D:\MS_Data\2025-04-29_OE240_Stellar_UTAustin\SearchResults\mergeResult.csv"
os.chdir(os.path.dirname(mergeResult_path))

mergeResult_df = pd.read_csv(mergeResult_path)
names = mergeResult_df["LipidMolec"]

for name in names:
    fixed_name = fix_name(name)
    print(fixed_name)
    print(get_class(fixed_name))