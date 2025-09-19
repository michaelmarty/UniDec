import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from matplotlib import rcParams
import os

#rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5
#rcParams['font.size']=8

dpath = "Z:\Group Share\\Melanie Odenkirk\\20220815_Ecoli_d7_lowerCone\\Exported_TLs"
os.chdir(dpath)

df = pd.read_csv(dpath +"\\220816_MS1_Ecoli_lipid_ND_res_outputs.csv")

filter_col1 = [col for col in df if col.endswith('Area')]
df_area= df[filter_col1]


extracted_col1 = df["Molecule Name"]
df_area = df_area.join(extracted_col1)
extracted_col2 = df["Molecule List Name"]
df_area = df_area.join(extracted_col2)
extracted_col3 = df["Precursor Adduct"]
df_area = df_area.join(extracted_col3)
df_area = df_area.set_index('Molecule Name')

print(df_area)