import csv
import os
import pandas
from matplotlib import pyplot as plt
import numpy as np


dpath = "Z:\Group Share\Melanie Odenkirk\\2022_Synapt_Experiments_Skyline_Processing\\Ecoli_transition_lists"
os.chdir(dpath)

df = pandas.read_csv(dpath +"\\220720_DIADDA_Skyline_Results.csv")

df_MHadduct = df.loc[df['Precursor Adduct'] == '[M-H]']

df_filter = df_MHadduct.drop(columns=['Molecule Formula', 'Precursor Ion Formula','Precursor Adduct', 'Precursor Neutral Formula',
                 'Precursor Charge', 'Explicit Retention Time','Explicit Retention Time Window', 'Product Mz',
                 'Product Charge', 'Product Ion Formula','Product Neutral Formula', 'Product Adduct',
                 'Fragment Ion Type', 'Total Area','Isotope Dist Index', 'Isotope Dist Rank',
                 'Isotope Dist Proportion', 'Full Scan Filter Width','Product Decoy Mz Shift'])

df_precursor = df_filter.loc[df_filter['Fragment Ion'] == 'precursor']

df_precursor["Lipid_ID_class"] = df_precursor['Molecule List Name'].astype(str) +"_"+ df["Molecule Name"]
df_index = df_precursor.set_index('Molecule Name')

df_MS1 = df_index.loc[:, df_index.columns.str.contains("MS1")]
df_TIC = df_index.loc[:, df_index.columns.str.contains("Total Ion Current")]
df_Area = df_index.loc[:, ~df_index.columns.str.contains("Total")]
df_Area = df_Area.drop(columns=['Molecule List Name', 'Lipid_ID_class', 'Precursor Mz',
                                'Fragment Ion'])
print(df_MS1)
df_log2MS1 = np.log2(df_MS1)

df_sort= df_MS1.sort_values(by=['20220713_Neg_DIA_LipidNDwithmasterISD_newbuffer_Quad440_run1 Total Area MS1'])
df_sortlog= df_log2MS1.sort_values(by=['20220713_Neg_DIA_LipidNDwithmasterISD_newbuffer_Quad440_run1 Total Area MS1'])
df_sortlog.plot(y=["20220713_Neg_DIA_LipidNDwithmasterISD_newbuffer_Quad440_run1 Total Area MS1"], kind="bar", figsize=(9, 8))

plt.show()