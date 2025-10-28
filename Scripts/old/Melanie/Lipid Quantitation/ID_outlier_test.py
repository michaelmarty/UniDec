
#########________________IMPORT Libraries_________############
import csv
import os
from typing import List, Any
import pandas
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt
import numpy

#########________________Read in and extract data_________############
dpath = "Z:\Group Share\\Melanie Odenkirk\\20220815_Ecoli_d7_lowerCone\\Exported_TLs"
os.chdir(dpath)

df = pandas.read_csv(dpath +"\\220816_MS1_Ecoli_lipid_ND_res_outputs_1.csv")



filter_col1 = [col for col in df if col.endswith('Retention Time')]
df_rt= df[filter_col1]
df_rt['mean_rt'] = df_rt.mean(axis=1)
df_rt = df_rt[['mean_rt']]


filter_col = [col for col in df if col.endswith('Mass Error PPM')]
df_ppm = df[filter_col]
df_ppm['mean_ppm'] = df_ppm.mean(axis=1)
df_ppm = df_ppm[['mean_ppm']]

#filter_col = [col for col in df if col.endswith('Collision Cross Section')]
#df_ccs = df[filter_col]
#df_ccs['mean_ccs'] = df_ccs.mean(axis=1)
#df_ccs = df_ccs[['mean_ccs']]

df_outlier = pd.merge(df_rt, df_ppm, left_index=True, right_index=True)

extracted_col1 = df["Molecule Name"]
df_outlier = df_outlier.join(extracted_col1)
extracted_col2 = df["Molecule List Name"]
df_outlier = df_outlier.join(extracted_col2)
df_stdev_limits = df_outlier.groupby("Molecule List Name").apply(numpy.std)
df_class_average = df_outlier.groupby("Molecule List Name").apply(numpy.mean)
df_outlier = df_outlier.set_index('Molecule Name')

df_limit = pandas.merge(df_stdev_limits, df_class_average, on = "Molecule List Name")
df_limit = pandas.merge(df_outlier,df_limit, left_on ="Molecule List Name", right_on = "Molecule List Name")
df_limit = df_limit.join(extracted_col1)
df_limit = df_limit.set_index('Molecule Name')

with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(df_limit)

df_limit.loc[df_limit['mean_rt'] > abs((df_limit['mean_rt_y']+abs(2*df_limit['mean_rt_x']))), 'RT_Outlier'] = 'True'
df_limit.loc[df_limit['mean_rt'] <= abs((df_limit['mean_rt_y']+abs(2*df_limit['mean_rt_x']))), 'RT_Outlier'] = 'False'

df_limit.loc[df_limit['mean_ppm'] > abs((df_limit['mean_ppm_y']+abs(2*df_limit['mean_ppm_x']))), 'PPM_Outlier'] = 'True'
df_limit.loc[df_limit['mean_ppm'] <= abs((df_limit['mean_ppm_y']+abs(2*df_limit['mean_ppm_x']))), 'PPM_Outlier'] = 'False'

print(df_limit.loc[df_limit['RT_Outlier'] == 'True'])
print(df_limit.loc[df_limit['PPM_Outlier'] == 'True'])

from pathlib import Path
filepath = Path('Z:/Group Share//Melanie Odenkirk//20220815_Ecoli_d7_lowerCone//Exported_TLs/names_test.csv')
df_limit.to_csv(filepath)


#########________________Plot data_________############
fig, axs = plt.subplots(ncols=2)

sb.boxplot(x = 'Molecule List Name',y = 'mean_rt', data = df_limit, ax=axs[0])
sb.set_palette("Set2")
sb.stripplot(x = "Molecule List Name",
              y = "mean_rt",
              hue = 'RT_Outlier',
              palette = ["k", "r"],
              data = df_limit, ax=axs[0])
sb.boxplot(x = 'Molecule List Name',y = 'mean_ppm', data = df_limit, ax=axs[1])
sb.set_palette("Set2")
sb.stripplot(x = "Molecule List Name",
              y = "mean_ppm",
             hue = 'PPM_Outlier',
             palette = ["k", "r"],
              data = df_limit, ax=axs[1])
plt.show()


#########________________Read in and extract data_________############
