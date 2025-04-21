import os
import pandas as pd
import glob
import numpy as np

from openpyxl import workbook

temps = [15, 20, 25, 30, 35]

globaldir = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_W14A Titrations All\\WT_W14A_"
fname = "Extract_grid_2D_Extract6.txt"
for t in temps:
    temp = t
    topdir = globaldir+str(t)+"C"
    print(topdir)
    files = []

    for path in glob.glob(f'{topdir}/*/**/', recursive=True):
        testpath = os.path.join(path,fname)
        if os.path.isfile(testpath):
            files.append(testpath)
    print(files)


combined_data = pd.DataFrame()
df = pd.read_csv(testpath, delimiter=' ')
print(df)

combined_data = combined_data.append(df, ignore_index=True)
exit()

excel_file_path = os.path.join("Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_W14A Titrations All", "2DCombiner.xlsx")

excel_write = pd.ExcelWriter(excel_file_path, engine='openpyxl')

combined_data.to_excel(excel_write, sheet_name='2DCombiner', index=False)

excel_write.save()
