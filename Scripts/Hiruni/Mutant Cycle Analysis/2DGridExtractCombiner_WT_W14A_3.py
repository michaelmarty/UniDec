import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rcParams
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import scipy.stats as stats
from matplotlib import rcParams
import statsmodels.api as sm
import scipy.optimize as opt
import scipy.stats as stats
import pandas as pd
import os
import pandas as pd
import glob
import numpy as np

from openpyxl import workbook
import os
import pandas as pd

import os
import pandas as pd
from openpyxl import load_workbook


root_folder = "Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final8\\Temperature Data"

dfs = []
combined_data = pd.DataFrame
file = []
name = []

for temp_folder in os.listdir(root_folder):
    temp_path = os.path.join(root_folder, temp_folder)

    if os.path.isdir(root_folder):
        triplicate_dfs = []

        triplicate_df = pd.read_csv(file_path, delimiter=" ")
        triplicate_dfs.append(triplicate_df)

        combined_df = pd.concat(triplicate_dfs)
        combined_df['Temperature'] = temp_folder
        dfs.append(combined_df)

        final_df = pd.concat(dfs, ignore_index=True)
        excel_file_path = 'combined_Data.xlsx'
        final_df.to_excel(excel_file_path, index=False)
        final_df.to_Excel(excel_file_path, index=False)
        print(f'Combined data saved tp {excel_file_path}')