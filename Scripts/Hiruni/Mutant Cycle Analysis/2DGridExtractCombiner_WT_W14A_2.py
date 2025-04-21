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

temps = [15, 20, 25, 30, 35]

globaldir = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\WT_W14A Titrations All\\WT_W14A_"
fname = "Extract_grid_2D_Extract6.txt"

for t in temps:
    temp = t
    topdir = globaldir + str(t) + "C"
    print(topdir)

    files = []
    data = []
    df = []

    for path in glob.glob(f'{topdir}/*/**/', recursive=True):
        testpath = os.path.join(path, fname)
        if os.path.isfile(testpath):
            files.append(testpath)
            print(files)

            for f in files:
                data = np.loadtxt(f)
                print(data)

                wts = []
                mts = []
                maxl = 7

                bm = data[:, 1] == 0
                mutant = data[bm][:, 2]
                wt = data[np.logical_not(bm)][:, 2]
                wt = wt[:maxl]
                mutant = mutant[:maxl]

                wts.append(wt)
                mts.append(mutant)

                wts = np.array(wts)
                mts = np.array(mts)

                wavg = np.mean(wts, axis=0)
                mavg = np.mean(mts, axis=0)

                print(mavg)

    #files = pd.read_csv(files, delimiter='')  # Adjust delimiter as needed
    #triplicate_dfs.append(triplicate_df)


