import matplotlib.pyplot as plt
import numpy as np
import os
from unidec.modules.hdf5_tools import *
from unidec.tools import datachop
import matplotlib as mpl

mpl.use('TkAgg')

dir = "C:\Data\Luis Genentech\AscQuad ILL22"
os.chdir(dir)

files = ["20231213_AscQuad_IL22Fc_Sialidase_PTR5_3uScans_100sCID.hdf5",
         "20231213_AscQuad_IL22Fc_Sialidase_PTR5_3uScans_100sCID_2nd.hdf5",
         "20231213_AscQuad_IL22Fc_Sialidase_PTR5_3uScans_100sCID_3rd.hdf5",]
         #"20231213_AscQuad_IL22Fc_Sialidase_PTR12_5uScans_100sCID.hdf5", ]

dataarrays = []

xvalues = [80000, 150000]
xvalues = [100000, 120000]

for file in files:
    data = get_mass_data(file)
    data = datachop(data, xvalues[0], xvalues[1])
    dataarrays.append(data)
    # Normalize
    data[:, 1] = data[:, 1] / np.max(data[:, 1])
    plt.plot(data[:, 0], data[:, 1], label=file[-22:])
plt.xlim(xvalues[0], xvalues[1])
plt.legend()

# Calculate correlation coefficient between data intensities
for i in range(0, len(dataarrays)):
    for j in range(i + 1, len(dataarrays)):
        min_len = min(len(dataarrays[i]), len(dataarrays[j]))
        corr = np.corrcoef(dataarrays[i][:min_len, 1], dataarrays[j][:min_len, 1])
        print(files[i], files[j], corr[0, 1])

plt.show()
