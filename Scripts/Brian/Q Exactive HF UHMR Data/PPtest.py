from unidec.modules.HTEng import *
import os
import matplotlib.pyplot as plt
import numpy as np

topdir = r"Z:\Group Share\BHT\Q Exactive HF Data\RPLC-MS\Agilent LC\20250828"
os.chdir(topdir)

file = "20250828_BHT_RPLC_C4_0.1mgmL_carbonicanhydrase_CDMS_6INJ_4_TIC.txt"
file = "20250828_BHT_RPLC_C4_0.1mgmL_carbonicanhydrase_CDMS_3INJ_1_TIC.txt"
file = "20250828_BHT_RPLC_C4_0.1mgmL_carbonicanhydrase_CDMS_6INJ_2_TIC.txt"

data = np.loadtxt(file, skiprows=1)
# Normalize data
data[:,1] = data[:,1] / np.max(data[:,1])
print(len(data))

eng = UniChromCDEng()
eng.config.demultiplexmode = "PP"
eng.ppskip = 0
eng.fulltime = data[:,0]
eng.fulltic = data[:,1]
eng.setup_demultiplex()
output, _ = eng.run_demultiplex(eng.fulltic, normalize=True)

plt.figure(figsize=(20,10))
plt.plot(data[:, 0], data[:, 1], label="Raw TIC")
eng.cycleindex = 0
plt.plot(data[eng.ppstartindex-eng.cycleindex:eng.ppstartindex+len(output)-eng.cycleindex,0], output, label="Demultiplexed")
plt.plot(data[:,0], eng.htkernel, label="HT Kernel")
plt.legend()
plt.show()