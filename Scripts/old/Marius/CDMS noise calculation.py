import numpy as np
import os

#import sklearn.cluster

path = "C:\\UniDecPastedSpectra\\02242021_MK_AMTB_100V_CD-MS_8min_aqu_15ms_unidecfiles"
os.chdir(path)
# Load numpy compressed file
data = np.load(self.path, allow_pickle=True)['data']
# Assume m/z is in column 1 and intensity column 2
mz = data[:, 0]
intensity = data[:, 1]
# Look for scans in column 3
self.baseline = np.minimum(data[:, 1])
self.noise = np.median(data[:, 1])
SNR = intensity - self.baseline/ self.noise - self.baseline
#SNR = (NL - B) / (N - B)


print(self.noise)
print(self.baseline)
print(SNR)