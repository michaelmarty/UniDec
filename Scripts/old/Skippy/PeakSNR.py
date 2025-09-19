# Import chromatogram data and calculate noise level and SNR using a rolling standard deviation over the peak

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage
#import unidec.tools as ud

file = "Z:\\Group Share\\Skippy\\Projects\\HT\\Good data sets\\Bgal bit len\\1 bit\\20240411 Bgal bit1 3_2024-05-10-11-55-20_unidecfiles\\TIC_chrom.txt"
pkRange = (8.5, 11.5)

data = np.loadtxt(file)
bounds = [np.argmax(np.array(data[:,0]) > pkRange[0]), np.argmin(np.array(data[:,0]) < pkRange[1])]
#print(bounds)

pkData = data[bounds[0]:bounds[1], :]

# calculate rolling standard deviation over pkData
window = 10
rolling_std = np.zeros(len(pkData))
for i in range(len(pkData)):
    if i < window:
        rolling_std[i] = np.std(pkData[:i, 1])
    else:
        rolling_std[i] = np.std(pkData[i-window:i, 1])
np.nan_to_num(rolling_std, copy=False, nan=0)
rolling_std *= 3
#print("rolling SD", rolling_std)

# smooth the data with a slavitsky-golay filter
sm_data = scipy.signal.savgol_filter(pkData[:, 1], 5, 3)

noise = np.mean(rolling_std)
print('Mean Noise:', noise)
snr = np.max(sm_data)/noise
print('SNR:', snr)

# plot chromatogram
fig, ax1 = plt.subplots()
#set figure size
fig.set_size_inches(4, 4)
ax1.plot(pkData[:, 0], pkData[:, 1])
ax1.plot(pkData[:, 0], sm_data, color='k')
ax1.set_xlabel("Time (min)")
ax1.set_ylabel("Intensity")
ax1.plot(pkData[:, 0], rolling_std, color='r')
plt.show()



'''
ax2 = ax1.twinx()
ax2.plot(pkData[:, 0], rolling_std, color='r')
ax2.set_ylabel("Noise level")



plt.show()
'''

''' 
plt.figure()
plt.plot(data[:, 0], data[:, 1], label="Chromatogram")
#plt.plot(data[bounds[0]:bounds[1], 0], data[bounds[0]:bounds[1], 1], color = 'r')
plt.plot(pkData[:, 0], pkData[:, 1], color = 'r', label="Peak")
plt.xlabel("Time (min)")
plt.ylabel("Intensity")
plt.show()
'''