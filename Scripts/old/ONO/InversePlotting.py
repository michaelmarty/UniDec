import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

# User inputs
chrom1 = 'TIC_decon.txt'
chrom2 = 'TIC_decon.txt'
label1 = 'Forward sequence'
label2 = 'Inverse sequence'


dir1 = 'Z:\\Group Share\\ONO\\CDMS\\5 bit HT data\\Oldre UniDec runs\\20240313 JDS Bgal 0o1g_l inj1o5 bit5 zp10 s1_unidecfiles'
dir2 = 'Z:\\Group Share\\ONO\\CDMS\\5 bit HT data\\Oldre UniDec runs\\20240313 JDS Bgal 0o1g_l inj1o5 bit-5 zp10 i1_unidecfiles'

dataLen = 460


# change the directory to the one with the data
os.chdir(dir1)

data = np.loadtxt(chrom1)
# plot the data
plt.subplot(211)
plt.plot(data[:, 0], data[:, 1], label = label1)

os.chdir(dir2)

data2 = np.loadtxt(chrom2)

# plot the data
plt.plot(data2[:, 0], data2[:, 1], label = label2)

plt.legend()

# average the data
data_avg = (data[:dataLen,1] + data2[:dataLen,1]) / 2

# average the data & evaluate with standard deviation
data_avg = (data[:dataLen,1] + data2[:dataLen,1]) / 2
st_dev = np.std(data_avg)
st_dev = round(st_dev, 2)

# subplot of averaged data
plt.subplot(212)
plt.plot(data[:dataLen, 0], data_avg, 'purple', label= 'Averaged, StDev: ' + str(st_dev))

#show the standard deviation at the top of the plot
#plt.title('Standard Deviation: ' + str(st_dev), fontsize=11)

plt.legend()
plt.show()