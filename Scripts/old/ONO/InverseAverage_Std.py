import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

# User inputs
chrom1 = 'TIC_decon.txt'
chrom2 = 'TIC_decon.txt'
trial = '5-bit, trial 1'

#dir1 = 'Z:\\Group Share\\ONO\\CDMS\\2024.04.30 5bit\\20240430_ ONO_Bgal_bit5_zp10_fwdSeq_3_unidecfiles'
#ir2 = 'Z:\\Group Share\\ONO\\CDMS\\2024.04.30 5bit\\20240430_ ONO_Bgal_bit5_zp10_revSeq_3_unidecfiles'
dir1 = 'Z:\\Group Share\\ONO\\CDMS\\2024.03.13 5bit\\20240313 JDS Bgal 0o1g_l inj1o5 bit5 zp10 s1_unidecfiles'
dir2 = 'Z:\\Group Share\\ONO\\CDMS\\2024.03.13 5bit\\20240313 JDS Bgal 0o1g_l inj1o5 bit-5 zp10 i1_unidecfiles'

dataLen = 470 # roughly the scan number divided by compression factor



def calc_rmsdz(x, y, startmin, endtime=None, data= None):
    startindex = np.argmin(np.abs(x - startmin))
    if endtime is not None:
        endindex = np.argmin(np.abs(x - endtime))
    else:
        endindex = len(x)
    ydat = y[startindex:endindex]
    if data is None:
        d = np.zeros_like(ydat)
    else:
        d = data[startindex:endindex]
    return np.sqrt(np.mean((ydat-d) ** 2))


# change the directory to the one with the data, load the data
os.chdir(dir1)
data = np.loadtxt(chrom1)
#print(data)

os.chdir(dir2)
data2 = np.loadtxt(chrom2)

# RMSD calculation
rmsd1 = calc_rmsdz(data[:dataLen, 0], data[:dataLen, 1], 11, 40)
rmsd2 = calc_rmsdz(data2[:dataLen, 0], data2[:dataLen, 1], 11, 40)

'''
# standard deviations
st_dev1 = np.std(data)
st_dev1 = round(st_dev1, 2)
st_dev2 = np.std(data2)
st_dev2 = round(st_dev2, 2)
'''

# plot the first data
plt.subplot(211)
plt.plot(data[:, 0], data[:, 1], linewidth=1.5, label='Forward, RMSD: ' + str(round(rmsd1,2)))
# limit x axis to 7.5-15
plt.xlim(0, 41)

#plt.subplot(212)
# plot the second data
plt.plot(data2[:, 0], data2[:, 1], 'red', linewidth=1.5, label='Reverse, RMSD: ' + str(round(rmsd2,2)))

plt.legend()

# average the data & evaluate with standard deviation
data_avg = (data[:dataLen,1] + data2[:dataLen,1]) / 2
rmsd = calc_rmsdz(data[:dataLen, 0], data_avg, 11, 40)
#print(rmsd)
#st_dev = np.std(data_avg)
#st_dev = round(st_dev, 2)


# subplot of averaged data
plt.subplot(212)
plt.plot(data[:dataLen, 0], data_avg, 'purple', linewidth=1.5, label= 'Averaged, RMSD: ' + str(round(rmsd,2)))
plt.xlim(0, 41)
plt.xlabel('Time (min)')
plt.ylabel('Intensity')
plt.legend()

#plt.suptitle(trial)

plt.show()