import time

from unidec.modules.HTEng import *

import unidec.tools as ud
# from modules.unidecstructure import DataContainer
import matplotlib.pyplot as plt
import scipy.fft as fft
import math
import os
import matplotlib

matplotlib.use('WxAgg')

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

dir = 'Z:\\Group Share\\ONO\\CDMS\\2024.03.13 5bit\\20240313 JDS Bgal 0o1g_l inj1o5 bit5 zp10 s1_unidecfiles'
os.chdir(dir)
fileName = 'TIC_chrom.txt'

configfile = ud.match_files(dir, "*conf.dat")[0]



data = np.loadtxt(fileName)
x = data[:,0]
y = data[:,1]



#plt.plot(x, y)
#plt.show()
#exit()

eng = HTEng()
eng.config.config_import(configfile)

print(len(x))

eng.fulltime = x
eng.fullscans = np.arange(len(x))

# eng.config.htbit = 5

eng.setup_ht()
dy, ddata = eng.htdecon(y)

rmsd = calc_rmsdz(x, dy, 9)

plt.plot(x, dy, 'k', linewidth=1.5, label= 'Raw RMSD: ' + str(round(rmsd,2)))
#plt.legend()

#plt.figure()
avgd, newdata = eng.masked_demultiplex_ht(y, 5, 1000)
#print(avgd)
rmsd = calc_rmsdz(x, avgd, 9)
#plt.plot(x, avgd/ np.amax(avgd) * np.amax(dy), 'r', linewidth=1.5, label='Masked Multiplexed\nRMSD: ' + str(round(rmsd,2)))
plt.plot(x, avgd/ np.amax(avgd) * np.amax(dy), 'r', linewidth=1.5, label='Masked Multiplexed \nRMSD: ' + str(round(rmsd,2)))

'''
# To use, must add return of "outdata" to HTEng.py. REMEMBER TO REMOVE OUTDATA WHEN DONE!!    
avgd, newdata, outdata = eng.masked_demultiplex_ht(y, 5, 100)
rmsd = calc_rmsdz(x, avgd, 9)
for i in range(4):
    plt.plot(x, outdata[i,:], linewidth=1.5, label='Artifact Chromatogram ' + str(i+1))
'''

plt.xlabel('Time (min)')
plt.ylabel('Intensity')
plt.xlim(0, 41)
plt.legend()
plt.show()

exit()

rmsds = []


# Plot different modes
plt.figure()
modes = ['rand', 'alt', 'con']
for i, e in enumerate(modes):
    avgd, newdata = eng.masked_demultiplex_ht(y, 5, 100, mode=e)
    rmsd = calc_rmsdz(x, avgd, 9)
    plt.subplot(3, 1, i+1)
    #plt.plot(x, dy, 'k', linewidth=1.5, label='Raw RMSD: ' + str(round(rmsd, 2)))
    plt.plot(x, avgd/ np.amax(avgd) * np.amax(dy), 'r', linewidth=1.5, label=str(e) + " RMSD: " + str(round(rmsd,2)))
    plt.legend()
    plt.ylabel('Intensity')

plt.xlabel('Time (min)')
#plt.show()
#exit()

'''
plt.figure()
plt.plot(x, dy, 'k', label= 'Raw RMSD: ' + str(round(rmsd,2)))
plt.xlabel('Time (min)')
plt.ylabel('Intensity')
'''


# Plot different window sizes
plt.figure()
plt.plot(x, dy, 'k', label= 'Raw')
winvals = [1, 10, 20]
# create color map from winvals
#colors = plt.cm.viridis(np.linspace(0, 1, len(winvals)))
for i, w in enumerate(winvals):

    avgd, newdata = eng.masked_demultiplex_ht(y, w, 100)
    avgd /= np.amax(avgd) / np.amax(dy)
    rmsd = calc_rmsdz(x, avgd, 9)
    rmsds.append(rmsd)
    rmsdtrue = calc_rmsdz(x, dy, startmin=7, endtime=8, data=avgd)
    plt.plot(x, avgd, linewidth=1.5, label='n: ' + str(w) + " RMSD: " + str(round(rmsdtrue,2)))
    #plt.plot(x, avgd, label='n: ' + str(w) + " RMSD: " + str(round(rmsd,2)) + " RMSDtrue: " + str(round(rmsdtrue,2)))

#plt.ylim(-100, 200)
plt.legend()
#plt.show()


# Plot with different number of generated artifact spectra
plt.figure()
plt.plot(x, dy, 'k', label= 'Raw')
num = [10, 50, 100, 1000]
for i, n in enumerate(num):

    avgd, newdata = eng.masked_demultiplex_ht(y, 5, n)
    avgd /= np.amax(avgd) / np.amax(dy)
    rmsd = calc_rmsdz(x, avgd, 9)
    rmsds.append(rmsd)
    rmsdtrue = calc_rmsdz(x, dy, startmin=7, endtime=8, data=avgd)
    plt.plot(x, avgd, linewidth=1.5, label='N: ' + str(n) + " RMSD: " + str(round(rmsdtrue,2)))

plt.legend()
plt.show()