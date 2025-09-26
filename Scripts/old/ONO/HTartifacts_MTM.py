import time

import numpy as np
import scipy.ndimage
import pandas as pd

from unidec.modules.HTEng import *

import unidec.tools as ud
# from modules.unidecstructure import DataContainer
import matplotlib.pyplot as plt
import scipy.fft as fft
import math
import os
import matplotlib

matplotlib.use('WxAgg')


def calc_rmsdz(x, y, startmin, endtime=None, data=None):
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
    return np.sqrt(np.mean((ydat - d) ** 2))


dir = 'Z:\\Group Share\\ONO\\CDMS\\5 bit HT data\\20240313 JDS Bgal 0o1g_l inj1o5 bit5 zp10 s1_unidecfiles'
dir = 'Z:\\Group Share\\Skippy\\Projects\\HT\\Good data sets\\Bgal bit len\\5 bit\\20240411 Bgal bit5 zp10 1_2024-04-26-03-44-37_unidecfiles'
os.chdir(dir)
fileName = 'TIC_chrom.txt'
fileName = 'TIC_decon.txt'

configfile = ud.match_files(dir, "*conf.dat")[0]

data = np.loadtxt(fileName)
x = data[:, 0]
y = data[:, 1]

# plt.plot(x, y)
# plt.show()
# exit()

eng = HTEng()
eng.config.config_import(configfile)
eng.config.HTksmooth = 0
eng.fulltime = x
eng.fullscans = np.arange(len(x))

eng.setup_ht()
dy, ddata = eng.htdecon(y)
indexes = eng.indexrange
ddata = ddata[indexes[0]: indexes[1]]
dy = dy[indexes[0]+eng.rollindex: indexes[1]+eng.rollindex]

kernel = eng.htkernel
sk = kernel / np.amax(kernel)
sk = sk == 1
sk = sk.astype(float)
#sk /= np.sum(sk)

dy2 = np.fft.irfft(np.fft.rfft(dy) * np.fft.rfft(sk)).real
if len(dy2) < len(ddata):
    dy2 = np.append(dy2, 0)

negs = np.array(dy2 < 0).astype(float)

dy3 = np.fft.irfft(np.fft.rfft(negs) * np.fft.rfft(sk)).real
if len(dy3) < len(ddata):
    dy3 = np.append(dy3, 0)

x2 = x[:len(kernel)]

#dy = dy/np.amax(dy) * np.amax(y)
#dy2 = dy2/np.amax(dy2) * np.amax(y)

plt.plot(x2, ddata, "b", label="Raw")
plt.plot(x2, dy, 'k')
#plt.plot(x2, kernel * np.amax(y) * sk, 'r', label='Kernel')
#plt.plot(x2, rmsd_shift/np.amax(rmsd_shift) * np.amax(y), 'g', label='RMSD Shift')
plt.plot(x2, dy2, 'r')
#plt.plot(x2, sk/np.amax(sk)*np.amax(y), "g")
#plt.plot(x2, negs/np.amax(negs) * np.amax(y), 'g')

plt.xlabel('Time (min)')
plt.ylabel('Intensity')
plt.show()

exit()
