import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy as sp
from scipy import signal

from numpy.fft import fft, ifft
import sys, os, time
import HTfunctions as HTf

mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['errorbar.capsize'] = 3
mpl.rcParams['patch.force_edgecolor'] = True
mpl.rcParams['patch.facecolor'] = 'b'
mpl.rcParams['lines.markersize'] = 7
mpl.rcParams['font.size'] = 12

# change directory to 'Z:\mtmarty\Group Share\Skippy\Projects\HT LC\Data'
os.chdir('Z:/mtmarty/Group Share/Skippy/Projects/HT LC/Data')

data = np.loadtxt('20230906 BSA 3 bit.csv', delimiter = ',')        #20230828 JDS Mb 100ul_min 5bit 1sInj 60sCyc 1.csv'
x = data[:,0]
y = data[:,1]

#seq = [1, 1, 1, 0, 1, 0, 0]
seqTxt = '1110100'               #'0000100101100111110001101110101'
seq = []
for e in seqTxt:
    seq.append(int(e))
seq = np.array(seq, dtype = float)
mtx = HTf.createSeqMTX(seq)
seqLen = len(seq)

# Smooth and plot raw data


ySmth = signal.savgol_filter(y, 15, 3)
#print(len(ySmth))

#interpLen = HTf.next_power_of_2(len(ySmth))
#print(interpLen)

osLen = HTf.next_multiple_of_x(seqLen, len(ySmth))
osFact = int(osLen/seqLen)

interp_func = sp.interpolate.interp1d(x, ySmth, kind='cubic')
interpX = np.linspace(x[0], x[-1], osLen)
interpY = interp_func(interpX)
timeSpacing = x[-1]/len(x)
print(osLen)
print(len(interpY))

mtx = HTf.createSeqMTX(seq)

invMtx = np.linalg.inv(mtx)
invSeq = invMtx[:, 0]
invExp = HTf.expandSequence(invSeq, osFactor = osFact, singlePulse = True)  # expand sequence
xSeq = np.arange(len(invExp))*timeSpacing
transformedMtx = HTf.periodic_corr(HTf.normVector(interpY), invExp)

plt.subplot(2, 1, 1)
plt.plot(interpX, HTf.normVector(interpY))
plt.yticks([0,0.5,1])
plt.subplot(2, 1, 2)
plt.plot(interpX, transformedMtx)
plt.xlabel('Time (min)')
plt.yticks([0,0.5,1])

plt.show()