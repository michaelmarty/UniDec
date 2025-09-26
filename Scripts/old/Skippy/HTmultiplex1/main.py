import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy import signal

from numpy.fft import fft, ifft
import sys, os, time
import HTfunctions as HTf

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 12

# User inputs
csvFiles = ['Raw_622.csv', 'Raw_922.csv', 'Raw_1822.csv']
fileName = csvFiles[2]  # Select file to use
seqTxt = '0000100101100111110001101110101'      #1,1,1,0,1,0,0
osFactor = 200  # Oversampling factor
timeSpacing = 0.165576
smoothFactors = [15, 3]  # Savitzky-Golay filter parameters

# Convert sequence to array
seq = []
for e in seqTxt:
    seq.append(int(e))
seq = np.array(seq, dtype = float)
mtx = HTf.createSeqMTX(seq)

# Import data
data = np.loadtxt(fileName, delimiter = ',')
x = data[:,0]
y = data[:,1]

# Smooth and plot raw data
plt.subplot(2, 1, 1)
#print(fileName)
ySmth = signal.savgol_filter(y, smoothFactors[0], smoothFactors[1])
plt.plot(x, y, alpha = 0.25, label = "%s"%fileName)
plt.plot(x, ySmth)
#plt.legend()
#plt.show()

invMtx = np.linalg.inv(mtx)
invSeq = invMtx[:, 0]
invExp = HTf.expandSequence(invSeq, osFactor = osFactor, singlePulse = True)  # expand sequence
xSeq = np.arange(len(invExp))*timeSpacing
print(len(xSeq))

transformedMtx = HTf.periodic_corr(HTf.normVector(ySmth), invExp)   #HT function

#expMtx = HTf.createSeqMTX(invExp)
#HTout = fwht(HTf.normVector(ySmth))
#print(shape(HTout))

plt.subplot(2, 1, 2)
plt.plot(xSeq, transformedMtx)
plt.show()
