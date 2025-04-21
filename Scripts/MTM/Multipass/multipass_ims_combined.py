import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
from copy import deepcopy

path = "C:\Data\Varun\STP_z14_allPasses_combined_raw.csv"
path2 = "C:\Data\Varun\STP_z15_allPasses_combined_raw.csv"

data = np.loadtxt(path, delimiter=",")
data2 = np.loadtxt(path2, delimiter=",")

cdata = ud.mergedata(data, data2)
cdata[:,1] += data[:,1]*0.2

ftdata = ud.fft(cdata)
b1 = ftdata[:,0] > 0.02
b2 = ftdata[:,0] < 0.06
b1 *= b2
ftdata = ftdata[b1]

maxindex = np.argmax(ftdata[:,1])
maxtime = 1/ftdata[maxindex,0]
print(maxtime)

ftdata2 = deepcopy(ftdata)
ftdata2[:,0] = 1/ftdata2[:,0]
ftdata2[:,1] -= np.amin(ftdata2[:,1])
ftdata2 = ud.normalize(ftdata2)

plt.figure()

plt.subplot(121)
#plt.title(path)
plt.plot(data[:,0], data[:,1]/np.amax(data[:,1]))
plt.plot(data2[:,0], data2[:,1]/np.amax(data2[:,1]))
plt.plot(cdata[:,0], cdata[:,1]/np.amax(cdata[:,1]))
plt.xlabel("Time (ms)")
plt.subplot(122)
plt.plot(data[:,0]-12, data[:,1]/np.amax(data[:,1]))
plt.plot(data2[:,0]-12, data2[:,1]/np.amax(data2[:,1])-1)
plt.plot(cdata[:,0]-12, cdata[:,1]/np.amax(cdata[:,1])-2)
plt.plot(ftdata2[:,0], ftdata2[:,1]-3)
plt.xlabel("Time per Cycle (ms)")
#plt.text(maxtime, ftdata[maxindex,1], str(maxtime))
plt.xlim(np.amin(ftdata2[:,0]), np.amax(ftdata2[:,0]))
#plt.vlines(maxtime, 0, ftdata[maxindex,1], linestyles='--', color="k")
plt.show()
