import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib
import unidec.tools as ud

matplotlib.use('WxAgg')
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 0.75
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 16
rcParams['font.sans-serif'] = "Arial"

path = "Z:\Group Share\Group\Papers\Flash Extraction MP Nanodisc\Data\AqpZ\GC_POPC_POPE\GC_mixed_datavg_full.txt"
#path = "Z:\Group Share\Group\Papers\Flash Extraction MP Nanodisc\Data\AqpZ\GC_POPC_POPG\GC_mixed_datavg_full.txt"

data = np.loadtxt(path)
data[:,1] = data[:,1]/np.amax(data[:,1])
m1 = 98950
m1 = ud.data_extract(data, m1, 5, window=200)
print(m1)



pc = 760
pe = 718
avg = (pc+pe)/2.
actual = 749


plt.plot(data[:,0], data[:,1])
plt.vlines(m1, 0, 1)
for n in range(1, 8):
    com = ud.data_extract(data, m1+n*avg, 5, window=100)
    plt.vlines(m1+n*pe, 0, 1, color="purple")
    plt.vlines(m1+n*pc, 0, 1, color="k")
    plt.vlines(m1+n*avg, 0, 1, color="r")
    #plt.vlines(m1+n*actual, 0, 1)
    plt.vlines(com, 0, 1, color="g")
plt.xlim(98500, 110000)
#plt.xlim(99600,99800)
plt.show()