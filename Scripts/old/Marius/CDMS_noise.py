import numpy as np
import matplotlib.pyplot as plt
from unidec.modules import fitting
import unidec.tools as ud
from unidec import modules as CDEng
import scipy
from matplotlib import rcParams

#path = "C:\\Data\\CDMS\\20210309_MK_ADH_pos_CDMS_512ms_5min_50ms_pressure01.RAW"
path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\Replicates\\D1T0_CDMS_1"
#path = "C:\\Data\\CDMS\\02242021_MK_BSA__CD-MS_Aqu_ACN_10min.RAW"
eng = CDEng.UniDecCD()
eng.open_file(path)

ints = eng.darray[:,1]
#Filters intensity and m/z range
eng.filter_int([0,100])
eng.filter_mz([5000,10000])
noisedata = eng.farray



hist, edges = np.histogram(noisedata[:,1], bins=100)
edges = edges[1:]-(edges[1]-edges[0])/2.

smoothhist = scipy.ndimage.gaussian_filter1d(hist, 1)

grad = np.gradient(smoothhist)
grad2 = np.gradient(grad)

gradmin = edges[np.argmin(grad)]
print("Gradient Minimum", gradmin)

maxindex = np.argmax(hist)
grad2maxindex = np.argmax(grad2[maxindex:])
grad2max = edges[maxindex+grad2maxindex]
print("Grad2 Max", grad2max)


x = np.array([np.amin(noisedata[:,1]),np.amax(noisedata[:,0])])
y = grad2max + 0 * x
#weights = np.ones_like(noisedata[:,1]) / len(noisedata[:,1])
#for i in noisedata[:,1]:
    #normalized_data = (noisedata[:,1] - np.amin(noisedata[:,1])) / (np.amax(noisedata[:,1]) - np.amin(noisedata[:,1])) #normalize intensities

rcParams['font.size']=12
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']= 42
rcParams["figure.figsize"] = [12.00, 5.00]

plt.subplot(121)
plt.plot(edges, hist/np.amax(hist))
plt.plot(edges, grad/np.amax(hist)*2)
plt.plot(edges, grad2/np.amax(hist)*10)
plt.vlines(grad2max, 0, 1, color='r')
#Normalizes by dividing by the maximum intensity
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#plt.hist(noisedata[:,1], bins=100, density=True)
#plt.xticks(10)
#plt.yticks(10)
#plt.show()

#exit()
plt.subplot(122)
plt.hist2d(noisedata[:,0], noisedata[:,1], bins=200)
plt.plot(x, y, color='r')
plt.show()




exit()

file="C:\Data\CDMS\\noise.txt"
file="C:\Data\CDMS\\bsa_long_avg.txt"
path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\10202020\\201019_mk_BSA_DPPC_nanodiscs4_singleion_longer_aquisition_50ms.raw"

mzrange = [2000, 15000]
irange = [50, 600]
#path = "Z:\\Group Share\Marius Kostelic\\CD-MS\\10192020 DPPC\\201019_mk_cz_DPPC_24hrs_pos_singleion_25injecttime_pressure05.raw"
#mzrange = [7500, 15000]
#irange = [50, 600]
cd = CDMS(path, mzrange=mzrange, irange=irange)

data = cd.mergeddata
#data = np.loadtxt(file,skiprows=8)
boo1 = data[:,1]>0
hist, edges = np.histogram(data[boo1,1], bins=10000)
centers = edges[1:]-(edges[1]-edges[0])/2.

com, std = ud.center_of_mass(np.transpose([centers, hist]), relative_cutoff=0.8)
print(com, std)

fits, fitdat = fitting.exp_fit(centers, hist)
print(fits)

plt.plot(centers, hist)
plt.plot(centers, fitdat)
plt.show()