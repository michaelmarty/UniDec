import os
from unidec import engine
import matplotlib.pyplot as plt
import numpy as np
# from sklearn.decomposition import PCA
from unidec.modules.isolated_packages import score_window as sw

directory = "C:\\Data\\AlgorithmTesting"
#directory = "C:\\Python\\UniDec3\\bin\\Example Data"
os.chdir(directory)

#file = "ADH.txt"
#file = "0.txt"
#file = "BSA.txt"
#file = "GroEL HCD unidec.txt"
file = "250313_AQPZ_POPC_100.txt"

'''
data = np.loadtxt(file)
noise = np.random.poisson(5,len(data))
ndat = np.transpose([data[:,0],noise])
file = "pnoise_5.txt"
np.savetxt(file,ndat)'''

eng = engine.UniDec()
eng.open_file(file, directory)
if True:
    eng.unidec_imports(everything=True)
else:
    eng.process_data()
    #eng.config.massub = 200000
    #eng.config.masslb = 50000
    eng.run_unidec()
#eng.config.peakthresh = 0.025
eng.pick_peaks()
# eng.pks_pca()
#eng.dscore(pow=1)
#eng.dscore(pow=2)
#eng.dscore(pow=3)
eng.dscore()

sw.score_plots(eng)

#eng.filter_peaks()

exit()
eng.dscore(xfwhm=3)

ztab = eng.data.ztab
for k, p in enumerate(eng.pks.peaks):
    m = p.zstack[0]
    ints = p.zstack[1:]
    summ = np.sum(ints, axis=0)
    sumz = np.sum(ints, axis=1)
    sumz /= np.amax(sumz)

    plt.figure()
    plt.subplot(131)
    plt.title("Combined Score: "+str(round(p.dscore*100,2))+"\nCharge Score: " + str(round(p.cs_score * 100, 2))
              )
    plt.plot(ztab, sumz)
    plt.plot(p.zdist[:, 0], p.zdist[:, 1], color="k")
    plt.subplot(132)
    for ival in ints:
        plt.plot(m, ival / np.amax(ints))
    plt.plot(p.mdist[:, 0], p.mdist[:, 1], color="k")
    plt.title("Mass:"+ str(p.mass)+"\nPeak Shape Score: " + str(round(p.mscore * 100, 2)))

    plt.subplot(133)
    for i, z in enumerate(ztab):
        v = p.mzstack[i]
        x = v[:, 0]
        y = v[:, 1]
        zvals = v[:, 2]
        plt.title("Uniqueness : " + str(round(p.uscore * 100, 2)) + "\nR Squared: " + str(round(p.rsquared * 100, 2))
                  )
        if True:
            try:
                x = x - np.amin(x)
                x = x / np.amax(x)
            except:
                pass
        plt.plot(x, y)#, color="k")
        plt.plot(x, zvals, color="k")

plt.show()
exit()

figManager = plt.get_current_fig_manager()
figManager.window.state('zoomed')

plt.show()
