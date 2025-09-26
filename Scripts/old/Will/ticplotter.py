from UniDec import unidec_modules as dr
import numpy as np
import matplotlib.pyplot as plt
import os

dir="Z:\wresager\Bleaching\Figs\RAW fig"
os.chdir(dir)

files=["20170904_WCR_RHO_BL_POPC_2.RAW","20170901_WCR_RHO_BL_POPE_3.RAW","20170830_WCR_RHO_BL_cont_3.RAW"]
labels=["POPC","POPE","APO"]

for i,file in enumerate(files):
    data=np.array(dr.DataImporter(file).get_tic())
    plt.plot(data[:,0],data[:,1]/np.amax(data[:,1]),label=labels[i])
plt.legend()
plt.show()