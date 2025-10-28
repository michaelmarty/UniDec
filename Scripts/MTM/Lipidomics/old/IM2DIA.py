from pyopenms import *
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from imstools import *

dir = "C:\Data\Lipidomics\\220808_SONAR_IMS_data\HDMSE"
file = "220808_Ecoli_Lightsplash_IMS_MSMS_ResMode_Neg.mzML"
#file = "test.mzML"

os.chdir(dir)

exp = MSExperiment()
MzMLFile().load(file, exp)
#dtbins = get_im_bins(exp)

filtered = MSExperiment()

for s in exp:
    rt = s.getRT()
    ms_level = s.getMSLevel()
    print(rt)
    #fda = s.getFloatDataArrays()[0]
    #dtdata = np.array(fda.get_data())
    #mzdata, intdata = np.array(s.get_peaks())
    '''
    merge_spec(mzdata, dtdata, intdata)
    exit()
    plt.plot(mzdata, intdata)
    plt.show()
    # help(fdas[0])
    exit()
    
    '''
    if ms_level == 1:
        s2 = merge_to_mz(s)
        s2 = centroid(s2)
        filtered.addSpectrum(s2)

    if ms_level == 2:
        s2 = merge_to_mz(s)
        s2 = centroid(s2)
        filtered.addSpectrum(s2)
#exit()
MzMLFile().store("collapsed_centroid.mzML", filtered)
exit()

