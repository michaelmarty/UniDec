from pyopenms import *
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from imstools import *
pathroot = "C:\Python\\feature-finder-ion-mobility"
sys.path.append(pathroot)
from feature_finder_im import *


dir = "C:\Data\Lipidomics\\220808_SONAR_IMS_data\HDMSE"
file = "220808_Ecoli_Lightsplash_IMS_MSMS_ResMode_Neg.mzML"
file = "test_ms1.mzML"

os.chdir(dir)

exp = MSExperiment()
MzMLFile().load(file, exp)


ff = FeatureFinderIonMobility()
nbins = 200
pptype = "custom"
peak_radius = 1
window_radius = 0.015
pp_mode = "int"  # ltr
ff_type = "centroided"  # centroided
dir = os.getcwd()
filter = 'none'
debug = False
bench = False

features = ff.run(exp, nbins, pptype, peak_radius, window_radius, pp_mode, ff_type, dir, filter, debug, bench)

FeatureXMLFile().store("Features.featureXML", features)
print('Found', features.size(), 'features')
exit()