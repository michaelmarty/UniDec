from pyopenms import *
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from pyopenms import *
from imstools import *

dir = "C:\Data\Lipidomics\\220808_SONAR_IMS_data\HDMSE"
file = "220808_Ecoli_Lightsplash_IMS_MSMS_ResMode_Neg.mzML"
# file = "test_ms1.mzML"

os.chdir(dir)

exp = MSExperiment()
MzMLFile().load(file, exp)
filtered = MSExperiment()
print("Loaded")
rtt = 7.59 * 60
mzt = 721.5
mzwindow = 0.1
mzslice = [mzt - mzwindow, mzt + mzwindow]
rtwindow = 2 * 60
rtslice = [rtt - rtwindow, rtt + rtwindow]

for s in exp:
    ms_level = s.getMSLevel()
    if ms_level == 1:
        rt = s.getRT()
        if rtslice[0] < rt < rtslice[1]:
            print(rt)
            s2 = merge_to_dt(s, mzslice=mzslice)
            filtered.addSpectrum(s2)
print("loop done")
chrom = filtered.calculateTIC().get_peaks()
dt, ints = merge_rt(filtered)
plt.subplot(211)
plt.plot(dt, ints)
plt.subplot(212)
plt.plot(chrom[0], chrom[1])
plt.show()

'''
file = "C:\Data\Lipidomics\\220808_SONAR_IMS_data\SONAR\\220808_Ecoli_Lightsplash_SONAR_Res_Neg_sep.mzML"
os.chdir(os.path.dirname(file))
sonar_cal = pd.read_csv("SWATH_exp.txt", delimiter="\t")
print(sonar_cal)
low = sonar_cal["m/z begin"].to_numpy()
high = sonar_cal["m/z end"].to_numpy()
bin = sonar_cal["CycleID"].to_numpy()[1:]
b1 = low < mzt
b2 = high > mzt
b3 = b1*b2
binmatch = bin[b3]
print(binmatch)
exit()

file = "test.mzml"
exp = MSExperiment()
MzMLFile().load(file, exp)
filtered = MSExperiment()

for s in exp:
    rt = s.getRT()
    id = s.getNativeID()
    ids = id.split()
    d = dict()
    for thing in ids:
        thing = thing.split("=")
        d[thing[0]] = int(thing[1])
    func = d["function"]
    scan = d["scan"]
    s1 = scan - 1
    e = s1 % 200 + 1
    c = math.floor(scan / 200)
    if rt > rtt and e in binmatch:
        print(c, e)
        filtered.addSpectrum(s)'''
