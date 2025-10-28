from CDMS3 import *

path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_BSA__CD-MS_Aqu_ACN_10min.RAW"
mzrange = [3000, 7000]
irange = [1, 600]
cd = CDMS(path, mzrange=mzrange, irange=irange, slope=18.5, intensity_threshold=100)

print("Median:", cd.median, "Slope:", cd.slope, "Slope/Median:", cd.slope/cd.median, "Median Small:", cd.median_small)