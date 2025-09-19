from CDMS3 import *
slopes = []
medians_noise = []
slope_median_ratios = []
median_slope_ratios = []
signal_to_noise = []
names =[]
path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_BSA_CD-MS_10_min_maxinject_1ms.RAW"
mzrange = [3000, 7000]
irange = [100, 600]
cd = CDMS(path, mzrange=mzrange, irange=irange, slope=1, intensity_threshold=200, injection_time=0.001)
name ="BSA data"
slopes.append(cd.slope)
medians_noise.append(cd.median)
slope_median_ratios.append(cd.slope/cd.median)
median_slope_ratios.append(cd.median/cd.slope)
signal_to_noise.append(cd.signaltonoise)
names.append(name)
#print(cd.signaltonosie_max)
#print("Median:", cd.median, "Slope:", cd.slope, "Slope/Median:", cd.slope/cd.median)
#print("S/N:", cd.signaltonoise[1])

path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_ADH_ACN_50_IST_10min_Aqu_CD-MS.RAW"
mzrange = [5000, 10000]
irange = [200, 800]
cd = CDMS(path, mzrange=mzrange, irange=irange, slope=1, intensity_threshold=300, injection_time=0.001)
name = "ADH ACN"
slopes.append(cd.slope)
medians_noise.append(cd.median)
slope_median_ratios.append(cd.slope/cd.median)
median_slope_ratios.append(cd.median/cd.slope)
signal_to_noise.append(cd.signaltonoise)
names.append(name)


path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_ADH_CDMS_15min_aqu.RAW"
mzrange = [5000, 10000]
irange = [200, 800]
cd = CDMS(path, mzrange=mzrange, irange=irange, slope=1, intensity_threshold=200, injection_time=0.001)
name = "ADH"
slopes.append(cd.slope)
medians_noise.append(cd.median)
slope_median_ratios.append(cd.slope/cd.median)
median_slope_ratios.append(cd.median/cd.slope)
signal_to_noise.append(cd.signaltonoise)
names.append(name)

path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_AMTB_100V_CD-MS_8min_aqu_15ms.RAW"
mzrange = [4000, 12000]
irange = [100, 900]
cd = CDMS(path, mzrange=mzrange, irange=irange, slope=1, intensity_threshold=200, injection_time=0.015)
name = "AMTB data"
slopes.append(cd.slope)
medians_noise.append(cd.median)
slope_median_ratios.append(cd.slope/cd.median)
median_slope_ratios.append(cd.median/cd.slope)
signal_to_noise.append(cd.signaltonoise)
names.append(name)
#print("Median:", cd.median, "Slope:", cd.slope, "Slope/Median:", cd.slope/cd.median, )


path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_BSA__CD-MS_Aqu_ACN_10min.RAW"
mzrange = [3000, 7000]
irange = [100, 600]
cd = CDMS(path, mzrange=mzrange, irange=irange, slope=1, intensity_threshold=200, injection_time=0.001)
name = "BSA ACN"
slopes.append(cd.slope)
medians_noise.append(cd.median)
slope_median_ratios.append(cd.slope/cd.median)
median_slope_ratios.append(cd.median/cd.slope)
signal_to_noise.append(cd.signaltonoise)
names.append(name)

path = "Z:\\Group Share\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_GroEL_0V_CD-MS_50_100V_10min_.RAW"
mzrange = [10000, 12000]
irange = [200, 2000]
cd = CDMS(path, mzrange=mzrange, irange=irange, slope=1, intensity_threshold=200, injection_time=0.050)
name = "GroEL"
slopes.append(cd.slope)
medians_noise.append(cd.median)
slope_median_ratios.append(cd.slope/cd.median)
median_slope_ratios.append(cd.median/cd.slope)
signal_to_noise.append(cd.signaltonoise)
names.append(name)

print("Slopes:", slopes)
print("Medians:", medians_noise)
print("Slope/Median:", slope_median_ratios)
print("Median/Slope:", median_slope_ratios)
print("S/N:", signal_to_noise)
print("names:", names)
