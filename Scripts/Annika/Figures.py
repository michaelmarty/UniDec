import os

import numpy as np
import unidec.tools as ud

import unidec.modules.waters_importer.WatersImporter as WDI
import matplotlib.pyplot as plt

directory = "Z:\Group Share\Annika\APCI Data\Data for Figures"

os.chdir(directory)

file1 = "brainpolarextAPCI.txt"

data1 = np.loadtxt(file1)

i1 = ud.integrate(data1, 7.0, 8.0)[0]

cholmass = 369.3085
mrange = 0.05
range2 = 0.5
trange = 1

path = "Z:\\Group Share\\Annika\\APCI Data\Brain Extract\\20230620_APCI_BrainPolarExt_brainmethod_10uLinj.raw"
importer = WDI.WatersDataImporter(path)
importer = WDI.WatersDataImporter(path)

xic = importer.get_eic(cholmass, mrange)

fit = ud.auto_peak_width(xic, singlepeak=True)
print(fit)

maxtime = xic[np.argmax(xic[:,1]),0]
print(maxtime)

d1 = WDI.WatersDataImporter(path).get_data(time_range=(maxtime - trange, maxtime + trange))
d = ud.datachop(d1, cholmass-range2, cholmass + range2)

cholmass1 = 369.3085
mrange1 = 0.05
range21 = 0.5
trange1 = 0.05

path1 = "Z:\\Group Share\\Annika\\APCI Data\\221006_BrainLipid_Premier_FastDDA_Exclu_R1_run2_recali_Pos.raw"

importer1 = WDI.WatersDataImporter(path1)

importer1 = WDI.WatersDataImporter(path1)

xic1 = importer1.get_eic(cholmass1, mrange1)

fit1 = ud.auto_peak_width(xic1, singlepeak=True)
print(fit1)

maxtime1 = xic1[np.argmax(xic1[:,1]),0]
print(maxtime1)

d2 = WDI.WatersDataImporter(path).get_data(time_range=(maxtime1 - trange1, maxtime1 + trange1))
d3 = ud.datachop(d2, cholmass1-range21, cholmass1 + range21)



plt.subplot(212)
plt.plot(d[:,0], d[:,1], color="k", label="APCI")
plt.legend()
plt.xlabel("m/z")
plt.ylabel("Intensity")


plt.subplot(211)
plt.plot(d3[:,0], d3[:,1], label="ESI" , color="red")
plt.legend()
plt.xlabel("m/z")
plt.ylim(0,100000)
plt.ylabel("Intensity")

plt.show()