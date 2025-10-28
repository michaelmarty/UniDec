import os

import numpy as np
import unidec.tools as ud

import unidec.modules.waters_importer.WatersImporter as WDI
import matplotlib.pyplot as plt

directory = "Z:\Group Share\Kimber\OzID - Current"

os.chdir(directory)

file1 = "spectestplz.txt"

data1 = np.loadtxt(file1)

i1 = ud.integrate(data1, 7.0, 8.0)[0]

cholmass = 369.3085
mrange = 0.05
range2 = 1
trange = 0.05

path = "Z:\\Group Share\\Annika\\APCI Data\\Brain Extract\\20230620_APCI_BrainPolarExt_brainmethod_10uLinj.raw"

importer = WDI.WatersDataImporter(path)

importer = WDI.WatersDataImporter(path)

xic = importer.get_eic(cholmass, mrange)

fit = ud.auto_peak_width(xic, singlepeak=True)
print(fit)

maxtime = xic[np.argmax(xic[:,1]),0]
print(maxtime)

d1 = WDI.WatersDataImporter(path).get_data(time_range=(maxtime - trange, maxtime + trange))
d = ud.datachop(d1, cholmass-range2, cholmass + range2)


#plt.plot(d[:, 0], d[:, 1])
#plt.show()

plt.subplot(211)
plt.plot(d[:,0], d[:,1], color="k")
plt.xlabel("m/z")
plt.ylabel("Intensity")


plt.subplot(212)
plt.plot(xic[:,0], xic[:,1], label="Cholesterol", color="red")
plt.legend()
plt.xlim(7,8)
plt.xlabel("Retention Time (min)")
plt.ylabel("Intensity")
plt.tight_layout()
plt.show()



plt.show()