import os
import numpy as np
import unidec.tools as ud

from UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt

path = "Z:\\Group Share\\Annika\\APCI Data\\Equisplash\\20230801_5_1_Eqd7_NotagPOPCCHOL_1_32uM_control.raw"

importer = ImporterFactory.create_importer(path)

cmass = 369.3185
d7cmass = 376.3511

popc_mass = 760.5529
d7pc_mass = 753.5815

importer2 = ImporterFactory.create_importer(path)
data = importer.get_eic(cmass, tolerance=0.05)
data2 = importer2.get_eic(d7cmass, tolerance=0.05)

importer3 = ImporterFactory.create_importer(path)
data3 = importer3.get_eic(d7pc_mass, tolerance=0.05)
importer4 = ImporterFactory.create_importer(path)
data4 = importer4.get_eic(popc_mass, tolerance=0.05)

startime = 6.4
endtime=7

data = ud.datachop(data, startime, endtime)
data2 = ud.datachop(data2, startime, endtime)

plt.subplot(211)
plt.plot(data[:,0], data[:,1], color="k", label = "Cholesterol")
plt.plot(data2[:,0], data2[:,1], color = "blue", label = "D7 Cholesterol")
plt.xlabel("Retention Time (mins)", fontdict={"size":14})
plt.ylabel("Intensity")
#plt.xlim(6.4,7)
plt.ylim(0,11000)
plt.legend()

plt.subplot(212)
plt.plot(data3[:,0], data3[:,1], color = "green", label = "D7 PC")
plt.plot(data4[:,0], data4[:,1], color = "red" , label = "POPC")
plt.xlabel("Retention Time (mins)")
plt.ylabel("Intensity")
plt.xlim(6.5,8.25)
plt.ylim(0,440000)
plt.legend()

plt.tight_layout()
plt.savefig("EIC.pdf")
plt.show()