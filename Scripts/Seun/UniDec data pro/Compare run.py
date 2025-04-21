import os
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt
import matplotlib

file1 = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250407\\20250407TDMS-TSSPH-20uL2.raw"
file2 = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250407\\20250407TDMS-TSSPH-40uL2.raw"
#file3 = "Z:\\Group Share\\Dash\\Synapt Runs 20241218\\20240108_3.raw"

#read data
d1 = ImporterFactory.create_importer(file1)
d2 = ImporterFactory.create_importer(file2)
#d3 = ImporterFactory.create_importer(file3)

tic1 = d1.get_tic()
tic2 = d2.get_tic()
#tic3 = d3.get_tic()



plt.figure(figsize=(16, 5))
plt.plot(tic1[:,0], tic1[:,1], label="20uL")
plt.plot(tic2[:,0], tic2[:,1], label="40uL")
#plt.plot(tic3, label="Sample 3")
plt.xlabel("Time (s)")
plt.ylabel("Intensity")
plt.title("TIC Comparison")
plt.legend()
plt.show()
