import numpy as np
import matplotlib.pyplot as plt
import os

path="Z:\\Group Share\\Hiruni Jayasekera\\HSJ_MAE AqpZ CL\\Final\\WT_R224A Titrations All\\WT_R224A_15C\\"
os.chdir(path)

#Can import 2 files only
File1 = "DG_resultsN2_Fixed.txt"
#File2 = "Z:\\Group Share\\Dash\\Text Files\\PLIMB-NQ-40V-DMPC\\20230131 TD NQ PDMPC18C 10824.txt"

#convert to numpy array
iso1 = np.loadtxt(File1)
#iso2 = np.loadtxt(File2)

fig, (ax1) = plt.subplots()
ax1.plot(iso1[:,0], iso1[:,1]/np.amax(iso1[:,1])*100, 'red')
#ax1.plot(iso2[:,0], iso2[:,1]/np.amax(iso2[:,1])*((-1)*100), 'blue')

np.savetxt("DG_resultsN2_Fixed.csv", iso1, delimiter=",")
plt.xlabel("m/z")
plt.ylabel("%")
plt.legend()
plt.show()