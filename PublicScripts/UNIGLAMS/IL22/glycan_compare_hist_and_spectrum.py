import numpy as np
import os
import matplotlib.pyplot as plt

topdir = "D:\Data\Luis Genentech\IL22_Tryp"
os.chdir(topdir)

histfile = "IL22_SA15_tryp_glycan_list_4thsiteTheoreticalMassBruteForceHistdat1.txt"
histdat = np.loadtxt(histfile)

massfile = "SA15_mass.txt"
massdat = np.loadtxt(massfile)

histfile2 = "IL22_SA8_tryp_glycan_list_4thsiteTheoreticalMassBruteForceHistdat3.txt"
histdat2 = np.loadtxt(histfile2)

histfile3 = "IL22_SA4_tryp_glycan_list_4thsiteTheoreticalMassBruteForceHistdat3.txt"
histdat3 = np.loadtxt(histfile3)

plt.figure(figsize=(8,8))

plt.subplot(221)
plt.title("SA4")
plt.plot(histdat3[:,0]/1000., histdat3[:,1]/np.amax(histdat3[:,1]), label="Bottom-Up")
plt.xlim(90,120)
plt.ylabel("Relative Intensity")

plt.subplot(222)
plt.title("SA8")
plt.plot(histdat2[:,0]/1000, histdat2[:,1]/np.amax(histdat2[:,1]), label="Bottom-Up")
plt.xlim(90,120)

plt.subplot(223)
plt.title("SA15")
plt.plot(histdat[:,0]/1000, histdat[:,1]/np.amax(histdat[:,1]), label="Bottom-Up")
plt.xlim(90,120)
plt.xlabel("Mass (Da")
plt.ylabel("Relative Intensity")

plt.subplot(224)
plt.plot(histdat[:,0]/1000, histdat[:,1]/np.amax(histdat[:,1]), label="Bottom-Up")
plt.plot(massdat[:,0]/1000, massdat[:,1]/np.amax(massdat[:,1]), label="UNIGLAMS")
plt.xlim(90,120)
plt.xlabel("Mass (Da")

plt.legend()
plt.savefig("IL22_hist_vs_mass.pdf")
plt.show()