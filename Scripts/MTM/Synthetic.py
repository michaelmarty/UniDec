import numpy as np
import matplotlib.pyplot as plt
from unidec import tools as ud
import os
import unidec.modules.MassSpecBuilder as MSBuild

directory = "C:\\Users\\marty\\OneDrive - University of Arizona\\Grants\\NSF_23\\Figures"
os.chdir(directory)

isomasslist =[5000, 7800]
chargelist = [2, 3]

data, zrange = MSBuild.simple_spectrum([50000], resolution=100, zwidth=1, padding=0)

print(zrange)
exit()
data2, zrange2 = MSBuild.isotopic_spectrum(isomasslist, chargelist, intlist=[600, 0.8], crop=True, mono=True,
                                           mz_bin_size=0.1, padding=0)
data2[:,1] *= 0.8


data3 = np.concatenate([[[2400,0]],data2, data])

print(data2)
plt.subplot(121)
#plt.plot(data[:, 0], data[:, 1])
#plt.plot(data2[:, 0], data2[:, 1])
plt.plot(data3[:, 0], data3[:, 1])
plt.ylim(0, 1)
# Turn off right and top axes
ax = plt.gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
# Set y-axis ticks
ax.set_yticks([])
plt.xlabel("m/z")



data, zrange = MSBuild.simple_spectrum([50000],resolution=100, zwidth=1, zlist=[1])

data2, zrange2 = MSBuild.isotopic_spectrum(isomasslist, [1, 1], intlist=[0.6, 0.8], crop=True, mono=True,
                                           mz_bin_size=0.1)
data2[:,1] *= 0.8


data3 = np.concatenate([[[2400,0]],data2, data])

plt.subplot(122)
#plt.plot(data[:, 0], data[:, 1])
#plt.plot(data2[:, 0], data2[:, 1])
plt.plot(data3[:, 0], data3[:, 1])
plt.ylim(0, 1)
# Turn off right and top axes
ax = plt.gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
# Set y-axis ticks
ax.set_yticks([])
plt.xlabel("Mass")

plt.tight_layout()

#plt.savefig("MixFigMass.pdf")


plt.show()

exit()
data2, zrange2 = MSBuild.isotopic_spectrum([5000], [2], intlist=[0.6], crop=True, mono=True,
                                           mz_bin_size=0.01, padding=0)

plt.figure()

plt.plot(data2[:, 0], data2[:, 1])
plt.ylim(0, 1)
# Turn off right and top axes
ax = plt.gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
# Set y-axis ticks
ax.set_yticks([])
plt.xlabel("Mass")

plt.tight_layout()

plt.savefig("Zoom1.pdf")

plt.show()
