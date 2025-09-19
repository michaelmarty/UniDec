import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

dir = "Z:\\Group Share\\Skippy\\Projects\\HT\\Good data sets\\Bgal GroEL\\20240412 Bgal GroEL bit5 zp10 2_2024-04-12-02-39-31_unidecfiles"

chromfiles = ud.match_files_recursive(dir, "chrom.txt")
cdata = [np.loadtxt(file) for file in chromfiles]
#print(cdata)
print(len(cdata))

deconfiles = ud.match_files_recursive(dir, "decon.txt")
ddata = [np.loadtxt(file) for file in deconfiles]
#print(ddata)
print(len(ddata))

plt.figure(figsize=(6.5, 3.5), layout='constrained')
colors = ["b", "g", "k"]

for i in range(len(cdata)):
    print(i)
    plt.subplot(211)
    plt.plot(cdata[i][:, 0], cdata[i][:, 1], label="Chromatogram", color=colors[i])


    #print(ddata[1].shape)
    plt.subplot(212)
    plt.plot(ddata[i+1][:, 0], ddata[i+1][:, 1], label="Deconvoluted", color=colors[i])

plt.subplot(211)
plt.legend(['GroEL', 'Bgal', 'TIC'])
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)

plt.ylabel("(# of Ions)", labelpad=9, loc='bottom')

plt.subplot(212)

ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)
plt.xlabel("Time (min)")
plt.ylabel("Intensity", labelpad=4, loc='top')

#plt.show()

plt.figure(figsize=(2, 1.5), layout='constrained')

for i in range(len(cdata)):
    plt.subplot(212)
    plt.plot(ddata[i+1][:, 0], ddata[i+1][:, 1], label="Deconvoluted", color=colors[i])

ax = plt.gca()
ax.spines[['right', 'top', 'left']].set_visible(False)
plt.xlim(6, 14)
#remove y axis tics
ax.set_yticks([])
plt.show()