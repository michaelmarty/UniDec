import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import matplotlib as mpl

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['font.size'] = 12

mpl.use("WxAgg")

dir = "Z:\\Group Share\\Skippy\\Projects\\HT\\Example data for MTM\\20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56_unidecfiles"

os.chdir(dir)

kfile = "20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56_htkernel.txt"
kdat = np.loadtxt(kfile)
kdat /= np.amax(kdat)

tfile = "TIC_chrom.txt"
tdat = np.loadtxt(tfile)
tdat[:,1] /= np.amax(tdat[:,1])

hfile = "TIC_HT_chrom.txt"
hdat = np.loadtxt(hfile)
hdat[:,1] /= np.amax(hdat[:,1])

padindex = 805
rollindex = -460

startindex = 350

for i in range(700, 1000):
    startindex = i

    kdat2 = np.concatenate([np.zeros(padindex), np.roll(kdat, startindex)])
    print(i)
    plt.figure(figsize=(8,6))
    plt.plot(kdat2, color="red")
    plt.plot(tdat[:,1], color="black")
    hindex = startindex + padindex + rollindex
    plt.plot(hdat[:hindex,1], color="blue")
    plt.xlabel("Scan Number")
    plt.ylabel("Normalized Intensity")
    plt.ylim(-1, 1)

    # Remove top and right axes
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    plt.savefig("C:\\Python\\Animations\\An" + str(startindex) + ".png", dpi=75, bbox_inches="tight")
    plt.close()
#plt.show()
