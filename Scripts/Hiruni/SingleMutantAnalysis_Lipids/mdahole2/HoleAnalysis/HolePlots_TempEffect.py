import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



dir = "/mnt/z/Group Share/Group/VMD/AqpZMutants/WT_AA_PE_PG_CL_1/namd/hole"
dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_PE_PG_CL_1/namd\\hole"
os.chdir(dir)
avgs4 = pd.read_excel("profile.xlsx")
avgs4 = avgs4.to_numpy()
mr4 = pd.read_excel("minrad.xlsx")

dir = "/mnt/z/Group Share/Group/VMD/AqpZMutants/WT_AA_PEPGCL_310K_1/namd/hole"
dir = "Z:\\Group Share/Group/VMD/AqpZMutants/WT_AA_PEPGCL_310K_1/namd\\hole"
os.chdir(dir)
avgs5 = pd.read_excel("profile.xlsx")
avgs5 = avgs5.to_numpy()
mr5 = pd.read_excel("minrad.xlsx")


#plt.fill_between(avgs[:, 0], avgs[:, 1] - avgs[:, 2], avgs[:, 1] + avgs[:, 2], alpha=0.75, label="TOCL")
#plt.fill_between(avgs2[:, 0], avgs2[:, 1] - avgs2[:, 2], avgs2[:, 1] + avgs2[:, 2], alpha=0.75, label="POCL2")
#plt.fill_between(avgs3[:, 0], avgs3[:, 1] - avgs3[:, 2], avgs3[:, 1] + avgs3[:, 2], alpha=0.75, label="POPG")
plt.fill_between(avgs4[:, 0], avgs4[:, 1] - avgs4[:, 2], avgs4[:, 1] + avgs4[:, 2], alpha=0.75, label="PEPGCL_303K")
plt.fill_between(avgs5[:, 0], avgs5[:, 1] - avgs5[:, 2], avgs5[:, 1] + avgs5[:, 2], alpha=0.75, label="PEPGCL_310K", color="r")
#plt.fill_between(avgs6[:, 0], avgs6[:, 1] - avgs6[:, 2], avgs6[:, 1] + avgs6[:, 2], alpha=0.75, label="PEPG")
#plt.plot(avgs[:, 0], avgs[:, 1], label="TOCL2")
#plt.plot(avgs2[:, 0], avgs2[:, 1], label="POCL2")
#plt.plot(avgs3[:, 0], avgs3[:, 1], label="POPG")
#plt.plot(avgs4[:, 0], avgs4[:, 1], label="PEPGCL")
#plt.plot(avgs5[:, 0], avgs5[:, 1], label="PEPGCL_310K")
#plt.plot(avgs6[:, 0], avgs6[:, 1], label="PEPG")

plt.legend()
plt.ylim(0, 7)
plt.ylabel("Radius (A)")
#plt.xlim(-20, 20)
plt.xlabel("Reaction Coordinate")
plt.savefig(os.path.join(dir, "HolePlot1.PDF"))

plt.show()

plt.figure()
mrs = [mr4, mr5]
colors = ["b", "r"]
labels = [ "PEPGCL_303K", "PEPGCL_310K"]
for i, mr in enumerate(mrs):
    mr = mr[mr['frame'] < 20]
    avg = np.average(mr.to_numpy()[:, 1:], axis=1)
    std = np.std(mr.to_numpy()[:, 1:], axis=1)
    plt.fill_between(mr['frame'], avg - std, avg + std, alpha=0.75, label=labels[i], color=colors[i])
    print(labels[i], np.average(avg), np.average(std))

plt.legend()
plt.savefig(os.path.join(dir, "HolePlot2.PDF"))
plt.show()

