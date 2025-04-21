import os
import matplotlib.pyplot as plt
#from tools import *
import h5py
import sys
sys.path.append("/mnt/c/Python/unidec")

topdir = "C:\\VMD\\AmtB_Nanodiscs"
dirs=["POPG_E3_4NH2_all_atom","POPC_E3_4NH2_all_atom"]#,"PGPC_E3_4NH2_all_atom"]
#dirs=["PGPC_E3_4NH2_all_atom"]
colors=["r","b","g"]

plt.figure()

for i,d in enumerate(dirs):
    path=os.path.join(topdir,d,"analysis.hdf5")
    hdf = h5py.File(path)

    agroup = hdf.get("positions")
    avg_lip = get_dataset(agroup, "avg_lipid")
    avg_p = get_dataset(agroup, "avg_phos")

    agroup = hdf.get("all_atoms")
    md = get_dataset(agroup, "min_dists")

    agroup = hdf.get("rmsf")
    resnames = get_dataset(agroup, "resnames")
    print(md.shape)

    for j,t in enumerate(avg_lip):
        for k,l in enumerate(t):

            if resnames[k]=="POPG":
                color="b"
            else:
                color="r"
            alpha=(1/float(i+2))**2.
            if md[j,k]<3.5:
                if l[2]>0:
                    plt.subplot(1, 2, 1, aspect="equal")
                    plt.plot(l[0],l[1],marker="o",color=color, alpha=alpha)
                else:
                    plt.subplot(1, 2, 2, aspect="equal")
                    plt.plot(l[0], l[1], marker="o", color=color, alpha=alpha)
plt.show()