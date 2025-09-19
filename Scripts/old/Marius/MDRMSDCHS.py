import numpy as np
import matplotlib.pyplot as plt
import os
from unidec import modules as hd
import unidec.tools as ud
import matplotlib as mpl

def get_molig(path):
    return hd.get_param(path, "molig")


def get_mass_defect_grid(path):
    massaxis = hd.get_data(path, "ms_dataset", "mass_axis")
    massgrid = hd.get_data(path, "ms_dataset", "mass_grid")
    num = len(massgrid) / len(massaxis)
    massgrid = massgrid.reshape(num, len(massaxis))
    molig = get_molig(path)
    data_array = []
    for i in range(0, num):
        data1d, data2d, m1grid, m2grid, igrid = ud.kendrick_analysis(np.transpose([massaxis, massgrid[i]]), molig)
        dat2 = data1d[:, 1]
        dat2 = dat2 - np.amin(dat2)
        dat2 /= np.amax(dat2)
        data_array.append(dat2)
    kdat = np.array(data_array)
    return kdat, data1d[:, 0]


def get_mdg_list(paths):
    kall = []
    dall = []
    for path in paths:
        k, d = get_mass_defect_grid(path)
        kall.append(k)
        dall.append(d)
    kall = np.array(kall)
    dall = np.array(dall)
    kavg = np.average(kall, axis=0)
    kstd = np.std(kall, axis=0)
    davg = np.average(dall, axis=0)
    return kavg, kstd, davg, kall

folder="C:\\Data\\"
file="test.txt"
path=os.path.join(folder,file)
print(path)


f1 = ["Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 1 20180515\\20180515_MMK_D1TOCHS_15_1(2)_pos_actual_Ramp_0_200_20.hdf5",
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 20180518 chs replicate 2\\20170518_MK_D1T0_POPC_CHS 15_2_pos(actual2)_Ramp_0_200_20.hdf5",
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 20180522 CHS REPLICATE 3 ACTUAL\\20180522_mk_D1TO_PCCHS15_4_9C_REPLIC3_POS_RAMP_0_200_20.hdf5"]

f2 = ["Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 1 20180515\\20180515_MMK_D1TOCHS_15_1_IMID_pos_Ramp_0_200_20.hdf5",
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 20180518 chs replicate 2\\20170518_MK_D1T0_POPC_CHS 15_2_pos_IMID_Ramp_0_200_20.hdf5",
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 20180522 CHS REPLICATE 3 ACTUAL\\20180522_mk_D1TO_PCCHS15_4_9C_REPLIC3_POS_IMID_RAMP_0_200_20.hdf5"]

f3 = [
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 1 20180515\\20180515_MMK_D1TOCHS_15_0_1_Neg_RAMP_0_200_20.hdf5",
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 20180518 chs replicate 2\\20170518_MK_D1T0_POPC_CHS 15_2_NEG_Ramp_0_200_20.hdf5",
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 20180522 CHS REPLICATE 3 ACTUAL\\20180522_mk_D1TO_PCCHS15_4_9C_REPLIC3_NEG_ramp_0_200_20_180522150604.hdf5"]

f4 = ["Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 1 20180515\\20180515_MMK_D1TOCHS_15_1(3)_NEG_IMID_Ramp_0_200_20",
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 20180518 chs replicate 2\\20170518_MK_D1T0_POPC_CHS 15_2_NEG_IMID_Ramp_0_200_20",
      "Z:\Group Share\Marius Kostelic\chs dicholestryl chol erg amb replicates\CHS 20180522 CHS REPLICATE 3 ACTUAL\\20180522_mk_D1TO_PCCHS15_4_9C_REPLIC3_NEG_IMID_RAMP_0_200_20"]
f5 = [""]

k1, k1std, davg1, kall1 = get_mdg_list(f1)
k2, k2std, davg2, kall2 = get_mdg_list(f2)

k3, k3std, davg3, kall3 = get_mdg_list(f3)
#k4, k4std, davg4, kall4 = get_mdg_list(f4)

#k1, d1 = get_mass_defect_grid(f1)
#k2, d2 = get_mass_defect_grid(f2)
k5 = k2-k1
#k3 = kall2[1] - kall1[1]

k6 = k2-k3
#k4 = kall2[0] - kall1[0]
volts = np.arange(0, 220, 20)
print(volts)

plt.figure()
plt.subplot(441)
plt.title("CHS positive mode")
plt.contourf(davg1, volts, k1, 50)
plt.subplot(443)
plt.title("CHS positive imid")
plt.contourf(davg1, volts, k2, 50)
plt.subplot(442)
plt.title("CHS negative mode")
plt.contourf(davg1, volts, k3, 50)
plt.subplot(444)
plt.title("CHS negative imid")
plt.contourf(davg1, volts, k3, 50)

smax=np.amax([k1std,k2std,k3std,k3std])
print(smax)

norm=mpl.colors.Normalize(vmin=-smax, vmax=smax)
plt.subplot(445)
plt.title("Diff CHS pos +/- imid")
plt.contourf(davg1, volts, k5, 50, cmap="seismic", norm=norm)
plt.subplot(446)
plt.title("Diff CHS neg +/- imid")
plt.contourf(davg1, volts, k6, 50, cmap="seismic", norm=norm)

'''
norm=mpl.colors.Normalize(vmin=0, vmax=smax)
plt.subplot(427)
plt.title("STD DEV CHS pos")
plt.contourf(davg1, volts, k1std, 50, cmap="Reds", norm=norm)
plt.subplot(428)
plt.title("STD DEV CHS neg")
plt.contourf(davg1, volts, k3std, 50, cmap="Reds", norm=norm)
'''

#norm=mpl.colors.Normalize(vmin=0, vmax=smax)
plt.subplot(447)
plt.title("Rep Diff CHS pos")
plt.contourf(davg1, volts, kall1[0]-kall1[1], 50, cmap="seismic", norm=norm)
plt.subplot(448)
plt.title("Rep Diff CHS neg")
plt.contourf(davg1, volts, kall3[0]-kall3[1], 50, cmap="seismic", norm=norm)

plt.tight_layout()
plt.show()
