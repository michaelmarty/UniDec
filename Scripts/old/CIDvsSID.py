__author__ = 'Michael.Marty'

import numpy as np
import matplotlib.pyplot as plt
import os

cidspectra=["MTM_150729_ND_DMPC_01_CID_50_01_rawdata","MTM_150729_ND_DMPC_01_CID_75_01_rawdata",
            "MTM_150729_ND_DMPC_01_CID_100_01_rawdata","MTM_150729_ND_DMPC_01_CID_125_01_rawdata"]
cidlabels=["50V CID","75V CID","100V CID","125V CID"]

sidspectra=["MTM_150729_ND_DMPC_01_CID_50_01_rawdata","MTM_150729_ND_DMPC_01_CID_50_SID_25_01_rawdata",
            "MTM_150729_ND_DMPC_01_CID_50_SID_50_02_rawdata","MTM_150729_ND_DMPC_01_CID_50_SID_75_02_rawdata",
            "MTM_150729_ND_DMPC_01_CID_50_SID_100_01_rawdata"]
sidlabels=["0V SID","25V SID","50V SID","75V SID","100V SID"]

topdir="C:\\MassLynx\\Mike.PRO\\Data\\150729_SID"
os.chdir(topdir)
plt.figure()
plt.subplot(311)
maxdat=[]
for i,file in enumerate(cidspectra):
    fname=file+"_unidecfiles\\"+file+"_input.dat"
    data=np.loadtxt(fname)
    data[:,1]=data[:,1]/np.amax(data[:,1])-i*0.25
    plt.plot(data[:,0],data[:,1],label=cidlabels[i])

    fname2=file+"_unidecfiles\\"+file+"_peakparam.dat"
    pdata=np.loadtxt(fname2)
    marg=np.argmax(pdata[:,4])
    ms=pdata[marg]
    maxdat.append(ms)
#plt.gca().set_xlim(6000,14000)
plt.legend(fontsize=10)

plt.subplot(312)
smaxdat=[]
for i,file in enumerate(sidspectra):
    fname=file+"_unidecfiles\\"+file+"_input.dat"
    data=np.loadtxt(fname)
    data[:,1]=data[:,1]/np.amax(data[:,1])-i*0.25
    plt.plot(data[:,0],data[:,1],label=sidlabels[i])

    fname2=file+"_unidecfiles\\"+file+"_peakparam.dat"
    pdata=np.loadtxt(fname2)
    marg=np.argmax(pdata[:,4])
    ms=pdata[marg]
    smaxdat.append(ms)

plt.legend(fontsize=10)
plt.xlabel("m/z (Th)")
maxdat=np.array(maxdat)
smaxdat=np.array(smaxdat)
plt.subplot(313)
plt.plot((maxdat[:,0]-44088)/678.,maxdat[:,2],label="CID",marker="o")
plt.plot((smaxdat[:,0]-44088)/678.,smaxdat[:,2],label="SID",marker="s")
plt.xlabel("# of DMPC")
plt.ylabel("Charge")
plt.legend(loc=2)
plt.show()
