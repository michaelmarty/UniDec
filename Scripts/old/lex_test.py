import sys

sys.path.append("C:\\Python\\UniDec")
import numpy as np
import matplotlib.pyplot as plt
import os


def create_comb(lmass,bmass,ztab,ltab,data):
    mz=[]
    for l in ltab:
        mz.append((bmass+lmass*l+ztab)/ztab)
    mz=np.ravel(mz)

    dcomb=deepcopy(data)
    dcomb[:,1]=0

    pos=[nearest(data[:,0],m) for m in mz]
    dcomb[pos,1]=1

    return mz, dcomb

plt.figure()
plt.subplot(121)

for i in range(0,2):
    if i==1:
        toppath = "C:\\Data\\Lipid exchange new\\20170417_WCR_ND_P1B_P2A_1_10_ramp_0_200_20"
        file="20170417_WCR_ND_P1B_P2A_1_10_ramp_0_200_20_100.txt"
    else:
        toppath = "C:\\Data\\Lipid exchange new\\20170415_WCR_ND_P1B_P2A_1_10_ramp_0_200_20"
        file = "20170415_WCR_ND_P1B_P2A_1_10_ramp_0_200_20_100.txt"
    path=os.path.join(toppath,file)
    data=np.loadtxt(path)
    data=normalize(data)

    ztab=np.arange(8.,15.)
    freq1=ztab/760.076
    freq2=ztab/755.

    ltab=np.arange(50,150)

    lmass=760.076

    lrange=np.arange(750,762,0.1)
    brange=np.arange(0.,1.01,0.01)
    sumsall=[]
    for l in lrange:
        lmass = l
        sums = []
        for b in brange:
            bmass = b*760. + 44088
            comb, dcomb=create_comb(lmass,bmass,ztab,ltab,data)
            sums.append(np.sum(dcomb[:,1]*data[:,1]*data[:,1]))
        sumsall.append(sums)

    s=np.transpose(sumsall)
    print(np.shape(s),brange.shape,lrange.shape)
    plt.subplot(1,2,i+1)
    plt.contourf(lrange,brange*760.,s,N=200)
plt.show()

exit()

junk, fext2, data2, data3=fft_process(data,diffrange=[740,770])
data2=limit_data(data2,np.amin(freq1),np.amax(freq1))
plt.figure()
plt.subplot(121)
plt.plot(data2[:,0],data2[:,1])
plt.ylim(0,0.2)
for f in freq1:
    plt.vlines(f,0,1)
for f in freq2:
    plt.vlines(f,0,1,color="r")
plt.subplot(122)
plt.plot(fext2[:,0],fext2[:,1])

plt.figure()
plt.plot(data[:,0],data[:,1])
#plt.plot(dcomb[:,0],data[:,1]*dcomb[:,1]-1)
#for c in comb:
#    plt.vlines(c,0,1)
plt.show()

