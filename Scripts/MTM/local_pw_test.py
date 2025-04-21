import sys

sys.path.append("C:\\Python\\UniDec")
import numpy as np
import matplotlib.pyplot as plt
import os


def mirrorer(datatop,index, window=100):
    d1=datatop[index:index+window,1]
    d2=datatop[index:index-window:-1,1]
    d3=d1*d2
    print(d3)
    return d3

toppath = "C:\\Data\\Lipid exchange new\\20170415_WCR_ND_P1B_P2A_1_10_ramp_0_200_20"
file = "20170415_WCR_ND_P1B_P2A_1_10_ramp_0_200_20_40.txt"
path=os.path.join(toppath,file)
data=np.loadtxt(path)
data=normalize(data)

index=100
for i in range(0,100):
    index=10000+i
    data2=mirrorer(data,index)
    data2/=np.amax(data2)
    data3=data[index-50:index+50,1]
    data3/=np.amax(data3)
    plt.plot(data2)
    plt.plot(data3)
    plt.show()

exit()


plt.plot(data[:,0],data[:,1])
plt.show()