import sys

sys.path.append("C:\\Python\\UniDec")
import numpy as np
import matplotlib.pyplot as plt
import os


def noise_grabber(datatop):
    filter=[-0.25,-0.25,1,-0.25,-0.25]
    out=cconv(datatop[:,1],filter)
    return out

toppath = "C:\\Data\\Lipid exchange new\\20170415_WCR_ND_P1B_P2A_1_10_ramp_0_200_20"
file = "20170415_WCR_ND_P1B_P2A_1_10_ramp_0_200_20_40.txt"
path=os.path.join(toppath,file)
data=np.loadtxt(path)
data=normalize(data)
data2=deepcopy(data)

gdat=gaussian_backgroud_subtract(data,200)
out=noise_grabber(data)

plt.plot(data2[:,0],data2[:,1])
plt.plot(gdat[:,0],out)
plt.show()