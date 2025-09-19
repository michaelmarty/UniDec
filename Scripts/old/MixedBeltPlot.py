import numpy as np
import os
import matplotlib.pyplot as plt
import fnmatch

topdir = "C:\\Data\\New\\"

belts = ["0", "0-1", "1",  "1-2","2", "0-2"]
colors = ["b","g","y","orange","r","purple"]

lims=[14420,14520]
max=0
volt=200
for i,b in enumerate(belts):
    folder = "POPC_D1T" + b + "m_HCD"
    print(folder)
    directory = os.path.join(topdir, folder)
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, '*_'+str(volt)+'.txt') or fnmatch.fnmatch(file, '*_'+str(volt)+'.0.txt'):
            data = np.loadtxt(os.path.join(directory, file))
            print(file)
            if max == 0:
                max=data[np.argmax(data[:,1]),0]
            lims=[max-25,max+85]

            b1=data[:,0]>lims[0]
            b2=data[:,0]<lims[1]
            b3=np.logical_and(b1,b2)
            #data=data[b3]

            data[:, 1] = data[:, 1] / np.amax(data[:, 1])
            data[:,1]-=0.25*i-1.25#i%2*(i+2)/2
            plt.plot(data[:, 0], data[:, 1], label=b, color=colors[i])

plt.legend(title="Extra T's")
plt.xlabel("m/z (Th)")
#plt.xlim(17460,17560)
plt.xlim(13500,21000)
#plt.savefig("C:\Users\Michael Marty\Box Sync\LabShared\MS_Data\mixed_belts_200v_zoom1.pdf")
plt.show()
