import numpy as np
import os
import matplotlib.pyplot as plt





#file2 = np.loadtxt('EcoliD1_Avg.txt')
#plt.plot(file2[:,0], file2[:,1], color="green")

#file3 = np.loadtxt('EcoliE3_Avg.txt')
#plt.plot(file3[:,0], file3[:,1], color="red")

#file4 = np.loadtxt('EcoliLipid_Demo.txt')
#plt.plot(file4[:, 0], file4[:, 1], color="blue")

#Time1 =[0, 0.5, 3,4, 5.4, 7, 7.01, 13,15, 15.1, 17, 18]
#MPA1 = [75,60,50,40,34,32.2,32.3,30,26.5,1,75,75]
#plt.plot(Time1, MPA1, color="red")

#Time2 = [0, 0.5, 17, 18]
#MPA2 = [75, 60, 75, 75]
#plt.plot(Time2, MPA2, color="blue")

path = "Z:\\mtmarty\\Group Share\\Guozhi Zhang\\Example Chrom"
os.chdir(path)
file1 = np.loadtxt('New Text Document.txt')
plt.plot(file1[:,0], file1[:,1], color="black")

#file2 = np.loadtxt('c1d400.txt')
#plt.plot(file2[:,0], file2[:,1], color="green")

#file3 = np.loadtxt('c2d250.txt')
#plt.plot(file3[:,0], file3[:,1], color="red")

#plt.title("TIC"#'Weighted Avg BPI', fontsize=25)
plt.xlabel('Time(min)', fontsize=12)
plt.ylabel("Base Peak Intensity", fontsize=12) #by total weight

#plt.legend(["1"], loc=1, fontsize=11)
plt.tight_layout()

plt.show()
