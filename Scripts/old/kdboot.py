__author__ = 'Michael.Marty'

import os
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

def model(lo,po,kd):
    return 0.5*(-1-((po-lo)/kd)+np.sqrt(4.*(lo/kd)+((lo-po)/kd-1.)**2.))

def errfun(kd,lo,po,rat,weightings=None):
    if weightings is None:
        return np.sum((rat-model(lo,po,kd))**2.)
    else:
        return np.sum(((rat-model(lo,po,kd))*weightings)**2.)

def fitmodel(kdguess,lo,po,rat,weightings=None):
    return opt.fmin(errfun,kdguess,args=(lo,po,rat,weightings),disp=0)[0]

def bootstrap(kdguess,lo,po,rat):
    numpts=1000
    hists=np.array([np.histogram(np.random.random_integers(0,len(lo)-1,size=len(lo)),bins=list(range(0,len(lo)+1)))[0] for i in range(0,numpts)])
    fits=[]
    for i in range(0,numpts):
        fits.append(fitmodel(kdguess,lo,po,rat,weightings=hists[i]))
    fits=np.array(fits)
    return fits


directory="C:\\cprog\\Joe\\Kd"
file="kD measurement.txt"

data=np.loadtxt(os.path.join(directory,file),skiprows=1)

po=data[:,0]
lo=data[:,1]
rat=data[:,7]

kdfit=fitmodel(0.5,lo,po,rat)
print("Best fit:",kdfit)

boots=bootstrap(kdfit,lo,po,rat)
bmean=np.mean(boots)
bstd=np.std(boots)
print("Bootstrapping mean and std deviation:",bmean,bstd)


plt.plot(lo,rat,marker="o")
plt.plot(lo,model(lo,po,kdfit),marker="o")
plt.show()