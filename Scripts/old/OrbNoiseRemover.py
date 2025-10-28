__author__ = 'Michael.Marty'


import numpy as np
import os
#import matplotlib.pyplot as plt

def LinearInterpolation(p0,p1,x):
    return p0[1]+(p1[1]-p0[1])*float(x-p0[0])/float(p1[0]-p0[0])

killrange=[8750.,8780.]
mid=(killrange[0]+killrange[1])/2.
span=(killrange[1]-killrange[0])/2.

path="C:\\MassLynx\\Mike.PRO\\Data\\150521\\mzML\\AqpZ_05_Ramp2"
os.chdir(path)
subdirs=[f for f in os.listdir(path) if os.path.isfile(f) and f[-4:]==".txt"]
for fname in subdirs:
    try:
        print(fname)
        data=np.loadtxt(fname)
        close=np.abs(data[:,0]-mid)<span
        slice=data[close]

        p1=slice[0]
        p2=slice[len(slice)-1]

        data[close,1]=[LinearInterpolation(p1,p2,x) for x in slice[:,0]]

        np.savetxt(fname,data)
    except:
        print("Failed",fname)