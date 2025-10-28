__author__ = 'Michael.Marty'

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import os
from scipy.stats import f

def CompareModels(e1,e2,p1,p2,n,thresh):
    plvl=0
    if(e1>e2 and p1>=p2):
        plvl=1
        model=1
    if(e1<e2 and p1<=p2):
        plvl=1
        model=0
    if(e1>e2 and p1<p2):
        v1=p2-p1
        v2=n-p2
        fval=(e1-e2)/(v1*(e2/v2))
        plvl=f.cdf(fval,v1,v2)
        if(plvl>=thresh):
            model=1
        else:
            model=0
    if(e2>e1 and p2<p1):
        v1=p1-p2
        v2=n-p1
        fval=(e2-e1)/(v1*(e1/v2))
        plvl=f.cdf(fval,v1,v2)
        if(plvl>=thresh):
            model=0
        else:
            model=1
    if(e1==e2 and p1>p2):
        plvl=1
        model=1
    if(e1==e2 and p1<p2):
        plvl=1
        model=0
    if(e1==e2 and p1==p2):
        plvl=1
        model=0
        print("Models identical")
    if(p1>n or p2> n):
        plvl=1
        print("Model Too Big")
    return plvl,model

def gaussian(x, height, center, width, offset):
    return height*np.exp(-(x - center)**2/(2*width**2)) + offset
def three_gaussians(x, h1, c1, w1, h2, c2, w2, h3, c3, w3, offset):
    if offset<0:
        return x*100000
    else:
        return (gaussian(x, h1, c1, w1, offset=0) +
            gaussian(x, h2, c2, w2, offset=0) +
            gaussian(x, h3, c3, w3, offset=0) + offset)

def four_gaussians(x, h1, c1, w1, h2, c2, w2, h3, c3, w3,h4,c4,w4, offset):
    if offset<0:
        return x*100000
    else:
        return (gaussian(x, h1, c1, w1, offset=0) +
            gaussian(x, h2, c2, w2, offset=0) +
            gaussian(x, h3, c3, w3, offset=0) +
                gaussian(x, h4, c4, w4, offset=0)+offset)

def two_gaussians(x, h1, c1, w1, h2, c2, w2, offset):
    return three_gaussians(x, h1, c1, w1, h2, c2, w2, 0,0,1, offset)

errfunc4 = lambda p, x, y: (four_gaussians(x, *p) - y)**2
errfunc3 = lambda p, x, y: (three_gaussians(x, *p) - y)**2
errfunc2 = lambda p, x, y: (two_gaussians(x, *p) - y)**2
errfunc1 = lambda p, x, y: (gaussian(x, *p) - y)**2

os.chdir("C:\\Data")

data = np.genfromtxt("C:\\Data\\lipidextracts.txt")
xaxis=np.arange(0,120)
proteins=["AmtB","AqpZ"]
pchoice=1

if pchoice==0:
    yaxis=np.average(data[:3],axis=0)
else:
    yaxis=np.average(data[3:],axis=0)
data=np.transpose([xaxis,yaxis])
#data=data[:100]

guess4 = [0.49, 9, 5, 0.6, 40, 10, 0.6, 80, 10,0.5,50,20, np.amin(yaxis)]
guess3 = [0.49, 9, 5, 0.6, 40, 10, 0.6, 80, 10, np.amin(yaxis)]
guess2 = [0.49, 9, 5, 0.5, 40, 20, np.amin(yaxis)]
guess1 = [1,data[np.argmax(data[:,1]),0],10,np.amin(yaxis)]
optim4, success = optimize.leastsq(errfunc4, guess4[:], args=(data[:,0], data[:,1]))
optim3, success = optimize.leastsq(errfunc3, guess3[:], args=(data[:,0], data[:,1]))
optim2, success = optimize.leastsq(errfunc2, guess2[:], args=(data[:,0], data[:,1]))
optim1, success = optimize.leastsq(errfunc1, guess1[:], args=(data[:,0], data[:,1]))
print(optim3[1], "+/-",optim3[2])
print(optim3[4], "+/-",optim3[5])
print(optim3[7], "+/-",optim3[8])
print(optim3[9])
print("")
print(optim2[1], "+/-",optim2[2])
print(optim2[4], "+/-",optim2[5])

np.savetxt(proteins[pchoice]+"_fits.txt",optim3)
np.savetxt(proteins[pchoice]+"_avg_lipid.txt",data)

e4=np.sum((four_gaussians(data[:,0], *optim4)-data[:,1])**2)
e3=np.sum((three_gaussians(data[:,0], *optim3)-data[:,1])**2)
e2=np.sum((two_gaussians(data[:,0], *optim2)-data[:,1])**2)
e1=np.sum((gaussian(data[:,0], *optim1)-data[:,1])**2)
p4=13
p3=10
p2=7
p1=4
#print CompareModels(e1,e2,p1,p2,len(data),0.99)
print(CompareModels(e2,e3,p2,p3,len(data),0.99))
print(CompareModels(e3,e4,p3,p4,len(data),0.99))


individuals=[]
for i in range(0,3):
    y=gaussian(data[:,0],optim3[i*3],optim3[i*3+1],optim3[i*3+2],0)
    individuals.append(y)
np.savetxt(proteins[pchoice]+"_gauss.txt",individuals)


plt.plot(data[:,0], data[:,1], lw=2, c='g', label='measurement')
plt.plot(data[:,0], three_gaussians(data[:,0], *optim3),
    lw=2, c='b', label='fit of 3 Gaussians')
plt.plot(data[:,0], two_gaussians(data[:,0], *optim2),
    lw=2, c='r', ls='--', label='fit of 2 Gaussians')

plt.plot(data[:,0], gaussian(data[:,0], *optim1),
    lw=2, c='y', ls=':', label='fit of 1 Gaussians')

plt.plot(data[:,0], four_gaussians(data[:,0], *optim4),
    lw=2, c='k', label='fit of 4 Gaussians')
plt.legend(loc='best')

magicnumbers=[78.8,50.8,7.5]
for v in magicnumbers:
    pass
    #plt.vlines(v,0,1)

plt.show()
#plt.savefig('result.png')