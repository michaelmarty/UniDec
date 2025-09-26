__author__ = 'Michael.Marty'

import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os
import subprocess
import string
from copy import deepcopy

def logistic(x,x0,k):
    return 1./(1.+np.exp(-k*(x-x0)))

def SkewedGauss(x,x0):
    m=x0[0]
    s=x0[1]
    a=x0[2]
    try:
        skew=x0[3]
    except:
        skew=0.
    if s<0 or skew>1.5:
        return 1e10
    else:
        return a*np.exp(-(x-m) * (x-m) / (2 * s * s))*logistic(x,m,skew)

def f0(x,x0):
    m=x0[0]
    s=x0[1]
    a=x0[2]
    if s<0 or a<0:
        return 1e10
    else:
        return a/(s*np.sqrt(2*np.pi))*np.exp(-(x-m) * (x-m) / (2 * s * s))#

def ExpDecay(x,x0):
    a=x0[1]
    b=x0[2]
    if a<0 or b<0:
        return 1e10
    else:
        return a*np.exp(-b*x)

def error2(x0,xdata,ydata,x1):
    test=np.array([(f0(x,[x1[0],x1[1],x0[0]])+ExpDecay(x,x0)) for x in xdata])
    #test=test/np.sum(test)
    return np.sum(np.sqrt((test-ydata)**2))

def gaussfit2(xdata,dist,mz,skewflag):
    mid=model(mz,p0mid)
    sig=model(mz,p0sig)
    x1=[mid,sig]
    a=np.exp(-mz/2000.)
    guess=np.array([0.1,a,1.])

    #fits=minimize(error2,guess,args=(xdata,dist,x1),method="Powell").x
    #fits=minimize(error2,guess,args=(xdata,dist,x1),method="SLSQP").x
    fits=minimize(error2,guess,args=(xdata,dist,x1),method="Nelder-Mead").x
    #np.mean(np.array([fits,fits1,fits2]),axis=0)
    #print fits
    return fits

def error(x0,xdata,ydata):
    test=np.array([SkewedGauss(x,x0) for x in xdata])
    return np.sum(np.sqrt((test-ydata)**2))

def gaussfit(xdata,dist,mz,skewflag):
    mid=xdata[np.argmax(dist)]
    a=np.amax(dist)
    a=0.1
    skew=0.0
    sig=2500.0
    a2=1.
    b=1.
    mid=model(mz,p0mid)
    sig=model(mz,p0sig)

    if skewflag==0:
        guess=np.array([mid,sig,a])
    else:
        guess=np.array([mid,sig,a,skew])
    #print guess

    #fits=minimize(error,guess,args=(xdata,dist),method="Powell").x
    #fits=minimize(error,guess,args=(xdata,dist),method="SLSQP").x
    fits=minimize(error,guess,args=(xdata,dist),method="Nelder-Mead").x
    #np.mean(np.array([fits,fits1,fits2]),axis=0)
    #print fits
    return fits

def model(mz,x):
    a=x[0]
    b=x[1]
    c=x[2]
    #d=x[3]
    #e=x[4]
    #f=x[5]
    #g=x[6]
    return a+b*(mz**c)#+d*np.exp(-mz*e)
    #return a+b*mz+c*mz**2+d*mz**3+e*mz**4+f*mz**5+g*mz**6

def model2(mz,x):
    a=x[0]
    b=x[1]
    return a*np.exp(-mz*b)


modelsize=3
method="Nelder-Mead"
#method="Powell"
#method="SLSQP"
p0all=[]
p0mid=[-3.32730372e-01,   6.14300712e-04,   1.00101841e+00]
p0sig=[ 0.08054858,  0.02428673,  0.50527904]
p0skew=[ -8.96170668e-02,   1.88356035e+00,   8.84233893e-04]

def errfun(x,xdata,ydata):
    yout=np.array([model(mz,x) for mz in xdata])
    return np.sum((yout-ydata)**2)

def errfun2(x,xdata,ydata):
    yout=np.array([model2(mz,x) for mz in xdata])
    return np.sum((yout-ydata)**2)


def sigfit(massrange,fit):
    p0=p0sig
    fits=minimize(errfun2,p0,args=(massrange,fit[:,1]),method=method).x
    print(fits)
    return fits

def midfit(massrange,fit):
    p0=p0mid
    fits=minimize(errfun,p0,args=(massrange,fit[:,0]),method=method).x
    print(fits)
    return fits

def skewfit(massrange,fit):
    p0=p0skew
    p0[0]=1.8
    p0[1]=1/2000.
    p0[2]=0.
    fits=minimize(errfun2,p0,args=(massrange,fit[:,2]),method=method).x
    #fits=np.append(fits,1.)
    print(fits)
    return fits

def safedivide(a,b):
    c=deepcopy(b)
    c[b!=0]=a[b!=0]/b[b!=0]
    return c

p0mid=[-3.32730372e-01,   6.14300712e-04,   1.00101841e+00]
p0sig=[ 0.08054858,  0.02428673,  0.50527904]
p0alpha=[  1.07875918e+00,   1.01017049e-03]
p0beta=[ 1.75546608e+00,   5.51128395e-04]
p0all=np.concatenate((p0alpha,p0beta,p0mid,p0sig))
#p0all=[  1.01796305e+00,   1.31066023e-03,   2.29009362e+00,   8.63732533e-04,-4.99937235e-01,   7.02808201e-04,   9.89579666e-01,   4.36275018e-01,9.75560419e-03,   5.86144765e-01]
p0all=[  1.00840852e+00,   1.25318718e-03,   2.37226341e+00,   8.19178000e-04, -4.37741951e-01,   6.64992972e-04,   9.94230511e-01,   4.64975237e-01,1.00529041e-02,   5.81240792e-01]
p0mid=p0all[4:7]
p0sig=p0all[7:10]
p0beta=p0all[2:4]
p0alpha=p0all[0:2]
def finalmodel(isorange,mz):
    mid=model(mz,p0mid)
    sig=model(mz,p0sig)
    alpha=model2(mz,p0alpha)
    beta=model2(mz,p0beta)
    dist=f0(isorange,[mid,sig,(1-alpha)])+ExpDecay(isorange,[0,alpha,beta])
    return dist/np.sum(dist)

def errfunfinal(x,massrange,isorange,intensities):
    error=0
    for i in range(0,len(massrange)):
        mz=massrange[i]
        mid=model(mz,x[4:7])
        sig=model(mz,x[7:10])
        alpha=model2(mz,x[0:2])
        beta=model2(mz,x[2:4])
        dist=f0(isorange,[mid,sig,(1-alpha)])+ExpDecay(isorange,[0,alpha,beta])
        error+=np.sum((dist-intensities[:,i])**2)
    return error

def MinimizeFinal(p0,massrange,isorange,intensities):
    fits=minimize(errfunfinal,p0,args=(massrange,isorange,intensities),method=method).x
    print(fits)
    return fits




class Isotopes:
    def __init__(self,*args,**kwargs):
        os.chdir("C:\\Python\\UniDec\\old")
        self.isoexe="C:\\MinGW\\isotopes\\isotopes.exe"
        self.isotab="C:\\MinGW\\isotopes\\isotopestable.txt"
        self.name="isotopegrid.txt"
        self.name2="isotopefits.txt"
        self.name3="isoparam.txt"
        self.massname="masses.txt"
        self.massavgine=111.1254
        self.avgine=np.array([4.9384,7.7583,1.3577,1.4773,0.0417])
        self.isotopes=np.array([[[12.,98.90],[13.0033548,1.10],[0,0],[0,0]],
                           [[1.0078250,99.985],[2.0141018,0.015],[0,0],[0,0]],
                           [[14.0030740,99.634],[15.0001089,0.367],[0,0],[0,0]],
                           [[15.9949146,99.772],[16.9991315,0.038],[17.9991605,0.2],[0,0]],
                           [[31.9720707,95.02],[32.9714585,0.75],[33.9678668,4.21],[35.9760809,0.02]]])
        self.nindex=[2,2,2,3,4]
        self.testmass=16420
        output=self.makeisotopes(self.testmass)
        print(output)
        exit()



        self.massrange=np.logspace(2.5,5.2,num=50)
        print("Mass Range: ",min(self.massrange),max(self.massrange))


        self.mlength=len(self.massrange)
        output=self.makeisotopes(max(self.massrange))
        self.maxiso=np.amax(output[:,0])+1
        self.isorange=np.arange(0,self.maxiso,1)
        self.intensity=None
        self.fits=None
        self.skewflag=2
        print("Number of Masses: ",self.mlength)

    def SaveIsotopes(self):
        temp=np.loadtxt(self.massname)
        if np.all(temp==self.massrange):
            self.LoadIsotopes()

        else:
            self.intensity=np.zeros((self.maxiso,self.mlength))
            for i in range(0,self.mlength):
                output=self.makeisotopes(self.massrange[i])
                for o in output:
                    self.intensity[o[0],i]=o[1]
            np.savetxt(self.name,self.intensity)
            np.savetxt(self.massname,self.massrange)
            print("Saved Intensity Grid: ",self.name)

    def LoadIsotopes(self):
        temp=np.loadtxt(self.name)
        if temp.shape==(self.maxiso,self.mlength):
            self.intensity=temp
            print("Intensity Grid Loaded")

    def LoadFits(self):
        temp=np.loadtxt(self.name2)
        if temp.shape==(self.mlength,4):
            self.fits=temp
            print("Fits Loaded")


    def makemass(self,testmass):
        num=testmass/self.massavgine*self.avgine
        intnum=[int(round(n)) for n in num]
        minmassint=(intnum*self.isotopes[:,0,0]).sum()
        formula=""
        if intnum[0] is not 0:
            formula=formula+"C"+str(intnum[0])
        if intnum[1] is not 0:
            formula=formula+"H"+str(intnum[1])
        if intnum[2] is not 0:
            formula=formula+"N"+str(intnum[2])
        if intnum[3] is not 0:
            formula=formula+"O"+str(intnum[3])
        if intnum[4] is not 0:
            formula=formula+"S"+str(intnum[4])
        return formula,minmassint

    def makeisotopes(self,testmass):
        formula,minmass=self.makemass(testmass)
        try:
            output=subprocess.check_output([self.isoexe,formula,self.isotab,"-f 0"]).splitlines()
        except:
            print(formula)
            print(testmass)
            exit()
        outarray=np.array([line.split("\t") for line in output]).astype(np.float)
        distribution=outarray[:,1]
        massdiff=outarray[:,0]-minmass
        deltas=[massdiff[i]-massdiff[i-1] for i in range(1,len(massdiff))]
        meandiff=np.mean(deltas)
        print(meandiff)
        try:
            diffpos=[int(round(m/meandiff)) for m in massdiff]
        except:
            print(formula)
            print(testmass)
            print(outarray)
            exit()
        array=np.array([diffpos,distribution]).transpose()
        return array


    def Plot(self):
        self.SaveIsotopes()
        self.FitIntensities()
        plt.figure()
        plt.contourf(self.massrange,self.isorange,np.sqrt(self.intensity),50)
        plt.plot(self.massrange,self.fits[:,0],color="k")
        plt.plot(self.massrange,self.fits[:,0]+self.fits[:,1],color="k")
        plt.plot(self.massrange,self.fits[:,0]-self.fits[:,1],color="k")
        plt.show()
        print("plotted")

    def FitFinal(self,*args):
        if self.intensity is None:
            self.SaveIsotopes()
        print(errfunfinal(p0all,self.massrange,self.isorange,self.intensity))
        self.fitsfinal=MinimizeFinal(p0all,self.massrange,self.isorange,self.intensity)
        print(errfunfinal(self.fitsfinal,self.massrange,self.isorange,self.intensity))


    def FitIntensities(self,*args):
        if self.intensity is None:
            self.SaveIsotopes()
        self.fits=np.zeros((self.mlength,1+self.skewflag))
        for i in range(0,self.mlength):
            dist=self.intensity[:,i]
            fit=gaussfit2(self.isorange,dist,self.massrange[i],self.skewflag)
            if self.skewflag==2:
                mid=model(self.massrange[i],p0mid)
                sig=model(self.massrange[i],p0sig)
                fitdat=[f0(x,[mid,sig,fit[0]])+ExpDecay(x,fit) for x in self.isorange]
            else:
                fitdat=[SkewedGauss(x,fit) for x in self.isorange]
            error=np.sum((dist-fitdat)**2)/np.sum(dist)
            print(error)
            #if error<0.1:
            self.fits[i]=fit
            if "plot" in args:
                plt.figure()
                plt.plot(self.isorange,dist)
                plt.plot(self.isorange,fitdat,marker=".")
                plt.show()
        np.savetxt(self.name2,self.fits)
        pass

    def Train(self):
        if self.fits is None:
            self.FitIntensities()
        plt.figure()
        midfits=midfit(self.massrange,self.fits)
        midfitdat=[model(mz,midfits) for mz in self.massrange]

        #plt.figure()
        plt.subplot(3,1,1)
        plt.plot(self.massrange,self.fits[:,0],marker=".")
        plt.plot(self.massrange,midfitdat,color="r")
        #plt.show()


        sigfits=sigfit(self.massrange,self.fits)
        sigfitdat=[model2(mz,sigfits) for mz in self.massrange]

        #plt.figure()
        plt.subplot(3,1,2)
        plt.plot(self.massrange,self.fits[:,1],marker=".")
        plt.plot(self.massrange,sigfitdat,color="r")


        if self.skewflag==1:
            skewfits=skewfit(self.massrange,self.fits)
            skewfitdat=[model(mz,skewfits) for mz in self.massrange]

            plt.subplot(3,1,3)
            plt.plot(self.massrange,self.fits[:,3],marker=".")
            plt.plot(self.massrange,skewfitdat,color="r")
            #plt.show()
        if self.skewflag==2:
            skewfits=skewfit(self.massrange,self.fits)
            skewfitdat=[model2(mz,skewfits) for mz in self.massrange]

            plt.subplot(3,1,3)
            plt.plot(self.massrange,self.fits[:,2],marker=".")
            plt.plot(self.massrange,skewfitdat,color="r")
            #plt.show()
        else:
            skewfits=midfits*0

        plt.show()
        self.trained=np.ravel(np.array([midfits,sigfits,skewfits]))
        print("Trained",self.trained)
        np.savetxt(self.name3,self.trained)

    def Test(self):
        self.numtest=10
        self.testmassrange=np.sort(10**(np.random.uniform(np.log10(min(self.massrange)),np.log10(max(self.massrange)),self.numtest)))
        print(self.testmassrange)
        plt.figure()
        self.trained=np.loadtxt(self.name3)
        for i in range(0,self.numtest):
            var=len(self.trained)/3
            mid=model(self.testmassrange[i],self.trained[0:var])
            sig=model(self.testmassrange[i],self.trained[var:var*2])
            if self.skewflag==1:
                skew=model(self.testmassrange[i],self.trained[var*2:var*3])
            else:
                skew=0
            testdist=np.array([SkewedGauss(x,[mid,sig,1,skew]) for x in self.isorange])

            testintensity=np.zeros(self.maxiso)
            output=self.makeisotopes(self.testmassrange[i])
            for o in output:
                testintensity[o[0]]=o[1]
            actualdist=testintensity
            testdist=testdist*np.amax(actualdist)/np.amax(testdist)
            #plt.figure()
            plt.plot(self.isorange,testdist,color="r")
            plt.plot(self.isorange,actualdist,color="b")
            #plt.show()

            error=0.
            maxval=max(actualdist)
            num=0.
            for i in range(0,len(self.isorange)):
                if actualdist[i]>maxval*0.001 and actualdist[i]!=0:
                    error+=np.abs(testdist[i]-actualdist[i])/maxval
                    num+=1.
            avgerror=error/num*100.

            #error=np.sum(safedivide(np.abs(testdist-actualdist),actualdist))/len(self.isorange)
            print(avgerror,error,num)

        plt.show()

    def TestFinal(self):
        self.numtest=50
        self.testmassrange=np.sort(10**(np.random.uniform(np.log10(min(self.massrange)),np.log10(max(self.massrange)),self.numtest)))
        print(self.testmassrange)
        plt.figure()
        toterror=0.
        errors=[]
        for i in range(0,self.numtest):
            testdist=np.array(finalmodel(self.isorange,self.testmassrange[i]))
            testintensity=np.zeros(self.maxiso)
            output=self.makeisotopes(self.testmassrange[i])
            for o in output:
                testintensity[o[0]]=o[1]
            actualdist=testintensity
            testdist=testdist*np.amax(actualdist)/np.amax(testdist)
            #plt.figure()
            plt.plot(self.isorange,testdist,color="r")
            plt.plot(self.isorange,actualdist,color="b")
            #plt.show()

            error=0.
            maxval=max(actualdist)
            num=0.
            for i in range(0,len(self.isorange)):
                if actualdist[i]>maxval*0.001 and actualdist[i]!=0:
                    error+=np.abs(testdist[i]-actualdist[i])/maxval
                    num+=1.
            avgerror=error/num*100.

            #error=np.sum(safedivide(np.abs(testdist-actualdist),actualdist))/len(self.isorange)
            print(avgerror,error,num)
            toterror+=avgerror
            errors.append(avgerror)
        print("Avg Error: ",toterror/self.numtest)
        print("Max Error: ",np.amax(errors))
        print("Dev: ",np.std(errors))
        plt.show()

iso=Isotopes()
#iso.LoadFits()
#iso.FitIntensities("plot")
#iso.Plot()
#exit()
#iso.FitFinal()
#iso.Train()
#iso.Test()
iso.TestFinal()




