__author__ = 'Michael.Marty'

import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math

def match(list,x,y):
    bool1=list[:,0]==x
    bool2=list[:,1]==y
    bool3=np.all(np.array([bool1,bool2]),axis=0)
    cut=list[bool3]
    val=np.sum(cut[:,2])
    return val

def MakeGrid(list,xvals,yvals):
    grid=np.zeros((len(xvals),len(yvals)))
    for i in list:
        xpos=np.searchsorted(xvals,i[0])
        ypos=np.searchsorted(yvals,i[1])
        grid[xpos,ypos]=i[2]
    return grid
    #return np.array([[match(list,x,y) for y in yvals]for x in xvals])


def plot2d(list):
    xvals=np.unique(list[:,0])
    yvals=np.unique(list[:,1])
    Z=MakeGrid(list,xvals,yvals)

    zvals=np.sum(Z,axis=1)
    integral=np.trapezoid(zvals,x=xvals)
    #print integral,np.sum(Z)/integral
    #print zvals.shape
    #plt.figure()
    #plt.contourf(xvals,yvals,Z.transpose(),100)
    #plt.show()
    return integral

def extract(data,ranges):
    length=len(ranges)/2
    integral=0
    booarr=[]
    for i in range(0,length):
        min=ranges[i*2]
        max=ranges[i*2+1]
        bool1=data[:,0]<max
        bool2=data[:,0]>min
        bool3=np.all(np.array([bool1,bool2]),axis=0)
        cut=data[bool3]
        integral+=np.trapezoid(cut[:,1],x=cut[:,0])
        booarr.append(bool3)
    #bool4=np.any(np.array(booarr),axis=0)
    #cut=data[bool4]
    return integral

def imextract(data,ranges,imranges):
    length=len(ranges)/2
    integral=0
    booarr=[]
    for i in range(0,length):
        min=ranges[i*2]
        max=ranges[i*2+1]
        immin=imranges[0]
        immax=imranges[1]
        bool1=data[:,0]<max
        bool2=data[:,0]>min
        bool3=data[:,1]<immax
        bool4=data[:,1]>immin
        bool5=np.all(np.array([bool1,bool2,bool3,bool4]),axis=0)
        cut=data[bool5]
        integral+=plot2d(cut)#np.sum(cut[:,2])#plot2d(cut)
        booarr.append(bool5)
    bool6=np.any(np.array(booarr),axis=0)
    cut=data[bool6]
    #integral=plot2d(cut)
    return integral

def model(a,b,k,t):
    #y=0.5*np.exp(-2.0*k*t)*(a-b+(a+b)*np.exp(2.0*k*t))
    #z=0.5*np.exp(-2.0*k*t)*(b-a+(a+b)*np.exp(2.0*k*t))
    y=(a-b)/2.*np.exp(-2.*k*t)+(a+b)/2.
    z=(b-a)/2.*np.exp(-2.*k*t)+(a+b)/2.
    return [y,z]

def exchange(times,initpspc,initpspc0,alpha,thalf,thalf2=None):
    if thalf2==None:
        thalf2=thalf
    outps=np.array([model(initpspc0,initpspc*alpha,thalf,t) for t in times])
    outpc=np.array([model(1.-initpspc0,(1.-initpspc)*alpha,thalf2,t) for t in times])
    norm0=outpc[:,0]+outps[:,0]
    norm1=outpc[:,1]+outps[:,1]
    return np.array([outpc[:,0]/norm0,outpc[:,1]/norm1,outps[:,0]/norm0,outps[:,1]/norm1])

def error(arr,times,collect,outdat=None):
    if len(arr)>4:
        thalf2=arr[4]
    else:
        thalf2=None
    fit=np.ravel(exchange(times,arr[0],arr[1],arr[2],arr[3],thalf2))
    dat=np.ravel(collect)
    bool1=dat>0
    if outdat is None:
        return np.sum((fit[bool1]-dat[bool1])**2)
    else:
        return fit

def minimize(times,collect):
    fit=opt.minimize(error,[0.5,0.0,1.0,0.05],args=(times,collect),method="Nelder-Mead").x
    return fit

def LipidRatios(fpos,fneg,pc,dpc,ps,dps,impc=None,imps=None,imflag=1,plot=None):
    if imflag==0:

        posdata=np.loadtxt(fpos+"_rawdata.txt")
        negdata=np.loadtxt(fneg+"_rawdata.txt")
        if plot is not None:
            plt.figure()
            plt.plot(posdata[:,0],posdata[:,1]/np.amax(posdata[:,1]))
            plt.plot(negdata[:,0],negdata[:,1]/np.amax(negdata[:,1]))
            #plt.show()
        pcext=extract(posdata,pc)
        dpcext=extract(posdata,dpc)
        psext=extract(negdata,ps)
        dpsext=extract(negdata,dps)

    else:
        posdata=np.loadtxt(fpos+"_im.txt")
        negdata=np.loadtxt(fneg+"_im.txt")


        if impc==None or imps==None:
            print("Need to define IM range")
            exit()
        pcext=imextract(posdata,pc,impc)
        dpcext=imextract(posdata,dpc,impc)
        psext=imextract(negdata,ps,imps)
        dpsext=imextract(negdata,dps,imps)

    pcrat=pcext/dpcext
    psrat=psext/dpsext
    tot=pcrat+psrat
    #print pcext,dpcext,pcrat,psext,dpsext,psrat
    fin=[pcrat/tot,psrat/tot]
    return fin

def FlipPol(fname):

    if "POS" in fname:
        old="POS"
        new="NEG"
        posname=fname
        negname=str(fname).replace(old,new)
    elif "NEG" in fname:
        old="NEG"
        new="POS"
        posname=str(fname).replace(old,new)
        negname=fname
    return [posname,negname]

'''
date="190215"
path="C:\\MassLynx\\Mike.PRO\\Data\\"+date
os.chdir(path)
tail="_QEX_S2_CA_CLIP_POS_48H_1_10.raw"
file=date+tail
fname=file[:-4]

names=FlipPol(fname)

imflag=0
if imflag==1:
    os.chdir("IMMS")
    print os.getcwd()

print LipidRatios(names[0],names[1],pc,dpc,ps,dps,impc=[95,125],imps=[75,95],imflag=imflag)
'''


if __name__=="__main__":

    date="270415"
    dps=[788,795.35]
    ps=[761,765]
    pc=[760,764,782,786]
    dpc=[787,795,809,817]
    imps=[90,150]
    impc=[85,135]
    imps=[60,120]
    impc=[60,100]

    path="C:\\MassLynx\\Mike.PRO\\Data\\"+date+"_stds"
    os.chdir(path)

    pols=['_POS','_NEG']
    conc=np.array([0,5,10,20,35,50])

    imflag=1
    if imflag==1:
        os.chdir("IMMS")
        print(os.getcwd())

    labels=["POPS","POPC"]
    labels2=["Clip","His"]
    labels3=["No Ca","Ca"]
    color=["r","b"]
    markers=["o","s",">","v"]
    style=["-","--","-.",":"]


    series=[]
    rconc=[]
    for t in range(0,len(conc)):

        c=conc[t]
        fpos="STD_"+str(c)+pols[0]
        fneg="STD_"+str(c)+pols[1]

        fin=LipidRatios(fpos,fneg,pc,dpc,ps,dps,impc=impc,imps=imps,imflag=imflag,plot=True)
        series.append(fin)

        rconc.append(c)
        print(fpos,c,fin)

    series=np.array(series)
    plt.figure(figsize=(10,6))
    plt.plot(rconc,series[:,0],color=color[1],marker=markers[0])
    plt.plot(rconc,1-np.array(rconc)*0.01,color=color[1],linestyle=style[1])
    plt.plot(rconc,series[:,1],color=color[0],marker=markers[0])
    plt.plot(rconc,np.array(rconc)*0.01,color=color[0],linestyle=style[1])
    plt.show()

