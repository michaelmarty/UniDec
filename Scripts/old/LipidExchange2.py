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
    integral=np.trapz(zvals,x=xvals)
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
        integral+=np.trapz(cut[:,1],x=cut[:,0])
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
            plt.show()
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

    date="180215"
    dps=[788,793]
    ps=[759.9,762.5]
    pc=[761.0,765.0,783,787.5]
    dpc=[789.5,797,811.6,819]
    ser=""
    imps=[100,125]
    impc=[50,75]

    #'''
    date="200215"
    dps=[788,795.35]
    ps=[761,765]
    pc=[760,764,782,786]
    dpc=[787,795,809,817]
    imps=[90,150]
    impc=[85,135]
    ser="S2_"
    '''
    #'''
    path="C:\\MassLynx\\Mike.PRO\\Data\\"+date
    os.chdir(path)

    experiments=['NOCA_','CA_']
    samples=['CLIP_','HIS_']
    pols=['POS_','NEG_']
    times=np.arange(0,50,6)

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

    print(times)
    plt.figure(figsize=(10,6))

    collectb=[]
    collecta=[]
    for x in range(0,len(experiments)):
        collect=[]
        for i in range(0,len(samples)):
            series=[]
            series2=[]
            rtime=[]
            d=date+"_"+ser
            for t in range(0,len(times)):
                exper=experiments[x]
                samp=samples[i]
                time=times[t]
                fpos=d+exper+samp+pols[0]+str(time)+"H"
                fneg=d+exper+samp+pols[1]+str(time)+"H"
                #try:
                fin=LipidRatios(fpos,fneg,pc,dpc,ps,dps,impc=impc,imps=imps,imflag=imflag)
                series.append(fin)
                series2.append(fin)
                rtime.append(time)
                print(fpos,exper,samp,time,fin)
                #except:
                #    print "Failed: ",exper,samp,time
                #    series2.append([0,0])
            series=np.array(series)
            plt.plot(rtime,series[:,0],color=color[1],marker=markers[2*x+i],linestyle=style[2*x+i],label=labels3[x]+' '+labels2[i]+' '+labels[1])
            plt.plot(rtime,series[:,1],color=color[0],marker=markers[2*x+i],linestyle=style[2*x+i],label=labels3[x]+' '+labels2[i]+' '+labels[0])
            collect.append(series2)
        if x==0:
            collecta=np.array(collect)
        else:
            collectb=np.array(collect)
        #plt.plot(times,np.sum(collect[:,:,0],axis=0),color="k",linestyle=style[x],label=labels3[x]+" Total PC")
        #plt.plot(times,np.sum(collect[:,:,1],axis=0),color="g",linestyle=style[x],label=labels3[x]+" Total PS")
        #print collect.shape
    ax=plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])

    collect2a=np.array([[collecta[0,t,0],collecta[1,t,0],collecta[0,t,1],collecta[1,t,1]] for t in range(0,len(times))]).transpose()
    collect2b=np.array([[collectb[0,t,0],collectb[1,t,0],collectb[0,t,1],collectb[1,t,1]] for t in range(0,len(times))]).transpose()
    fit=minimize(times,collect2a)
    fit2=minimize(times,collect2b)
    if len(fit)>3:
        thalf2=fit[3]
        thalf2b=fit2[3]
    else:
        thalf2=fit[2]
        thalf2b=fit2[2]

    print("Fit: ",fit)
    a=fit[0]*fit[2]
    b=fit[1]
    k=fit[3]
    print("t infinity",(a+b)/2.)
    print("t half", 1/(k)*1/2.*np.log(2.))
    print("Error: ",error(fit,times,collect2a))

    print("Fit: ",fit2)
    a=fit2[0]*fit2[2]
    b=fit2[1]
    k=fit2[3]
    print("t infinity",(a+b)/2.)
    print("t half", 1/(k)*1/2.*np.log(2.))
    print("Error: ",error(fit2,times,collect2b))

    fitdat=exchange(times,fit[0],fit[1],fit[2],fit[3],thalf2=thalf2)
    fitdatb=exchange(times,fit2[0],fit2[1],fit2[2],fit[3],thalf2=thalf2b)
    #fitdat=error(fit,times,collect2a,outdat="Yes!")
    #fitdatb=error(fit2,times,collect2b,outdat="Yes!")

    plt.legend(bbox_to_anchor=(1.45, 0.5),loc=7,borderaxespad=0.)
    plt.xlabel("Time (hrs)")
    plt.ylabel("Relative Intensity")
    plt.savefig("LipidExchange_"+date+"_"+str(imflag)+".png")

    plt.plot(times,fitdat[0],color='y',label="No Ca Fit")
    plt.plot(times,fitdat[1],color='y')
    plt.plot(times,fitdat[2],color='y')
    plt.plot(times,fitdat[3],color='y')

    plt.plot(times,fitdatb[0],color='g',label="Ca Fit")
    plt.plot(times,fitdatb[1],color='g')
    plt.plot(times,fitdatb[2],color='g')
    plt.plot(times,fitdatb[3],color='g')

    plt.legend(bbox_to_anchor=(1.45, 0.5),loc=7,borderaxespad=0.)
    plt.savefig("LipidExchange_"+date+"_"+str(imflag)+"_fits.png")
    plt.show()