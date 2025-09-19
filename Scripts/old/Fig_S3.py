__author__ = 'Michael.Marty'

import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
import matplotlib as mpl
from matplotlib import rcParams

from UniDec.unidec_modules.isolated_packages import option_d

rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5
#rcParams['font.size']=18

def mergedata(data1,data2):
    f = interp1d(data2[:, 0], data2[:, 1],bounds_error=False,fill_value=0)
    intx=data1[:,0]
    inty=f(intx)
    newdat = np.column_stack((intx, inty))
    return newdat

topdir="C:\\SVN\\Data"
proteins=["mz pos","mass neg"]
pmass=[126783,98848]
lmass=[760.076,760.076]

datafiles2=[["AmtB_DMPC\\AmtB_09_2","AmtB_DMPC\\AmtB_09_2"],["AmtB_DMPC\\AmtB_09_neg","AmtB_DMPC\\AmtB_09_neg"]]
ncols=len(datafiles2[0])+1
ratios=np.ones(ncols)
ratios[-1]=0.05   ### ratio of column with, 0.05 is the last column for colorbar
ylabels=["Positive Mode\nCollision Voltage (V)","Negative Mode\nCollision Voltage (V)"]

gs=gridspec.GridSpec(2,ncols,width_ratios=ratios)
gs.update(wspace=0.1,hspace=0.1)

fig=plt.figure(figsize=(6.3,5))

#magicnumber=[[12.882,20.045,72.040,54,10.32],[7.498,16.385,78.802,44,5.58]]
#magicnumber=[[10.32,20.045,72.040,54],[5.58,16.385,78.802,44]]
#magicnumber=[[7.498,78.002],[12.882,72.040]]

cmap= option_d.viridis
#cmap="jet"

size=8
for i,protein in enumerate(proteins):
    group=datafiles2[i]
    for j,folder in enumerate(group):

        if j==0:
            maximum=20000
            minorticks=np.arange(0,maximum,5000) ###change later
            majorticks=np.arange(50,190,50)
            #MS PLOTS

            plt.subplot(gs[i*ncols+j])
            ax=plt.gca()
            ax.set_xticks(minorticks)
            ax.set_yticks(majorticks)
            ax.set_xlim(2500,maximum)
            if i==1:
                plt.xlabel("m/z",fontsize=size)
                #plt.text(18000,170,"Negative",color="w",fontsize=12)
            else:
                pass
                #plt.text(18000,170,"Positive",color="w",fontsize=12)
        else:
            maximum=220
            minorticks=np.arange(0,maximum,25) ###change later
            majorticks=np.arange(50,190,50)
            #MS PLOTS

            plt.subplot(gs[i*ncols+j])
            ax=plt.gca()
            ax.set_xticks(minorticks,minor=True)
            ax.set_yticks(majorticks)
            ax.set_xlim(20,maximum)
            if i==1:
                plt.xlabel("Mass (kDa)",fontsize=size)
                #plt.text(200,170,"Negative",color="w",fontsize=12)
            else:
                pass
                #plt.text(200,170,"Positive",color="w",fontsize=12)

        if i<1:
            ax.set_xticklabels([])
        if not j==0:
            ax.set_yticklabels([])
        else:
            #ax.set_yticklabels([0,"",1],fontsize=size)
            plt.ylabel(ylabels[i],fontsize=size)
        plt.setp(ax.get_yticklabels(),fontsize=size)
        plt.setp(ax.get_xticklabels(),fontsize=size,color="k")
        #labels=["POPC","POPC E3","DMPC","DMPC E3"]  ###x axis label
        #styles=["-",":","--","-."]
        #if i==0:
        #    plt.title(labels[j],fontsize=size)

        if j==0:
            path=os.path.join(topdir,os.path.join(folder,"ExtractFull2D_Th_.txt"))
            print(path)
            data=np.loadtxt(path)
            data[:,2]=data[:,2]/np.amax(data[:,2])
            X = data[:, 0]    ###m/z or mass
        else:
            path=os.path.join(topdir,os.path.join(folder,"ExtractFull2D_Da_.txt"))
            print(path)
            data=np.loadtxt(path)
            data[:,2]=data[:,2]/np.amax(data[:,2])
            X = data[:, 0]/1000.    ###m/z or mass
        Y = data[:, 1]        ###volt
        if np.max(Y)<100:      ###some of them use 0-19 for CV, not 10-190
            Y = Y*10
        C = data[:, 2]    ###intensity

        xvals = np.unique(X)
        yvals = np.unique(Y)
        xlen = len(xvals)
        ylen = len(yvals)
        newgrid = np.reshape(C, (xlen, ylen))
        print(newgrid.shape)
        newx = np.reshape(X, (xlen, ylen))
        newy = np.reshape(Y, (xlen, ylen))
        cax=plt.contourf(newx, newy, newgrid, 100, cmap=cmap)
        ticcol="w"
        ax=plt.gca()
        for line in ax.xaxis.get_ticklines():
            line.set_color(ticcol)
        for line in ax.yaxis.get_ticklines():
            line.set_color(ticcol)
        #ax.tick_params(axis="x",which="both",colors=ticcol)
        #for v in magicnumber[i]:
        #    pass
        #    plt.vlines(v,0,1)

    plt.subplot(gs[i*ncols+ncols-1])
    ax=plt.gca()
    cbar=mpl.colorbar.ColorbarBase(ax,cmap=cmap,orientation="vertical",ticks=[0,0.5,1])
    #cbar=plt.colorbar(use_gridspec=False,ticks=[0,np.amax(newgrid)/2,np.amax(newgrid)])
    cbar.ax.get_yaxis().set_tick_params(direction='out')
    cbar.ax.set_yticklabels(["0","%","100"],fontsize=size)
    plt.subplots_adjust(left=0.125, right=0.925)

plt.savefig("C:\\Users\\kin.hoi\\Dropbox\\MSNanodisc\\Figures\\Figure_S3.png")
#plt.savefig("C:\\Users\\kin.hoi\\Dropbox\\MSNanodisc\\Figures\\Figure_S2.eps")
plt.show()

