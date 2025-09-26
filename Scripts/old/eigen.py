__author__ = 'Michael.Marty'

import numpy as np
#from copy import deepcopy
#import random
import sklearn
import matplotlib.cm as cm
import wx
import time
from ctypes import *
from UniDec.unidec_modules import PlottingWindow
from UniDec.unidec_modules import unidectools as ud

try:
    libs=cdll.LoadLibrary("libmypfunc.dll")
except:
    print("Failed to load libmypfunc, components might be slow")

def ndis(x, y, s):
    return np.exp(-(x - y) * (x - y) / (2.0 * s * s))


class EigenWindow(wx.Dialog):
    def __init__(self, *args, **kwargs):
        super(EigenWindow, self).__init__(*args, **kwargs)
        self.SetSize((1400, 900))
        self.SetTitle("Components Tools")

    def InitUI(self,*args,**kwargs):
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.plot1= PlottingWindow.Plot1d(self, figsize=(10, 8))

        self.vbox.Add(self.plot1,1,wx.EXPAND)

        self.SetSizer(self.vbox)
        self.Center()

    def Eigen(self,massdat,ztab,massgrid,config,pks):
        self.eig=Eigen(massdat,ztab,massgrid,config,pks,plot=self.plot1)

    def Cluster(self):
        self.eig.Cluster()

class Eigen:
    def __init__(self,massdat,ztab,massgrid,config,pks,plot,*args,**kwargs):
        
        print("Determining Clusters")
        tstart=time.perf_counter()
        self.peakwidth=config.mzsig
        self.massdat=massdat
        self.massdat[:,1]=self.massdat[:,1]/np.amax(self.massdat[:,1])
        self.psfun=config.psfun
        self.cutoff=config.peakthresh
        self.plot=plot
        self.pks=pks
        self.mgrid,self.zgrid = np.meshgrid(massdat[:,0],ztab,indexing='ij')
        self.mzgrid=(self.mgrid+config.adductmass*self.zgrid)/self.zgrid
        self.massgrid=np.reshape(massgrid,self.mzgrid.shape)

        self.blackouts=self.SetBlackouts()

        self.totalnum=10
        #print self.blackouts
        try:
            self.maskresults=self.CMask()
        except:
            print("Failed C speedup")
            self.maskresults=np.array([self.InverseMaskContinuous(self.blackouts[i%len(self.blackouts)]) for i in range(0,len(self.blackouts)*self.totalnum)])
        tend=time.perf_counter()
        print("Masking: %.2gs" % (tend - tstart))
        scaler=config.resolution
        self.totals=np.array([[massdat[np.argmax(x),0],1/scaler*np.amax(x)*(np.amax(massdat[:,0])-np.amin(massdat[:,0]))+np.amin(massdat[:,0])] for x in self.maskresults])
        self.totals2=np.array([[massdat[np.argmax(x),0],np.amax(x)] for x in self.maskresults])


        bool=self.totals2[:,1]>0*np.amax(self.totals2[:,1])
        self.subresults=self.maskresults[bool]
        self.totals=self.totals[bool]
        self.totals2=self.totals2[bool]

        ax1=self.plot.figure.add_subplot(211)
        ax2=self.plot.figure.add_subplot(212,sharex=ax1)
        ax2.scatter(self.totals2[:,0],self.totals2[:,1])
        self.plot.repaint()
        tend=time.perf_counter()
        print("Shadowing done: %.2gs" % (tend - tstart))

    def CMask(self):

        width=self.peakwidth*1
        rands=np.random.uniform(-width,width,len(self.blackouts)*self.totalnum)
        brand=np.array([self.blackouts[i%len(self.blackouts)]+rands[i] for i in range(0,len(rands))])

        bnum=len(brand)
        bnum2=bnum*len(self.massdat)

        mzgrid=np.ravel(self.mzgrid)
        massgrid=np.ravel(self.massgrid)

        cmzgrid=(c_double*len(mzgrid))(*mzgrid)
        cmassgrid=(c_double*len(massgrid))(*massgrid)
        crand=(c_double*bnum)(*brand)
        cout=(c_double*bnum2)()

        gdim1=len(self.massdat)
        gdim2=len(self.mzgrid[0])
        psfun=int(self.psfun)
        libs.inversemask(byref(cmzgrid),
                         byref(cmassgrid),
                         byref(crand),
                         byref(cout),
                         (c_double)(width),
                         (c_int)(psfun),
                         (c_int)(bnum),
                         (c_int)(gdim1),
                        (c_int)(gdim2))

        output=np.frombuffer(cout)
        return output.reshape((bnum,len(self.massdat)))

    def Cluster(self):
        tstart=time.perf_counter()
        bandwidth=sklearn.cluster.estimate_bandwidth(self.totals,quantile=0.1)
        cluster=sklearn.cluster.mean_shift(self.totals,bandwidth=bandwidth,bin_seeding=True,min_bin_freq=int(max(1,self.totalnum/3)),cluster_all=False)
        nclust=len(cluster[0])
        #cluster=sklearn.cluster.AgglomerativeClustering(nclust,linkage="ward").fit_predict(self.totals)
        #cluster=[np.unique(cluster),cluster]
        #nclust=len(cluster[0])
        #cluster=sklearn.cluster.k_means(self.totals,nclust)
        print("Clusters: ",cluster[0])
        print("Number of clusters: ",nclust)

        colors = cm.rainbow(np.linspace(0, 1, nclust))
        self.plot.clear_plot()
        ax1=self.plot.figure.add_subplot(211)
        ax2=self.plot.figure.add_subplot(212,sharex=ax1)
        ax1.plot(self.massdat[:,0],self.massdat[:,1],color="k")
        maxval=np.amax([np.sum(self.subresults[cluster[1]==i],axis=0) for i in np.unique(cluster[1])])
        for i in np.unique(cluster[1]):
            if i!=-1:
                subtotals=self.totals2[cluster[1]==i]
                subspectrum=np.sum(self.subresults[cluster[1]==i],axis=0)
                scale=1/maxval*np.amax(self.massdat[:,1])
                ax1.plot(self.massdat[:,0],subspectrum*scale,color=colors[i])
                ax2.scatter(subtotals[:,0],subtotals[:,1],color=colors[i])


        self.plot.setup_zoom([ax1,ax2], 'box')

        self.plot.repaint()
        tend=time.perf_counter()
        print("Done clustering: %.2gs" % (tend - tstart))

    def SetBlackouts(self):
        boo1=self.massgrid>0
        usedmz=self.mzgrid[boo1]
        mzrange=np.arange(np.amin(usedmz),np.amax(usedmz),self.peakwidth)
        return mzrange

    def InverseMaskContinuous(self,blackout):
        width=self.peakwidth*1
        blackout=blackout+np.random.uniform(-width,width)
        bool=np.abs(self.mzgrid-blackout)<width*5
        #prob=ndis(self.mzgrid[bool],blackout,self.peakwidth*1)
        prob=ud.make_peak_shape(self.mzgrid[bool],self.psfun,width,blackout)
        grid=np.zeros_like(self.massgrid)
        grid[bool]=self.massgrid[bool]*prob
        out=np.sum(grid,axis=1)
        return out