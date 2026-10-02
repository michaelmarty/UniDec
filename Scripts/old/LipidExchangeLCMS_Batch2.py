__author__ = 'Michael.Marty'

import numpy as np
import matplotlib.pyplot as plt
import pymzml
import os

ps=[760.,765.]
dps=[788.,794.]
pc=[804.,809.]
dpc=[832.,839.]

ps=[760.5129,761.5152,762.5176,763.52,764.5223]

pc=[804.5755,805.5778,806.5802,807.5826,808.5849]

dps=[788.68,789.69,790.70,791.70,792.71,793.71,794.71]

dpc=[832.75,833.75,834.76,835.77,836.77,837.77,838.77]


peaks=[ps,dps,pc,dpc]
sigma=0.02


exchangetimes=[0,4,8,12,18,24,36,48]
samples=["++","+-","-+","--"]

header="MTM_150520_"
dir="C:\\MassLynx\\Mike.PRO\\Data\\150520\\mzml"

#plt.figure()
aext=[]
for t in exchangetimes:
    taext=[]
    for s in samples:
        file=header+str(t)+s+".mzML"
        path=os.path.join(dir,file)
        msrun=pymzml.run.Reader(path)
        times=[]
        cints=[]
        #plt.figure()
        for spectrum in msrun:
            try:
                scantime=float(spectrum["scan time"])
                #print scantime,time,time+timestep
            except:
                scantime=-1

            if scantime>14 and scantime<30:
                times.append(scantime)
                #impdat=np.transpose([spectrum.mz,spectrum.i])
                mz=np.array(spectrum.mz)
                intensity=np.array(spectrum.i)

                ints=[]
                for i in range(0,len(peaks)):
                    pcol=peaks[i]
                    inttot=0
                    for p in pcol:
                        mztest=np.abs(mz-p)
                        btest=mztest<sigma
                        inttot+=np.trapezoid(intensity[btest],mz[btest])
                        #plt.plot(mz[btest],intensity[btest])
                    ints.append(inttot)


                cints.append(ints)
        #plt.show()

        cints=np.transpose(cints)
        times=np.array(times)

        maxes=np.trapezoid(cints,times,axis=1)
        psr=maxes[0]/maxes[1]
        pcr=maxes[2]/maxes[3]
        #print psr,pcr
        tot=psr+pcr
        psp=psr/tot*100
        pcp=pcr/tot*100
        print("POPS: %.2f  POPC: %.2f " %(psp,pcp))
        taext.append([psp,pcp])

        #for i in xrange(0,len(samples)):
        #    plt.plot(times,cints[i])
    aext.append(taext)
    #plt.show()

aext=np.array(aext)
aext=np.transpose(aext,axes=(1,0,2))

markers=["o","s",">","v"]
style=["-","--","-.",":"]
plt.figure()
for i,s in enumerate(samples):
    plt.plot(exchangetimes,aext[i,:,0],label="PS"+s,marker=markers[i],linestyle=style[i],color="r")
    plt.plot(exchangetimes,aext[i,:,1],label="PC"+s,marker=markers[i],linestyle=style[i],color="b")
ax=plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
#plt.plot(stds,np.array(stds),label="Ideal",linestyle="--")
plt.legend(bbox_to_anchor=(1.45, 0.5),loc=7,borderaxespad=0.)
plt.xlabel("Time (hr)")
plt.ylabel("Relative Intensity")
plt.show()