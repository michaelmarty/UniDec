__author__ = 'Michael.Marty'

import numpy as np
import os
import matplotlib.pyplot as plt
from LipidExchange2 import *




path="C:\\MassLynx\\Mike.PRO\\Data\\130215"
os.chdir(path)

date="130215_"
experiments=['NOCA_','CA_']
samples=['CLIP_','HIS_']
pols=['POS_','NEG_']
times=np.arange(0,50,24)

dps=[787,795]
ps=[761.0,764.0]
pc=[795,800]
dpc=[823,831]


labels=["POPS","POPC"]
labels2=["Clip","His"]
labels3=["No Ca","Ca"]
color=["r","b"]
markers=["o","s",">","v"]
style=["-","--","-.",":"]

print(times)
plt.figure(figsize=(10,6))
for x in range(0,len(experiments)):
    collect=[]
    for i in range(0,len(samples)):
        series=[]
        series2=[]
        rtime=[]
        for t in range(0,len(times)):
            exper=experiments[x]
            samp=samples[i]
            time=times[t]
            if 1==1:
                #posdata=np.loadtxt(date+exper+samp+pols[0]+str(time)+"H_rawdata.txt")
                negdata=np.loadtxt(date+exper+samp+str(time)+"H_rawdata.txt")

                #posdata[:,1]=posdata[:,1]/np.amax(posdata[:,1])
                negdata[:,1]=negdata[:,1]/np.amax(negdata[:,1])

                pcext=extract(negdata,pc)
                dpcext=extract(negdata,dpc)
                pcrat=pcext/dpcext

                psext=extract(negdata,ps)
                dpsext=extract(negdata,dps)
                psrat=psext/dpsext

                tot=pcrat+psrat
                fin=[pcrat/tot,psrat/tot]
                #collect[i,t]=fin
                series.append(fin)
                series2.append(fin)
                rtime.append(time)
                print(exper,samp,time,fin)
            #except:
                #print "Failed: ",exper,samp,time
                #series2.append([0,0])
        series=np.array(series)
        plt.plot(rtime,series[:,0],color=color[1],marker=markers[2*x+i],linestyle=style[2*x+i],label=labels3[x]+' '+labels2[i]+' '+labels[1])
        plt.plot(rtime,series[:,1],color=color[0],marker=markers[2*x+i],linestyle=style[2*x+i],label=labels3[x]+' '+labels2[i]+' '+labels[0])
        collect.append(series2)
    collect=np.array(collect)
    #plt.plot(times,np.sum(collect[:,:,0],axis=0),color="k",linestyle=style[x],label=labels3[x]+" Total PC")
    #plt.plot(times,np.sum(collect[:,:,1],axis=0),color="g",linestyle=style[x],label=labels3[x]+" Total PS")
ax=plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

plt.legend(bbox_to_anchor=(1.3, 0.5),loc=7,borderaxespad=0.)
plt.xlabel("Time (hrs)")
plt.ylabel("Relative Intensity")


fit=exchange(times,0.4,0.5,0.05)
plt.plot(times,fit[0])
plt.plot(times,fit[1])
plt.plot(times,fit[2])
plt.plot(times,fit[3])

plt.show()