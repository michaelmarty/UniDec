__author__ = 'Michael.Marty'

import numpy as np
import matplotlib.pyplot as plt
import pymzml
import os

ps=[760.,765.]
dps=[788.,794.]
pc=[804.,809.]
dpc=[832.,839.]

stds=[0,5,10,20,35,50]

ranges=[ps,dps,pc,dpc]

mids=[(r[0]+r[1])/2. for r in ranges]
spans=[(r[1]-r[0])/2. for r in ranges]



header="MTM150515_PSPC_"
dir="C:\\MassLynx\\Mike.PRO\\Data\\150515\\mzML"

#plt.figure()
hext=[]
aext=[]
uext=[]
for s in stds:
    file=header+str(s)+".mzML"


    path=os.path.join(dir,file)
    msrun=pymzml.run.Reader(path)
    times=[]
    cints=[]
    plt.figure()
    for spectrum in msrun:
        try:
            scantime=float(spectrum["scan time"])
            #print scantime,time,time+timestep
        except:
            scantime=-1

        if scantime>14 and scantime<22:
            times.append(scantime)
            #impdat=np.transpose([spectrum.mz,spectrum.i])
            mz=np.array(spectrum.mz)
            intensity=np.array(spectrum.i)

            ints=[]
            for i in range(0,len(ranges)):
                mid=mids[i]
                span=spans[i]
                mztest=np.abs(mz-mid)
                btest=mztest<span
                inttot=np.sum(intensity[btest])
                ints.append(inttot)
                #plt.plot(mz[btest],intensity[btest],color="b")
            cints.append(ints)

    cints=np.transpose(cints)
    times=np.array(times)

    maxes=np.amax(cints,axis=1)
    psr=maxes[0]/maxes[1]
    pcr=maxes[2]/maxes[3]
    #print psr,pcr
    tot=psr+pcr
    psp=psr/tot*100
    pcp=pcr/tot*100
    print("POPS: %.2f  POPC: %.2f " %(psp,pcp))
    hext.append([psp,pcp])

    maxes=np.sum(cints,axis=1)
    psr=maxes[0]/maxes[1]
    pcr=maxes[2]/maxes[3]
    #print psr,pcr
    tot=psr+pcr
    psp=psr/tot*100
    pcp=pcr/tot*100
    print("POPS: %.2f  POPC: %.2f " %(psp,pcp))

    aext.append([psp,pcp])

    maxes=np.sum(cints,axis=1)
    print(maxes[1]/maxes[3])
    psr=maxes[0]
    pcr=maxes[2]
    #print psr,pcr
    tot=psr+pcr
    psp=psr/tot*100
    pcp=pcr/tot*100
    print("POPS: %.2f  POPC: %.2f " %(psp,pcp))

    uext.append([psp,pcp])

    for i in range(0,len(ranges)):
        plt.plot(times,cints[i])
plt.gca().set_xlabel("Time (min)")
plt.show()

hext=np.array(hext)
aext=np.array(aext)
uext=np.array(uext)

dext=np.array([0.0075,0.10655,0.223953,0.3602313,0.5521191,0.647349])*100

da=np.average(dext[1:]/stds[1:])
ha=np.average(hext[1:,0]/stds[1:])
aa=np.average(aext[1:,0]/stds[1:])

print(da,ha,aa)
print(np.average([da,ha,aa]))

plt.figure()
plt.plot(stds,hext[:,0],label="LCMS Height",marker="o")
plt.plot(stds,aext[:,0],label="LCMS Area",marker="o")
plt.plot(stds,uext[:,0],label="LCMS Absolute Area",marker="o")
plt.plot(stds,dext,label="Direct MS",marker="o")
plt.plot(stds,stds,label="Ideal",linestyle="--")
#plt.plot(stds,np.array(stds),label="Ideal",linestyle="--")
plt.legend(loc=4)
plt.xlabel("%PS Standard")
plt.ylabel("%PS Measured")
plt.show()