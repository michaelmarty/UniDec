import sys
import os
sys.path.append("C:\\Python\\UniDec")
import fnmatch

from matplotlib import rcParams, pyplot as plt

#rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth']=0.5
rcParams['axes.linewidth']=0.5
#rcParams['font.size']=18
# matplotlib.rc('font', family='sans-serif')
# matplotlib.rc('font', serif='Helvetica')

def get_files(directory, lip):
    string = '*BL_' + str(lip) + "_*"+".hdf5"
    print(string)
    path=None

    paths= []

    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            print(file)
            path = os.path.join(directory, file)
            paths.append(path)
    return paths

topdir = "Z:\wresager\\Bleaching\\New Data"
lip = "cont"
lmass= 750.
basemass=41970.
paths = get_files(topdir,lip)
os.chdir(topdir)
ptargets=[41980, 41740, 42143]

window=20
exglob=[]
fitglob=[]
for n,file in enumerate(paths):
    print(file)
    peaks=get_data(file,"peaks","peakdata")
    extracts=get_data(file,"peaks","extracts")
    times=get_metadata(file,"timemid")


    fint=[]
    exdata=[]
    for i,t in enumerate(ptargets):
        for j,p in enumerate(peaks[:,0]):
            if np.abs(t-p)<window:
                ex=extracts[:,j]
                exdata.append(ex)
                fit, fitdat = exp_fit(times, ex, bguess=0)
                print(t,p,fit)
                fint.append(fit)
                plt.show()
    fint=np.array(fint)
    print(fint)
    fint/=fint[0,0]
    print(fint)

    fitglob.append(fint)
    exglob.append(exdata)
    print(np.shape(exdata))
exglob=np.array(exglob)
fitglob=np.array(fitglob)
print(fitglob.shape)
print(fitglob)
kavg=np.average(fitglob, axis=0)
kstd=np.std(fitglob, axis=0)
print(kavg.shape)
print(kavg[:,0])
print(kstd[:,0])
peak= ("Rho", "-Palm", "+Mannose")
colors=["b","orange","purple"]
y_pos = np.arange(len(peak))
plt.bar(y_pos, kavg[:,0], align='center', alpha=1.0, color=colors, yerr=kstd[:,0])
plt.xticks(y_pos, peak)
plt.title('Normalized Rate of Aggregation')
plt.ylabel("Normalized K value")
plt.show()
avg=np.average(exglob, axis=0)
std=np.std(exglob,axis=0)
marker=["o","s","<"]
print(np.shape(avg))
for i,a in enumerate(avg):
    fit, fitdat=exp_fit(times,a,bguess=0)

    plt.plot(times,fitdat, color=colors[i])
    plt.scatter(times,a, color=colors[i],marker=marker[i])
    plt.title('Peak Decay Over Time')
    plt.xlabel("Time(min)")
    plt.ylabel("Intesity")
#plt.gca().set_xlim([1,5])
#plt.gca().set_ylim([0,1])
plt.show()




exit()
plt.errorbar(times, a, marker="o", linestyle=" ", yerr=std[i])
chems = ["AA", "IM", "TEA"]
concs = [20, 30, 40, 50, 60]
concs = [20,40,60]
xvals=[0,20,40,60,80,100,120,140,160,180,200]

colors=["r","b","g"]
styles=["-","--","-.",":",":"]

for i, conc in enumerate(concs):
    for j, chem in enumerate(chems):
        p = get_files(topdir, conc, chem)
        if p is not None:
            cdat=[]
            for x in range(0,len(xvals)):
                dataset="ms_dataset/"+str(x)
                data=get_data(p,dataset,"charge_data")
                data[:,1]/=np.amax(data[:,1])
                data[data[:,1]<0.1,1]=0
                if not isempty(data):
                    xval=data[np.argmax(data[:,1]),0]
                    xval=center_of_mass(data)[0]
                    cdat.append(xval)
            print(cdat)
            plt.plot(xvals,cdat,label=chem+" "+str(conc), color=colors[j], linestyle=styles[i])
            #plt.plot(data[:,0],data[:,1],color=colors[j], linestyle=styles[i])
            plt.xlabel("Collision Voltage (V)")
            plt.ylabel("Peak Charge")
            #plt.legend()
plt.show()
