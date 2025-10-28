import sys

sys.path.append("C:\\Python\\UniDec")
import fnmatch

topdir = "Z:\Group Share\\Deseree\\Ciara"
#topdir="C:\Xcalibur\Data\Ciara CR"


def get_files(directory, conc, chem):
    string = '*_' + str(conc) + "_" + str(chem) + "*"+"IST"+"*"+".hdf5"
    print(string)
    path=None


    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            print(file)
            path = os.path.join(directory, file)
    return path


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
