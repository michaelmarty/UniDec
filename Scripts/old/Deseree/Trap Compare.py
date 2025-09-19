import sys

sys.path.append("C:\\Python\\UniDec")
import fnmatch

topdir = "Z:\Group Share\\Deseree\\Ciara"
#topdir="C:\Xcalibur\Data\Ciara CR"


def get_files(directory, conc, chem,trap):
    string = '*_' + str(conc) + "_" + str(chem) + "*"+str(trap)+"*"+".hdf5"
    print(string)
    path=None


    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            print(file)
            path = os.path.join(directory, file)
    return path


chems = ["AA", "IM", "TEA"]
concs = [20, 30, 40, 50, 60]
concs = [40]
xvals=[0,20,40,60,80,100,120,140,160,180,200]
trap=["HCD","IST"]

colors=["r","b","g"]
styles=["-","--","-.",":",":"]

type=0

for i, conc in enumerate(concs):
    for j, chem in enumerate(chems):
        for k, t in enumerate(trap):
            p = get_files(topdir, conc, chem,t)
            if p is not None:
                if type is 0:
                    ydat=[]
                    for x in range(0,len(xvals)):
                        dataset="ms_dataset/"+str(x)
                        data=get_data(p,dataset,"charge_data")
                        data[:,1]/=np.amax(data[:,1])
                        data[data[:,1]<0.8,1]=0
                        if not isempty(data):
                            xval=data[np.argmax(data[:,1]),0]
                            xval=center_of_mass(data)[0]
                            ydat.append(xval)

                if type is 1:
                    ydat = []
                    data = get_data(p, "peaks", "extracts")
                    if not isempty(data):
                        data = data[:, 0]
                        ydat=data

                print(ydat)
                plt.plot(xvals,ydat,label=chem+" "+str(conc), color=colors[j], linestyle=styles[k])
                #plt.plot(data[:,0],data[:,1],color=colors[j], linestyle=styles[i])
                plt.xlabel("Collision Voltage (V)")
                if type is 1:
                    plt.ylabel("Average Mass")
                elif type is 0:
                    plt.ylabel("Peak Charge")
                #plt.legend()
plt.show()
