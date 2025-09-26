import sys

sys.path.append("C:\\Python\\UniDec")
import fnmatch

topdir = "Z:\Group Share\Deseree\Ciara\N2 and Ar Runs for Comparison"
#topdir="C:\Xcalibur\Data\Ciara CR"


def get_files(directory, conc, chem, gas):
    string = '*_' + str(conc) + "_" + str(chem) + "*"+"IST"+"_" + str(gas) + "*"+".hdf5"
    print(string)
    path=None
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            print(file)
            path = os.path.join(directory, file)
    return path


chems = ["AA", "TEA"]
concs = [40]
gases = ["AR","N"]
xvals=[0,20,40,60,80,100,120,140,160,180,200]

colors=["r","b"]
styles=["-","--","-.",":",":"]

for i, conc in enumerate(concs):
    for j, chem in enumerate(chems):
        for k, gas in enumerate(gases):
            p = get_files(topdir, conc, chem, gas)
            if p is not None:
                #data=get_data(p,"ms_dataset/1","charge_data")
                data=get_data(p,"peaks","extracts")
                if not isempty(data):
                    data=data[:,0]
                    print(data)
                    plt.plot(xvals,data,label=chem+" "+str(conc)+" "+str(gas), color=colors[k], linestyle=styles[j])
                    #plt.plot(data[:,0],data[:,1],color=colors[j], linestyle=styles[i])
                    plt.xlabel("Collision Voltage (V)")
                    plt.ylabel("Average Mass (Da)")
                    #plt.legend()
plt.show()
