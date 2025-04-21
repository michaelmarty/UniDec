import numpy as np
import unidec.tools as ud
import matplotlib.pyplot as plt

path = "C:\\Users\\Michael Marty\\OneDrive - University of Arizona\\Teaching\\325 Fall 2019\\Data\\20190328 JGR DESALT COLUMN OX1.raw"

outfile = "C:\\Users\\Michael Marty\\OneDrive - University of Arizona\\Teaching\\325 Fall 2019\\Data\\lcms.csv"
outfile2 = "C:\\Users\\Michael Marty\\OneDrive - University of Arizona\\Teaching\\325 Fall 2019\\Data\\lcms.xlsx"

d = WatersDataImporter(path)
scans = d.get_scans_from_times([0,20])
print(scans)

darray = []
tarray = []
template = d.data[0]
for i, dat in enumerate(d.data):
    if i < scans[1]:
        time = d.times[i]
        tarray.append(time)
        newdat = ud.mergedata(template, d.data[i])
        darray.append(newdat[:, 1])
tarray = np.array(tarray)
darray = np.array(darray)

tarray2 = []
darray2 = []
skip=10
for i in range(0, len(tarray), skip):
    t = np.mean(tarray[i:i+skip])
    d = np.sum(darray[i:i+skip], axis=0)
    tarray2.append(t)
    darray2.append(d)

tarray = np.array(tarray2)
darray = np.array(darray2)
mzaxis = np.array(template[:, 0])

plt.plot(mzaxis, np.sum(darray, axis=0))
plt.show()
#exit()

print(darray.shape, mzaxis.shape, tarray.shape)
tarray = np.insert(tarray, 0, 0)
print(tarray.shape)
outarray = np.vstack(([mzaxis], darray))
print(outarray.shape)
tarray2=np.array([tarray]).transpose()
print(tarray2.shape)
outarray = np.concatenate((tarray2, outarray), axis=1)
print(outarray.shape)

#df = pd.DataFrame(outarray)
#df.to_excel(outfile2)

#outarray.tofile(outfile, sep=",")
np.savetxt(outfile, outarray, delimiter=",")
