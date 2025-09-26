import h5py
from UniDec.unidec_modules.hdf5_tools import get_dataset
import matplotlib.pyplot as plt

file = "Z:/Group Share/Matt/tmp.hdf5"
hdf = h5py.File(file)
msdata = hdf.get("ms_dataset/10")
rawdata = get_dataset(msdata, "raw_data")
data2 = get_dataset(msdata, "processed_data")
massdat = get_dataset(msdata, "mass_data")
compiledstuff = hdf.get("ms_dataset")
massaxis = get_dataset(compiledstuff, "mass_axis")
masssum = get_dataset(compiledstuff, "mass_sum")
peakstuff = hdf.get("peaks")
extracts = get_dataset(peakstuff, "extracts")
massgrid1 = get_dataset(msdata, "mass_grid")
massgrid2 = get_dataset(compiledstuff, "mass_grid")

font = {'fontname' : 'Times New Roman',
        'size'   : 60}

plt.plot(rawdata[:, 0], rawdata[:, 1], 'k', lw=5)
plt.ylabel("Intensity", **font)
plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
plt.xticks([2000, 2500, 3000, 3500, 4000], **font)
plt.yticks([0, 1300000, 2500000], **font)
plt.xlim(2000,4000)
plt.show()
plt.plot(data2[:, 0], data2[:, 1], 'k', lw=5)
plt.ylabel("Intensity", **font)
plt.gca().get_yaxis().get_major_formatter().set_powerlimits((0, 0))
plt.xticks([2000, 2500, 3000, 3500, 4000], **font)
plt.yticks([0, 0.5, 1], **font)
plt.xlim(2000,4000)
plt.show()
masskda = [x/1000 for x in massdat[:, 0]]
plt.plot(masskda, massdat[:, 1], 'k', lw=5)
plt.ylabel("Intensity", **font)
plt.xlim(21, 23)
plt.xticks([21, 22, 23], **font)
plt.yticks([0, 0.5, 1.0, 1.5], **font)
plt.show()
massaxiskda = [x/1000 for x in massaxis]
plt.plot(massaxiskda, masssum, 'k', lw=5)
plt.ylabel("Intensity", **font)
plt.xlim(21, 23)
plt.xticks([21, 22, 23], **font)
plt.yticks([0,0.5,1], **font)
plt.show()
plt.plot(list(range(0, 100, 5)), extracts[:, 0], 'k--', list(range(0, 100, 5)), extracts[:, 1], 'k-', lw=5)
plt.xticks([0,20,40,60,80], **font)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0], **font)
plt.ylabel("Intensity", **font)
plt.show()