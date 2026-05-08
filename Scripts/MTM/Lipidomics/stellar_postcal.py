import os
import numpy as np
import unidec.UniDecImporter.ImporterFactory as imp
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.fft as fft


def linearize(d, min, binsize, newx):
    newdata = np.zeros_like(newx)
    indexes = d[:,0] - min
    # Round to the nearest binsize
    indexes = np.round(indexes / binsize).astype(int)
    newdata[indexes] = d[:,1]
    return newdata


def apply_linearize(data, minx=100, maxx=2001, binsize=0.1, normalize=True):
    newx = np.arange(minx, maxx, binsize)
    outdata = np.zeros((len(data), len(newx)))
    for i, d in enumerate(data):
        d2 = linearize(d, minx, binsize, newx)
        if np.amax(d2) == 0:
            continue
        if normalize:
            d2 = d2/np.amax(d2)
        outdata[i, :] = d2
    return newx, outdata

def apply_correction(data, maxshift=3):
    # data[1] = np.roll(data[0], 10)
    sumdata = data.sum(axis=0)
    fftsum = fft.rfft(sumdata)
    fftsum = np.conj(fftsum)

    fftdata = np.zeros((len(data), len(sumdata)))
    for i, d in enumerate(data):
        f1 = fft.rfft(d)
        f2 = fft.irfft(f1 * fftsum)
        fftdata[i, :len(f2)] = f2
    #     plt.plot(np.abs(fftdata[i]))
    # plt.show()
    fftdata[:, maxshift:len(sumdata)-maxshift] = 0
    maxes = np.argmax(fftdata, axis=1)
    maxes[maxes>len(sumdata)/2] = maxes[maxes>len(sumdata)/2] - len(sumdata)
    print(maxes)
    print("Number above max shift:", np.sum(np.abs(maxes) > maxshift))

    outdata = data.copy()
    for i, d in enumerate(data):
        if maxes[i] != 0 and np.abs(maxes[i])<maxshift:
            outdata[i] = np.roll(d, -maxes[i])

    return outdata




mpl.use("WxAgg")

file = r"C:\Data\Lipids\20260226_AS_004_HEK_UntDDA_BothPolarities_01.raw"
os.chdir(os.path.dirname(file))

importer = imp.create_importer(file)
print(importer.scans)
# s1 = importer.get_single_scan(1)
# s2 = importer.get_single_scan(2)
# x, ydata = apply_linearize([s1, s2], binsize=1)

data = importer.get_all_scans()
data = data[::1]
print("Got scans")
x, ydata = apply_linearize(data, binsize=0.1, normalize=True)
print("Linearized")
ydata2 = apply_correction(ydata, maxshift=3)
print("Plotting")
plt.plot(x, np.sum(ydata, axis=0))
plt.plot( x, np.sum(ydata2, axis=0))
# for y in ydata:
#     plt.plot(x, y)
plt.show()

# Seems to work decently but slow
# Idea to do a unit comparison of all vs all to find closest matches
# Then use higher res linearization and pairwise compare to find exact match

