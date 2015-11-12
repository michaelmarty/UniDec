__author__ = 'Michael.Marty'

import numpy as np
from copy import deepcopy
from unidectools import *

file = "C:\\cprog\\UniDecDemo\\MTM_ADH_5_rawdata.txt"

#file = ".\\TestSpectra\\Traptavidin.txt"

import matplotlib.pyplot as plt

data = np.loadtxt(file, skiprows=10)
# boo1 = data[:, 0] > 4000
# boo2 = data[:, 0] < 7000
# boo3 = np.all([boo1, boo2], axis=0)
# data = data[boo3]
# data=gsmooth(data,10)
# data=datacompsub(data,100)
data[:, 1] = data[:, 1] / np.amax(data[:, 1])
data[:, 1] = np.clip(data[:, 1], 0, 1)

if False:
    plt.plot(data[:, 0], data[:, 1])
    plt.show()

fwhm, psfun, mid = auto_peak_width(data)
print fwhm
thresh = auto_noise_level(data) + 0.01
avgbs = average_bin_size(data)
window = fwhm / avgbs / 4.
print window, avgbs

cwt, wav = single_cwt(data[:, 1], window)

cwtdat = np.transpose([data[:, 0], cwt])

peaks = peakdetect(cwtdat, None, window * 2., thresh)
for p in peaks:
    p[1] = data[np.where(data[:, 0] == p[0]), 1]

zm = 1. / peaks[:, 0]

ivals = peaks[:, 1]
maxpos = np.argmax(ivals)
mins = []
for maxpos in xrange(0, len(ivals)):
    x = np.abs(zm - zm[maxpos])
    x = np.sort(x)
    diffs = np.diff(x)
    diffs = diffs / np.amax(diffs)
    minpos = np.argmin(diffs)
    minval = np.average(x[minpos:minpos + 2])
    mass = 1. / minval
    mz = peaks[maxpos, 0]
    charge = np.round(mass / mz)
    mass = mz * charge
    tolerance = charge * fwhm
    vec = np.array([mass, ivals[maxpos], mz, charge, diffs[minpos], tolerance, zm[maxpos], False])
    if charge > 100:
        vec = vec * 0
    mins.append(vec)
    if False:
        print vec
        plt.scatter(x, ivals)
        plt.plot(x[:-1], diffs)
        plt.xlim((np.amin(x), np.amax(x)))
        plt.show()

mins = np.array(mins)
sindex = np.argsort(mins[:, 4])
mins = mins[sindex]

peakassignments = []
for i in xrange(0, len(mins)):
    row = mins[i]
    swaps = []
    if row[0] > 0:
        topmass = row[0]
        topcharge = row[3]
        tolerance = row[5]
        print tolerance
        # print topmass,topcharge
        zms = np.abs(row[6] - mins[:, 6])
        for j, zm in enumerate(zms):
            if zm > 0 and mins[j, 2] > 0 and not mins[j, 7]:
                mass = 1. / zm

                charge = np.round(mass / mins[j, 2])
                newmass = mins[j, 2] * charge
                # print newmass
                if np.abs(newmass - topmass) < tolerance and np.abs(charge - topcharge) == 1:

                    mins[j, 0] = newmass
                    mins[j, 3] = charge
                    mins[j, 7] = True
                    mins[j, 6] = charge / newmass

                    print topmass, newmass

                    if not row[7]:
                        peakassignments.append([row, mins[j]])
                        row[7] = True
                    else:
                        for row2 in peakassignments:
                            for p in np.array(row2)[:, 0]:
                                if p == row[0]:
                                    row2.append(mins[j])

                    swaps.append(j)
    if False:
        for j in swaps:
            newspot = 0
            for k, row3 in enumerate(mins):
                if not row3[7]:
                    newspot = k
                    break

            if j > newspot:
                swaprows = deepcopy(mins[newspot:j])
                mins[newspot] = mins[j]
                mins[newspot + 1:j + 1] = swaprows
                print "swap", i, newspot, j

                # print i,j

    row[7] = True

boo1 = mins[:, 7] == 1

x = mins[boo1, 0]
y = mins[boo1, 4]
if True:
    plt.figure()
    plt.scatter(x, y)
    plt.xlim((np.amin(x), np.amax(x)))
    plt.ylim((np.amin(y), np.amax(y)))
    # plt.show()






# plt.plot(data[:,0],cwt)
plt.figure()
plt.plot(data[:, 0], data[:, 1])
plt.hlines(thresh, data[0, 0], np.amax(data[:, 0]), "r")

for p in peaks:
    plt.plot(p[0], p[1], marker="o", color="gray")

for m in peakassignments:
    m = np.array(m)
    x = m[:, 2]
    y = m[:, 1]
    sindex = np.argsort(x)
    x = x[sindex]
    y = y[sindex]
    print np.average(m[:, 0]), x, y
    plt.plot(x, y, marker="o")

# plt.plot(data[:, 0], fitdat)
plt.show()