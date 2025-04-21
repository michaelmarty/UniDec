import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
#import skimage
import math
import os
from scipy import signal

#def denoise(data, weight=0.1):
#    return skimage.restoration.denoise_tv_chambolle(data, weight)

def autocorr(datatop):
    """
    Note: ASSUMES LINEARIZED DATA
    :param datatop: 1D data
    :param config: Config file (optional file)
    :return: Autocorr spectrum, peaks in autocorrelation.
    """
    corry = signal.fftconvolve(datatop[:, 1], datatop[:, 1][::-1], mode='same')
    corry /= np.amax(corry)
    maxpos1 = np.argmax(datatop[:, 1])
    start = np.amax([maxpos1 - len(datatop) / 10, 0])
    end = np.amin([len(datatop) - 1, maxpos1 + len(datatop) / 10])
    cutdat = datatop[int(start):int(end)]
    if len(cutdat) < 20:
        cutdat = datatop
    xdiff = np.mean(cutdat[1:, 0] - cutdat[:len(cutdat) - 1, 0])
    corrx = np.arange(0.0, len(corry)) * xdiff
    maxpos = np.argmax(corry)
    corrx = corrx - corrx[maxpos]
    autocorr = np.transpose([corrx, corry])
    boo1 = autocorr[:, 0] > xdiff
    cpeaks = calc_local_maxes(autocorr[boo1])
    return autocorr, cpeaks

def create_kernel(dlen, steplen):
    k = np.zeros(dlen)
    k[:steplen] = 0.5
    k[:-steplen - 1:-1] = -0.5
    return k


def cconv2(a, b):
    A = fftpack.fft(a)
    B = fftpack.fft(b)
    C = A * B
    c = np.abs(fftpack.ifft(C))
    return c


def calc_local_maxes(data, window=5, threshold=0):
    peaks = []
    length = len(data)
    for i in range(0, length):
        if data[i, 1] > threshold:
            start = i - window
            end = i + window
            if start < 0:
                start = 0
            if end > length:
                end = length
            start = int(start)
            end = int(end) + 1
            testmax = np.amax(data[start:end, 1])
            if data[i, 1] == testmax and np.all(data[i, 1] != data[start:i, 1]):
                peaks.append([data[i, 0], data[i, 1], i])
    return np.array(peaks)


def filter_ints(data, min, max):
    b1 = data[:, 1] > min
    b2 = data[:, 1] < max
    b3 = np.logical_and(b1, b2)
    return data[b3]


def get_noise_level(data, n=10):
    m = math.floor(len(data) / 10)
    newdat = np.array(data[:n*m]).reshape((n, m))
    stds = np.std(newdat, axis=0)
    return np.median(stds)


def get_avg_ts(data):
    diffs = np.diff(data)
    return np.mean(diffs)


def make_segments(data, starts, ends):
    segments = []
    for i, s in enumerate(starts):
        e = ends[i]
        seg = data[s:e, 1]
        m = np.mean(seg)
        std = np.std(seg)
        ts = data[s, 0]
        te = data[e, 0]
        segments.append([m, std, ts, te])
    return np.array(segments)


file = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Lab Data\\Zhihan\\ascii76.asc"
file = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Lab Data\\Zhihan\\VpuTrace.txt"

# Master Controllers
kernel_width = 5  # In Number of Data Points
noise_mult = 6
peak_window = 20  # In Number of Data Points
seg_pad = 3  # In Number of Data Points
max_val = 100 #1e-11
diff_cutoff = 0.0 # In multiple of the noise level

if os.path.splitext(file)[1]==".asc":
    data = np.loadtxt(file, delimiter=",", usecols=(1,2))
else:
    data = np.loadtxt(file, skiprows=1)
#data[:,2] = denoise(data[:,2], weight=2)

nl = get_noise_level(data[:, 1])
print("Noise = ", nl)

ts = get_avg_ts(data[:, 0])
print("Data Length = ", len(data))
print("Avg Time Step = ", ts)

# Convolve the data with a step function to find where steps are
k = create_kernel(len(data), kernel_width)
conv = cconv2(data[:, 1], k)
cnl = np.mean(conv)
print("Convolution Noise = ", cnl)
cdata = np.transpose([data[:, 0], conv])
localmaxes = calc_local_maxes(cdata, window=peak_window)
localmaxes = filter_ints(localmaxes, cnl * noise_mult, max_val)

# Make Segments and get difference between those
starts = localmaxes[:-1, 2].astype(int) + seg_pad
ends = localmaxes[1:, 2].astype(int) - seg_pad
segs = make_segments(data, starts, ends)
diffs = np.abs(np.diff(segs[:, 0]))
try:
    segs = make_segments(data, starts, ends)
    diffs = np.abs(np.diff(segs[:, 0]))
    diffs = diffs[diffs>nl*diff_cutoff]
except:
    segs = []
    diffs =[]

# Histogram of the raw data
dhist, edges = np.histogram(data[:,1], bins="auto")
edges = edges[1:]-(edges[1]-edges[0])/2.
histdat = np.transpose([edges, dhist])
autocorrdat, acpeaks = autocorr(histdat)
autocorrdat = autocorrdat[autocorrdat[:,0]>0]

plt.subplot(221)
# plt.plot(data[:,1], data[:,2])
# data[:,2] = denoise(data[:,2])
plt.plot(data[:, 0], data[:, 1])
plt.plot(data[:, 0], conv)
plt.plot(localmaxes[:, 0], localmaxes[:, 1], marker="o", linestyle=" ")
plt.xlabel("Time")
plt.ylabel("Current")
for s in segs:
    plt.hlines(s[0], s[2], s[3], linewidth=3, color="k")

plt.subplot(222)
if len(diffs)>0:
    plt.hist(diffs, bins="auto")
plt.xlabel("Step Size (pA)")
plt.ylabel("Count")
plt.subplot(223)

plt.hist(data[:,1], bins="auto")
plt.plot(histdat[:,0], histdat[:,1])
plt.xlabel("Current (pA)")
plt.ylabel("Count")

plt.subplot(224)
plt.plot(autocorrdat[:,0], autocorrdat[:,1])
plt.xlabel("Current Autocorrelation")

plt.show()
