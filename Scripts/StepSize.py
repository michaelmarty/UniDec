import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
#import skimage
import math


#def denoise(data, weight=0.1):
#    return skimage.restoration.denoise_tv_chambolle(data, weight)


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
    newdat = np.array(data).reshape((n, m))
    stds = np.std(newdat, axis=0)
    return np.median(stds)


def get_avg_ts(data):
    diffs = np.diff(data)
    return np.mean(diffs)


def make_segments(data, starts, ends):
    segments = []
    for i, s in enumerate(starts):
        e = ends[i]
        seg = data[s:e, 2]
        m = np.mean(seg)
        std = np.std(seg)
        ts = data[s, 1]
        te = data[e, 1]
        segments.append([m, std, ts, te])
    return np.array(segments)

def get_step_hist(diffs, binsize=0.1e-11):
    min = np.amin(diffs)
    max = np.amax(diffs)
    bins = np.arange(min, max, binsize)
    hist, edges = np.histogram(diffs, bins)
    



file = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Lab Data\\Zhihan\\ascii76.asc"

# Master Controllers
kernel_width = 5  # In Number of Data Points
noise_mult = 4
peak_window = 20  # In Number of Data Points
seg_pad = 3  # In Number of Data Points
max_val = 1e-11

data = np.loadtxt(file, delimiter=",")

#data[:,2] = denoise(data[:,2], weight=2)

nl = get_noise_level(data[:, 2])
print("Noise = ", nl)

ts = get_avg_ts(data[:, 1])
print("Data Length = ", len(data))
print("Avg Time Step = ", ts)

k = create_kernel(len(data), kernel_width)

conv = cconv2(data[:, 2], k)
cnl = np.mean(conv)
print("Convolution Noise = ", cnl)

cdata = np.transpose([data[:, 1], conv])

localmaxes = calc_local_maxes(cdata, window=peak_window)
localmaxes = filter_ints(localmaxes, cnl * noise_mult, max_val)

starts = localmaxes[:-1, 2].astype(int) + seg_pad
ends = localmaxes[1:, 2].astype(int) - seg_pad

segs = make_segments(data, starts, ends)
diffs = np.abs(np.diff(segs[:, 0]))

plt.subplot(121)
# plt.plot(data[:,1], data[:,2])
# data[:,2] = denoise(data[:,2])
plt.plot(data[:, 1], data[:, 2])
plt.plot(data[:, 1], conv)
plt.plot(localmaxes[:, 0], localmaxes[:, 1], marker="o", linestyle=" ")

for s in segs:
    plt.hlines(s[0], s[2], s[3], linewidth=3, color="k")

plt.subplot(122)
plt.hist(diffs, bins="auto")
plt.show()
