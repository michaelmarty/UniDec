import numpy as np
import os
import matplotlib.pyplot as plt
import unidec_modules.unidectools as ud
from copy import deepcopy
from scipy import fftpack
import scipy.stats as stats

'''
def weiner(data, kernel):
    G = fftpack.fft(kernel)
    H = fftpack.fft(data)
    G2 = (1 / (1 + 1 / (H * (np.abs(G) ** 2)))) / G
    F = H * G2
    f = fftpack.ifft(F)
    f = np.abs(f)
    f = f / np.amax(f)
    return f'''


def cconv2(a, b):
    A = fftpack.fft(a)
    B = fftpack.fft(b)
    C = A * B
    c = np.abs(fftpack.ifft(C))
    return c


def make_kernel(data):
    try:
        loc = np.argmax(data[:, 1])
        part1 = data[loc:, 1]
        part2 = data[:loc, 1]
    except:
        loc = np.argmax(data)
        part1 = data[loc:]
        part2 = data[:loc]
    output = np.append(part1, part2, axis=0)
    output /= np.amax(output)
    return output


def import_kernel_and_run(data, kfile):
    dd = DoubleDec()
    dd.kimport(kfile)
    dd.dd_run(data=data)
    return dd.dec2


def batch_dd(files, kfile):
    for f in files:
        outfile = os.path.splitext(f)[0]+"dd.txt"
        dd = DoubleDec()
        dd.dd_import(f, kfile)
        dd.dd_run()
        np.savetxt(outfile, dd.dec2)
        print(f, kfile, outfile)
    pass


class DoubleDec:
    def __init__(self, *args, **kwargs):
        self.data = None
        self.kernel = None
        self.dec = None
        self.dec2 = None
        self.extracts = None
        pass

    def dd_prep(self):
        try:
            self.kernel[:, 1] /= np.amax(self.kernel[:, 1])
        except:
            self.kernel /= np.amax(self.kernel)
        self.data[:, 1] /= np.amax(self.data[:, 1])

        if len(self.kernel) > len(self.data):
            self.data = ud.pad_data_length(self.data, pad_until_length=len(self.kernel))
        elif len(self.data) > len(self.kernel):
            self.kernel = ud.pad_data_length(self.kernel, pad_until_length=len(self.data))
        self.kernel = make_kernel(self.kernel)

    def dd_core(self, data, kernel):
        ckernel = kernel[::-1]
        I = deepcopy(data)
        i = 0
        diff = 1
        while i < 50 and diff > 0.0001:
            newI = I * cconv2(ud.safedivide(data, cconv2(kernel, I)), ckernel)
            diff = np.sum((I - newI) ** 2) / np.sum(I)
            I = newI
            i += 1
        I /= np.amax(I)
        return I

    def dd_run(self, data=None, kernel=None):
        if data is not None:
            self.data = data
        if kernel is not None:
            self.kernel = kernel
        self.dd_prep()
        self.dec = self.dd_core(self.data[:, 1], self.kernel)
        self.dec2 = np.transpose([self.data[:, 0], self.dec])

    def dd_import(self, datapath, kernelpath):
        self.kimport(kernelpath)
        self.data = np.loadtxt(datapath)

    def kimport(self, kernelpath):
        self.kernel = np.loadtxt(kernelpath)
        self.kdata = deepcopy(self.kernel)

    def Extract(self, data, basemass=41983, m1=762, m2=63, m1range=[0, 7], m2range=[0, 7], exmethod=1, window=10):
        self.nm1 = np.arange(m1range[0], m1range[1])
        self.nm2 = np.arange(m2range[0], m2range[1])
        self.m2grid, self.m1grid = np.meshgrid(self.nm2, self.nm1)
        self.masses = basemass + m1 * self.m1grid + m2 * self.m2grid
        self.extracts = np.array(
            [[ud.data_extract(data, self.masses[i, j], extract_method=exmethod, window=window) for i in self.nm1] for j
             in self.nm2])
        self.extracts /= np.amax(self.extracts)
        return self.extracts

    def WeightedAvgs(self, cutoff=0):
        self.im2 = np.sum(self.extracts, axis=1)
        self.im1 = np.sum(self.extracts, axis=0)

        self.m2array = np.array([ud.weighted_avg(self.nm2, e) for e in np.transpose(self.extracts)])
        self.m1array = np.array([ud.weighted_avg(self.nm1, e) for e in self.extracts])

        b1 = self.m2array != 0
        b2 = self.m1array != 0
        # print(self.m2array[b1], self.nm1[b1], self.m1array[b2], self.nm2[b2])
        print(stats.pearsonr(self.m2array[b1], self.nm1[b1]), stats.pearsonr(self.m1array[b2], self.nm2[b2]))

        b1 = self.im1>cutoff*np.amax(self.im1)
        b2 = self.im2>cutoff * np.amax(self.im2)
        self.wm2 = ud.weighted_avg(self.nm2[b2], self.im2[b2])
        self.wm1 = ud.weighted_avg(self.nm1[b1], self.im1[b1])
        print(self.wm1, self.wm2)
        return self.wm1, self.wm2

    def PlotPeaks(self):
        peaks = ud.peakdetect(self.dec2, window=10, threshold=0.1)
        spacing = 0.4
        simdat = deepcopy(self.data)
        simdat[:, 1] = 0
        for i, p in enumerate(peaks):
            sticks = self.dec2[:, 0] == p[0]

            sticks = sticks.astype(np.float) * p[1]
            sticks = cconv2(sticks, self.kernel)

            simdat[:, 1] += sticks
            # plt.plot(simdat[:,0], sticks-(i+1)*spacing)

        sim2 = cconv2(self.dec, self.kernel)
        sim2 /= np.amax(sim2)
        # plt.plot(simdat[:,0], simdat[:,1], color="r")
        plt.plot(self.data[:, 0], self.data[:, 1], color="k")
        plt.plot(self.data[:, 0], self.dec + spacing * 2, color="b")
        plt.plot(self.data[:, 0], sim2, color="r")
        # plt.plot(dec[:,0], kernel)

        # plt.xlim(41500, 43000)
        plt.show()

    def plot2(self):
        plt.subplot(121)
        plt.plot(self.dec2[:, 0], self.dec2[:, 1])
        plt.plot(self.data[:, 0], self.data[:, 1])
        plt.subplot(122)
        plt.imshow(self.extracts)
        plt.show()

    def waplot(self):
        plt.plot(self.nm1, self.m2array)
        plt.plot(self.nm2, self.m1array)
        plt.show()


if __name__ == "__main__":
    dir = "D:\\Data\\rhodopsin"

    f1 = "0.txt"
    f2 = "250.txt"

    p1 = os.path.join(dir, f1)
    p2 = os.path.join(dir, f2)

    # p1 = "D:\Data\WangLabAnalysis\\20200430_jat_wanglabsample_concthenbiospin_unidecfiles\\20200430_jat_wanglabsample_concthenbiospin_mass.txt"
    # p1 = "D:\Data\WangLabAnalysis\wanglabsample_protease_3_unidecfiles\wanglabsample_protease_3_mass.txt"
    # p2 = "D:\Data\WangLabAnalysis\Titration\\20200502_jat_wanglabsample_240_500uMadded_1_unidecfiles\\20200502_jat_wanglabsample_240_500uMadded_1_mass.txt"
    # p2 ="D:\Data\WangLabAnalysis\\200515\\244\\20200514_JAT_WangLabSample_244_175uMadded_200515115154_unidecfiles\\20200514_JAT_WangLabSample_244_175uMadded_200515115154_mass.txt"

    dd = DoubleDec()
    dd.dd_import(p2, p1)

    # plt.plot(dd.kernel[:,0], dd.kernel[:,1])
    # plt.show()
    dd.dd_run()
    extracts = dd.Extract(dd.dec2)
    dd.WeightedAvgs()
    dd.PlotPeaks()
    # dd.plot2()
