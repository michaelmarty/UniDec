import numpy as np
import multiprocessing
from unidec import engine, tools as ud
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib import rcParams
import os

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
__author__ = 'Michael.Marty'

import warnings

warnings.filterwarnings("ignore")

luminance_cutoff = 135

dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\UniScore\\Figures\\Synthetic"
os.chdir(dir)

save = False
show = True


# Test for error function:
# 1. Peak Width
# 2. Noise
# 3. Number of Charge States
# 4. Baseline
# 5. Ambiguity
# Invariant towards:
# 1. Data density
# 2. Peak Shape
# 3. Different peak widths
# 4. Different Intensities
# 5. MZsig

class UniDecTest(engine.UniDec):
    def test_spectra(self, mlist=None, ilist=None, res=1000, noise=0.0, pad=1000, window=None,
                     plot=False, restricted=True, massbins=10, psfun=0, mzsig=None, **kwargs):
        if mlist is None:
            mlist = [200000, 210000]
        if ilist is None:
            ilist = [50, 100]

        m1 = mlist[0]
        m2 = mlist[1]
        if mzsig is None:
            mzsig = m1 / float(res) / float(ud.predict_charge(m1)) / 10.

        self.open_test_spectrum(masslist=mlist, intlist=ilist, noise=noise, resolution=res, psfun=psfun, silent=True,
                                **kwargs)
        self.config.massbins = massbins
        self.config.mzsig = mzsig
        self.config.psfun = psfun
        self.config.linflag = 2
        self.config.mzbins = 0
        self.config.peaknorm = 2
        if window is None:
            self.config.peakwindow = pad
        else:
            self.config.peakwindow = window
        if restricted:
            self.config.masslb = m1 - pad
            self.config.massub = m2 + pad
            self.config.endz = np.amax(self.data.ztab)
            self.config.startz = np.amin(self.data.ztab)
        self.process_data(silent=True)
        self.run_unidec(silent=True)
        self.config.peakthresh = 0.1
        self.pick_peaks()
        self.autointegrate()
        # self.get_errors(**kwargs)
        self.dscore(xfwhm=3)
        ''''''
        if plot:
            self.make_plot(massrange=[m1 - pad, m2 + pad])
            # self.config.print_config()
        pass

    def plot_test(self, testdat, testmdat, peaks):
        plt.rcParams.update({'font.size': 10})
        colormap = cm.get_cmap('RdYlGn')
        norm = Normalize(vmin=0, vmax=1)

        colormap2 = cm.get_cmap('viridis')
        norm2 = Normalize(vmin=0, vmax=len(testdat) - 1)

        plt.figure(figsize=(8, 6))
        offset = 100

        xmin = 0
        xmax = 10000000
        for i, d in enumerate(testdat):
            m = testmdat[i]
            minval = np.amin(m[:, 0])
            maxval = np.amax(m[:, 0])
            if minval > xmin:
                xmin = minval
            if maxval < xmax:
                xmax = maxval

        upos = xmin + (xmax - xmin) * 0.03
        abcindex = ord("a")
        for i, d in enumerate(testdat):
            color2 = colormap2(norm2(i))
            lum = ud.get_luminance(color2)
            if lum > luminance_cutoff / 255.:
                color2 = np.array(color2) * 0.8
                color2[3] = 1

            p1 = plt.subplot(121)
            plt.plot(d[:, 0], d[:, 1] * 100 - i * offset, color=color2)
            plt.xlabel("$\it{m/z}$")
            p1.get_yaxis().set_ticks([0, 50, 100])
            p1.get_yaxis().set_ticklabels(["0", '%', "100"])
            p1.spines['top'].set_visible(False)
            p1.spines['right'].set_visible(False)
            plt.ylim(100 - (len(testdat) * 100), 100)
            plt.xlim(np.amin(d[:, 0]), np.amax(d[:, 0]))

            p2 = plt.subplot(122)
            plt.xlabel("Mass (kDa)")
            p2.get_yaxis().set_ticks([0, 50, 100])
            p2.get_yaxis().set_ticklabels(["0", '%', "100"])
            p2.spines['top'].set_visible(False)
            p2.spines['right'].set_visible(False)

            m = testmdat[i]
            plt.plot(m[:, 0] / 1000, m[:, 1] - i * offset, color=color2)

            for p in peaks[i].peaks:
                # plt.plot(p.mass, p.integral - i * offset, color="g", marker='o')
                text = str(int(round(p.dscore * 100)))
                color = colormap(norm(p.dscore))
                lum = ud.get_luminance(color)
                if lum > luminance_cutoff / 255.:
                    # print(p.mass, lum, color)
                    color = np.array(color) * 0.8
                    color[3] = 1
                plt.text(p.mass / 1000, 20 + p.height - i * offset, text, horizontalalignment="center",
                         verticalalignment="top", color=color)
            text = "UniScore=" + str(int(round(peaks[i].uniscore * 100)))
            color = colormap(norm(peaks[i].uniscore))
            lum = ud.get_luminance(color)
            if lum > luminance_cutoff / 255.:
                # print(p.mass, lum, color)
                color = np.array(color) * 0.8
                color[3] = 1
            p2.text(upos / 1000, 90 - i * offset, text, horizontalalignment="left",
                    verticalalignment="top", color=color)

            label = chr(abcindex).upper()
            abcindex += 1
            p1.text(np.amax(d[:, 0]), 92 - i * offset, label, horizontalalignment="right",
                    verticalalignment="top", color="k", fontsize=12)

            label2 = chr(abcindex).upper()
            abcindex += 1

            p2.text(xmax / 1000, 92 - i * offset, label2, horizontalalignment="right",
                    verticalalignment="top", color="k", fontsize=12)

        plt.xlim(xmin / 1000, xmax / 1000)
        plt.ylim(100 - (len(testdat) * 100), 100)
        plt.tight_layout()
        # plt.show()

    def test_FDR(self, **kwargs):
        rlist = [1000, 2000]
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 210000]
        sigs = [0]
        num = 50
        dscores = []
        for sig in sigs:
            for n in range(0, num):
                self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.001, pad=2000,
                                  window=1000,
                                  plot=False, mzsig=sig, scramble=True,
                                  psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(self.pks)
                for p in self.pks.peaks:
                    dscores.append(p.dscore)
        # self.plot_repeats(peaks)
        #self.plot_test(testdat, testmdat, peaks)

        #os.chdir(dir)
        #print("Saving", os.getcwd())
        #np.savetxt("noise50its.txt", dscores)

        analyzeFDR(dscores)

        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("fdr.png")
            plt.savefig("fdr.pdf")
        if show:
            plt.show()


def analyzeFDR(d=None):
    if d is None:
        os.chdir(dir)
        d = np.loadtxt("noise100its.txt")

    sd = np.sort(d)
    l = len(sd)

    plt.subplot(121)
    plt.plot(sd)

    fdrs = [0,0.001, 0.01, 0.05]
    for fdr in fdrs:
        f = int(l*(1-fdr))+1
        if f >= l:
            f = l-1
        print("FDR:", fdr, "Value:", sd[f])
        plt.plot(f, sd[f], marker="o", color="red")


    plt.subplot(122)
    bins = np.arange(0,1,0.05)
    plt.hist(d, bins=bins)



if __name__ == "__main__":

    if True:
        multiprocessing.freeze_support()
        eng = UniDecTest()

        eng.test_FDR()
