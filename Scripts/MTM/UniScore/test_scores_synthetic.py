import numpy as np
import multiprocessing
import unidec.modules.isolated_packages.score_window as sw
from unidec import engine, tools as ud
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib import rcParams
import os

import matplotlib as mpl

mpl.use("TkAgg")

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
__author__ = 'Michael.Marty'

import warnings

warnings.filterwarnings("ignore")

luminance_cutoff = 135

dir = "C:\\Users\\Michael Marty\\University of Arizona\\martylab - Documents\\Papers\\UniScore\\Figures\\Synthetic"
os.chdir(dir)

save = True
show = True
scorewindow = False


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
                     plot=False, restricted=True, massbins=10, psfun=0, mzsig=None, massrange=None, zrange=None,
                     **kwargs):
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

        if massrange is not None:
            self.config.masslb = massrange[0]
            self.config.massub = massrange[1]
            print(massrange)

        if zrange is not None:
            self.config.endz = zrange[1]
            self.config.startz = zrange[0]
            print(zrange)

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

        if scorewindow:
            sw.score_window(self)
        pass

    def plot_test(self, testdat, testmdat, peaks, fitdat=None):
        plt.rcParams.update({'font.size': 10})
        colormap = cm.get_cmap('RdYlGn')
        norm = Normalize(vmin=0, vmax=1)

        colormap2 = cm.get_cmap('viridis')
        norm2 = Normalize(vmin=0, vmax=len(testdat) - 1)

        plt.figure(figsize=(8, 6))
        offset = 100

        xmin = 0
        xmin2 = 1000000
        xmax = 10000000
        for i, d in enumerate(testdat):
            m = testmdat[i]
            minval = np.amin(m[:, 0])
            maxval = np.amax(m[:, 0])
            if minval > xmin:
                xmin = minval
            if maxval < xmax:
                xmax = maxval
            if minval < xmin2:
                xmin2 = minval

        upos = xmin2 + (xmax - xmin2) * 0.03
        abcindex = ord("a")
        for i, d in enumerate(testdat):
            color2 = colormap2(norm2(i))
            lum = ud.get_luminance(color2)
            if lum > luminance_cutoff / 255.:
                color2 = np.array(color2) * 0.8
                color2[3] = 1

            p1 = plt.subplot(121)
            plt.plot(d[:, 0], d[:, 1] * 100 - i * offset, color=color2)
            if fitdat is not None:
                plt.plot(d[:, 0], fitdat[i] * 100 - i * offset, color="r")
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

        # plt.xlim(xmin / 1000, xmax / 1000)
        plt.ylim(100 - (len(testdat) * 100), 100)
        plt.tight_layout()
        # plt.show()

    def test_noise(self, **kwargs):
        noiselist = [0, 0.001, 0.01, 0.1, 0.3]
        testdat = []
        testmdat = []
        peaks = []
        for n in noiselist:
            self.test_spectra(res=1000, noise=n, pad=2000, window=2000, plot=False, **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(self.pks)
        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("noisetest.png")
            plt.savefig("noisetest.pdf")
        if show:
            plt.show()

    def test_width(self, **kwargs):
        reslist = [300, 500, 1000, 2000, 5000]
        testdat = []
        testmdat = []
        peaks = []
        for res in reslist:
            self.test_spectra(res=res, noise=0.01, pad=2000, window=2000, plot=False, **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(self.pks)

        os.chdir(dir)
        np.savetxt("5000test.txt", self.data.data2)

        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("widthtest.png")
            plt.savefig("widthtest.pdf")
        if show:
            plt.show()

    def test_zwidth(self, **kwargs):
        zwlist = [0.5, 1, 3]
        testdat = []
        testmdat = []
        peaks = []
        for zw in zwlist:
            self.test_spectra(res=1000, noise=0.01, pad=20000, window=2000, plot=False, zwidth=zw, **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(self.pks)
        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("zwidthtest.png")
            plt.savefig("zwidthtest.pdf")
        if show:
            plt.show()

    def test_baseline(self, **kwargs):
        blist = [0, 0.05, 0.1, 0.25]
        testdat = []
        testmdat = []
        peaks = []
        for b in blist:
            self.test_spectra(res=1000, noise=0.01, pad=20000, window=2000, plot=False, baseline=b, **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(self.pks)
        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("baselinetest.png")
            plt.savefig("baselinetest.pdf")
        if show:
            plt.show()

    def test_ambiguity(self, **kwargs):
        m1 = 200000
        m2list = [202000, 210000, 215000, 205000, 220000]

        testdat = []
        testmdat = []
        peaks = []
        for m2 in m2list:
            mlist = [m1, m2]
            self.test_spectra(mlist=mlist, res=300, noise=0.0, pad=25000, window=1000, plot=False, **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            peaks.append(self.pks)

        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("ambiguitytest.png")
            plt.savefig("ambiguitytest.pdf")
        if show:
            plt.show()

    def test_density(self, **kwargs):
        massbinlist = [1, 10, 100]

        testdat = []
        testmdat = []
        peaks = []
        num = 1
        for b in massbinlist:
            for n in range(num):
                self.test_spectra(res=1000, noise=0.01, pad=20000, window=2000, plot=False, massbins=b, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(self.pks)
        self.plot_test(testdat, testmdat, peaks)

        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("densitytest.png")
            plt.savefig("densitytest.pdf")
        if show:
            plt.show()

    def test_psfun(self, **kwargs):
        psfuns = [0, 2, 1]
        testdat = []
        testmdat = []
        peaks = []
        num = 1
        for p in psfuns:
            for n in range(num):
                self.test_spectra(res=1000, noise=0.0, pad=20000, window=2000, plot=False, psfun=p, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(self.pks)

        self.plot_test(testdat, testmdat, peaks)

        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("psfun.png")
            plt.savefig("psfun.pdf")
        if show:
            plt.show()

    def test_widths(self, **kwargs):
        r1 = 1000
        r2list = [300, 500, 800, 1200, 3000]
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 202000]
        num = 1
        for r2 in r2list:
            rlist = [r1, r2]
            for n in range(num):
                self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.01, pad=3000,
                                  window=1000, plot=False,
                                  psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(self.pks)

        self.plot_test(testdat, testmdat, peaks)

        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("widthtest2.png")
            plt.savefig("widthtest2.pdf")
        if show:
            plt.show()

    def test_sigs(self, **kwargs):
        rlist = [500, 500]
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 202000]
        sigs = [0, 10, 25]
        num = 1
        for sig in sigs:
            for n in range(0, num):
                self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.02, pad=20000,
                                  window=1000,
                                  plot=False, mzsig=sig,
                                  psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(self.pks)
        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("deconsigtest.png")
            plt.savefig("deconsigtest.pdf")
        if show:
            plt.show()

    def test_sigs_narrow(self, **kwargs):
        rlist = [500, 500]
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 202000]
        sigs = [0, 10, 25]
        num = 1
        for sig in sigs:
            for n in range(0, num):
                self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.02, pad=2000,
                                  window=1000,
                                  plot=False, mzsig=sig,
                                  psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(self.pks)
        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("deconsigtestnarrow.png")
            plt.savefig("deconsigtestnarrow.pdf")
        if show:
            plt.show()

    def test_repeat(self, **kwargs):
        rlist = [1000, 2000]
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 210000]
        sigs = [1]
        num = 5
        for sig in sigs:
            for n in range(0, num):
                self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.1, pad=2000,
                                  window=1000,
                                  plot=False, mzsig=sig,
                                  psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(self.pks)
        # self.plot_repeats(peaks)
        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("repeat.png")
            plt.savefig("repeat.pdf")
        if show:
            plt.show()

    def test_repeat_low(self, **kwargs):
        rlist = [1000, 2000]
        testdat = []
        testmdat = []
        peaks = []
        mlist = [200000, 210000]
        sigs = [1]
        num = 5
        for sig in sigs:
            for n in range(0, num):
                self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.01, pad=2000,
                                  window=1000,
                                  plot=False, mzsig=sig,
                                  psfun=0, plot_corr=False, plot_fits=False, **kwargs)
                testdat.append(self.data.data2)
                testmdat.append(self.data.massdat)
                peaks.append(self.pks)
        # self.plot_repeats(peaks)
        self.plot_test(testdat, testmdat, peaks)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("repeatlow.png")
            plt.savefig("repeatlow.pdf")
        if show:
            plt.show()

    def test_missing(self, **kwargs):
        rlist = [1000, 2000]
        testdat = []
        testmdat = []
        fitdat = []
        peaks = []
        mlist = [200000, 210000]
        massranges = [[198000, 212000], [198000, 212000], [198000, 212000], [208000, 212000]]
        zrange = [[20, 40], [30, 40], [32, 40], [20, 40]]

        for i, x in enumerate(massranges):
            self.test_spectra(mlist=mlist, rlist=rlist, res=1000, noise=0.01, pad=20000,
                              window=1000,
                              plot=False, mzsig=1, massrange=x, zrange=zrange[i],
                              psfun=0, plot_corr=False, plot_fits=False, **kwargs)
            testdat.append(self.data.data2)
            testmdat.append(self.data.massdat)
            fitdat.append(self.data.fitdat)
            peaks.append(self.pks)
        # self.plot_repeats(peaks)
        self.plot_test(testdat, testmdat, peaks, fitdat=fitdat)
        if save:
            os.chdir(dir)
            print("Saving", os.getcwd())
            plt.savefig("missing.png")
            plt.savefig("missing.pdf")
        if show:
            plt.show()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    eng = UniDecTest()

    # eng.test_noise()
    # eng.test_baseline()
    # eng.test_ambiguity()
    # eng.test_zwidth()
    # eng.test_width(mz_bin_size=0.1)
    #eng.test_repeat_low()
    # eng.test_repeat()
    # eng.test_sigs()
    # eng.test_sigs_narrow()
    eng.test_missing()

    # eng.test_density()
    # eng.test_psfun()
    # eng.test_widths()


