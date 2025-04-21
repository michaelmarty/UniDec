from unidec import engine
import numpy as np
import unidec.tools as ud
import multiprocessing


class UniDecTest(engine.UniDec):
    def test_spectra(self, mlist=None, ilist=None, res=1000, noise=0.0, pad=2000, window=None,
                     plot=False, restricted=True, massbins=10, psfun=0, mzsig=None, **kwargs):
        if mlist is None:
            mlist = [200000, 210000]
        if ilist is None:
            ilist = [100, 100]

        m1 = mlist[0]
        m2 = mlist[1]
        if mzsig is None:
            mzsig = m1 / float(res) / float(ud.predict_charge(m1)) / 2.

        self.open_test_spectrum(masslist=mlist, intlist=ilist, noise=noise, resolution=res, psfun=psfun, silent=True,
                                **kwargs)
        self.config.massbins = massbins
        self.config.mzsig = mzsig
        self.config.psfun = psfun
        self.config.linflag = 2
        self.config.mzbins = 0
        self.config.peaknorm = 2
        self.config.poolflag = 2
        #self.config.massbins = 1000
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
        self.get_spectrum_peaks()
        self.run_unidec(silent=True)
        self.config.peakthresh = 0.1
        self.pick_peaks()
        if plot:
            self.make_plot(massrange=[m1 - pad, m2 + pad])


if __name__ == "__main__":
    multiprocessing.freeze_support()
    eng = UniDecTest()
    eng.test_spectra(mlist=[200000, 202000], res=1000, plot=True, restricted=True, pad=2000)
