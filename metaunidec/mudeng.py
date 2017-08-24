import os
import time
import shutil
from copy import deepcopy
import subprocess
import zipfile
import fnmatch
import string
import numpy as np
from scipy.interpolate import interp1d
from scipy import signal
import unidec_modules.unidectools as ud
from unidec_modules import unidecstructure, peakstructure, unidec_enginebase
from mudstruct import MetaDataSet

__author__ = 'Michael.Marty'


def metaunidec_call(config, *args, **kwargs):
    if "path" in kwargs:
        path=kwargs["path"]
        del kwargs["path"]
    else:
        path=config.hdf_file
    call = [config.UniDecPath, str(path)]
    if len(args) > 0:
        for arg in args:
            call.append(arg)
    if len(kwargs) > 0:
        for key in kwargs.keys():
            call.append("-" + str(key))
            call.append(kwargs[key])

    out = subprocess.call(call)
    print call, out
    return out


class MetaUniDec(unidec_enginebase.UniDecEngine):
    def __init__(self):
        """
        UniDec Engine

        Consists of three main subclasses: Config, DataContiner, Peaks

        :return: None
        """
        unidec_enginebase.UniDecEngine.__init__(self)
        self.data = MetaDataSet(self)
        self.config.filetype = 1
        self.config.metamode = -1
        self.config.linflag = 2

    def open(self, path):
        self.data = MetaDataSet(self)
        self.pks = peakstructure.Peaks()
        if path is None:
            path = self.config.hdf_file
        else:
            self.config.hdf_file = path
            self.config.filename = os.path.split(path)[1]
            self.config.outfname = os.path.splitext(self.config.filename)[0]
            self.config.outfname = os.path.join("UniDec_Figures_and_Files", self.config.outfname)
            dir = os.path.dirname(path)
            os.chdir(dir)
            dirnew = os.path.split(self.config.outfname)[0]
            if not os.path.isdir(dirnew):
                os.mkdir(dirnew)
            self.config.default_file_names()
        self.config.read_hdf5(path)
        self.data.import_hdf5(path)
        self.update_history()

    def process_data(self):
        self.pks.peaks = []
        self.config.write_hdf5()
        self.out=metaunidec_call(self.config, "-proc")
        self.data.import_hdf5()
        self.update_history()

    def run_unidec(self):
        if not self.check_badness():
            self.pks.peaks = []
            self.config.write_hdf5()
            self.out =metaunidec_call(self.config)
            self.data.import_hdf5()
            self.update_history()

    def make_grids(self):
        self.out =metaunidec_call(self.config,"-grids")

    def sum_masses(self):
        self.data.import_grids_and_peaks()

    def pick_peaks(self):
        self.config.write_hdf5()
        self.sum_masses()
        self.pks = peakstructure.Peaks()
        self.pks.add_peaks(self.data.peaks, massbins=self.config.massbins)
        self.pks.default_params(cmap=self.config.peakcmap)
        for i,p in enumerate(self.pks.peaks):
            p.extracts = self.data.exgrid[i]
        self.update_history()
        self.export_params()

    def export_params(self, e=None):
        peakparams = []
        for p in self.pks.peaks:
            peakparams.append([str(p.mass), str(p.height), str(p.area), str(p.label)])
        outfile = self.config.outfname + "_peaks.txt"
        np.savetxt(outfile, np.array(peakparams), delimiter=",", fmt="%s")

        peakexts=[]
        for p in self.pks.peaks:
            peakexts.append(np.concatenate(([p.mass],p.extracts)))
        outfile = self.config.outfname + "_extracts.txt"
        np.savetxt(outfile, np.array(peakexts))
        print "Peak info saved to:", outfile

    def export_spectra(self, e=None):
        for s in self.data.spectra:
            outfile = self.config.outfname + "_" +str(s.var1)+".txt"
            np.savetxt(outfile, s.rawdata)
            print outfile
            self.config.config_export(self.config.outfname + "_config.dat")

    def batch_set_config(self, paths):
        for p in paths:
            try:
                self.config.write_hdf5(p)
                print "Assigned Config to:",p
            except Exception, e:
                print e

    def batch_run_unidec(self, paths):
        for p in paths:
            try:
                tstart=time.clock()
                metaunidec_call(self.config,"-all", path=p)
                print "Run:", p, " Time:  %.3gs" % (time.clock() - tstart)
            except Exception, e:
                print e

    def batch_extract(self,paths):
        for p in paths:
            try:
                print "Extracting:", p
                self.open(p)
                self.pick_peaks()
            except Exception, e:
                print e

    def fit_data(self, fit="sig"):
        print "Fitting: ", fit
        xvals = self.data.var1
        self.data.fitgrid = []
        self.data.fits = []
        for i, y in enumerate(self.data.exgrid):
            if fit == "exp":
                fits, fitdat = ud.exp_fit(xvals, y)
            elif fit == "lin":
                fits, fitdat = ud.lin_fit(xvals, y)
            elif fit == "sig":
                fits, fitdat = ud.sig_fit(xvals, y)
            else:
                print "ERROR: Unsupported fit type"
                break
            print fits
            self.data.fitgrid.append(fitdat)
            self.data.fits.append(fits)
        self.data.export_fits()


if __name__ == '__main__':

    eng = MetaUniDec()
    '''
    testpath = "C:\Python\UniDec\unidec_src\UniDec\\x64\Release\\test.hdf5"
    eng.data.new_file(testpath)
    data1 = [1, 2, 3]
    data2 = [4, 5, 6]
    data3 = [7, 8, 9]
    eng.data.add_data(data1)
    eng.data.add_data(data2)
    eng.data.add_data(data3)
    eng.data.remove_data([0, 2])
    exit()
    '''

    testdir = "C:\Python\UniDec\unidec_src\UniDec\\x64\Release"
    testfile = "JAW.hdf5"
    testpath = os.path.join(testdir, testfile)
    eng.open(testpath)
    eng.run_unidec()
    eng.pick_peaks()
