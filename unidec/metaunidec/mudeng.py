import os
import numpy as np

import unidec.tools as ud
from unidec.modules import peakstructure, unidec_enginebase
from unidec.metaunidec.mudstruct import MetaDataSet
import unidec.modules.mzmlparse_auto as automzml
import time

try:
    from pyimzml.ImzMLWriter import ImzMLWriter
    from unidec.metaunidec.imzml_reader import imzml_to_hdf5
except:
    print("pyimzML not found. Imaging features won't work.")

__author__ = 'Michael.Marty'


def metaunidec_call(config, *args, **kwargs):
    if "path" in kwargs:
        path = kwargs["path"]
        del kwargs["path"]
    else:
        path = config.hdf_file
    call = [config.UniDecPath, str(path)]
    if len(args) > 0:
        for arg in args:
            call.append(arg)
    if len(kwargs) > 0:
        for key in list(kwargs.keys()):
            call.append("-" + str(key))
            call.append(kwargs[key])
    tstart = time.perf_counter()
    out = ud.exe_call(call)
    tend = time.perf_counter()
    # print(call, out)
    print("Execution Time:", (tend - tstart))
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

    def open(self, path, speedy=False):
        st = time.perf_counter()
        self.clear()
        if path is None:
            path = self.config.hdf_file
        else:
            self.setup_filenames(path)
        self.config.read_hdf5(path)
        self.data.import_hdf5(path, speedy=speedy)
        self.update_history()
        print("Opening Time:", time.perf_counter()-st)

    def clear(self):
        self.data = MetaDataSet(self)
        self.pks = peakstructure.Peaks()

    def setup_filenames(self, path):
        self.config.hdf_file = path
        self.config.dirname, self.config.filename = os.path.split(path)
        outfname = os.path.splitext(self.config.filename)[0]
        dir = os.path.dirname(path)
        dirnew = os.path.join(dir, "UniDec_Figures_and_Files")
        if not os.path.isdir(dirnew):
            os.mkdir(dirnew)
        self.config.udir = dirnew
        self.config.outfname = os.path.join(dirnew, outfname)
        self.config.default_file_names()

    def process_data(self):
        self.pks.peaks = []
        self.config.write_hdf5()
        self.out = metaunidec_call(self.config, "-proc")
        self.data.import_hdf5()
        self.update_history()

    def run_unidec(self):
        if not self.check_badness():
            self.pks.peaks = []
            self.config.write_hdf5()
            self.out = metaunidec_call(self.config)
            self.data.import_hdf5()
            self.update_history()

    def make_grids(self):
        self.out = metaunidec_call(self.config, "-grids")

    def sum_masses(self, refresh=True):
        self.data.import_grids_and_peaks(refresh=refresh)

    def pick_peaks(self, refresh=True):
        if refresh:
            self.config.write_hdf5()
        s = time.perf_counter()
        self.sum_masses(refresh=refresh)

        scores_included = False
        try:
            if len(self.data.peaks[0]) == 3:
                scores_included = True
        except:
            pass

        self.pks = peakstructure.Peaks()
        self.pks.add_peaks(self.data.peaks, massbins=self.config.massbins, scores_included=scores_included)
        self.pks.default_params(cmap=self.config.peakcmap)

        if len(self.pks.peaks) > 0:
            ud.peaks_error_FWHM(self.pks, self.data.massdat)

        self.peaks_error_replicates(self.pks, self.data.spectra, self.config)

        for i, p in enumerate(self.pks.peaks):
            p.extracts = self.data.exgrid[i]

        self.update_history()
        self.export_params()


    def pick_scanpeaks(self):
        self.config.write_hdf5()
        self.out = metaunidec_call(self.config, "-scanpeaks")
        self.data.import_hdf5()
        self.sum_masses()

        combined_peaks = self.combine_scanpeaks()

        self.pks = peakstructure.Peaks()
        self.pks.add_peaks(combined_peaks, massbins=self.config.massbins * 0.1, scores_included=True)
        self.pks.default_params(cmap=self.config.peakcmap)

        self.scanpeaks_extracts()
        ud.peaks_error_FWHM(self.pks, self.data.massdat)
        self.peaks_error_replicates(self.pks, self.data.spectra, self.config)

        self.update_history()
        self.export_params()
        pass

    def filter_peaks(self, minscore=0.4):
        newpeaks = []
        indexes = []
        for i, p in enumerate(self.pks.peaks):
            if p.dscore > minscore:
                newpeaks.append([p.mass, p.height, p.dscore])
                indexes.append(i)
        newpeaks = np.array(newpeaks)
        indexes = np.array(indexes)
        if len(newpeaks) > 0:
            self.pks = peakstructure.Peaks()
            self.pks.add_peaks(newpeaks, massbins=self.config.massbins, scores_included=True)
            self.pks.default_params(cmap=self.config.peakcmap)

            self.data.exgrid = self.data.exgrid[indexes]
            self.normalize_exgrid()

            ud.peaks_error_FWHM(self.pks, self.data.massdat)
            self.peaks_error_replicates(self.pks, self.data.spectra, self.config)

            self.update_history()
            self.export_params()
        else:
            print("No Peaks Found with Scores Above", minscore)
            self.pks = peakstructure.Peaks()

    def combine_scanpeaks(self):
        allpeaks = np.array([])
        window = self.config.massbins * 2
        for s in self.data.spectra:
            for i, p in enumerate(s.peaks):
                if len(allpeaks) > 0:
                    masses = allpeaks[:, 0]
                    adiff = np.abs(masses - p[0])
                    if np.any(adiff < window):
                        index = np.argmin(adiff)
                        allpeaks[index, 1] += p[1]
                        weighted_dscore = (p[1] * p[2] + allpeaks[index, 1] * allpeaks[index, 2]) / (
                                p[1] + allpeaks[index, 1])
                        allpeaks[index, 2] = weighted_dscore

                        weighted_mass = (p[1] * p[0] + allpeaks[index, 1] * allpeaks[index, 0]) / (
                                p[1] + allpeaks[index, 1])
                        weighted_mass = ud.round_to_nearest(weighted_mass, self.config.massbins * 0.1)
                        allpeaks[index, 0] = weighted_mass

                    else:
                        index = len(allpeaks)
                        allpeaks = np.append(allpeaks, [p], axis=0)
                else:
                    index = 0
                    allpeaks = np.array([p])
                s.pks.peaks[i].index = index

        if self.config.peaknorm == 1:
            allpeaks[:, 1] /= np.amax(allpeaks[:, 1])

        if self.config.peaknorm == 2:
            allpeaks[:, 1] /= np.sum(allpeaks[:, 1])

        return allpeaks

    def scanpeaks_extracts(self):
        if ud.isempty(self.data.exgrid):
            print("Empty extract grid, running UniDec...")
            self.sum_masses()
            self.out = metaunidec_call(self.config, "-scanpeaks")
            self.data.import_hdf5()

        self.data.exgrid = np.zeros((len(self.pks.peaks), len(self.data.spectra)))
        for i, p in enumerate(self.pks.peaks):
            for j, s in enumerate(self.data.spectra):
                ints = 0
                for p2 in s.pks.peaks:
                    if p2.index == i:
                        ints += p2.height
                self.data.exgrid[i, j] = ints

        self.normalize_exgrid()

    def normalize_exgrid(self):
        print(self.config.exnorm)
        if self.config.exnorm == 3:
            for i in range(0, len(self.data.exgrid)):
                if np.amax(self.data.exgrid[i]) != 0:
                    self.data.exgrid[i] /= np.amax(self.data.exgrid[i])

        if self.config.exnorm == 4:
            for i in range(0, len(self.data.exgrid)):
                if np.amax(self.data.exgrid[i]) != 0:
                    self.data.exgrid[i] /= np.sum(self.data.exgrid[i])

        if self.config.exnorm == 1:
            for i in range(0, len(self.data.exgrid[0])):
                if np.amax(self.data.exgrid[:, i]) != 0:
                    self.data.exgrid[:, i] /= np.amax(self.data.exgrid[:, i])

        if self.config.exnorm == 2:
            for i in range(0, len(self.data.exgrid[0])):
                if np.amax(self.data.exgrid[:, i]) != 0:
                    self.data.exgrid[:, i] /= np.sum(self.data.exgrid[:, i])

        for i, p in enumerate(self.pks.peaks):
            p.extracts = self.data.exgrid[i]

    def peaks_heights(self):
        self.sum_masses()
        for p in self.pks.peaks:
            p.mztab = []
            p.mztab2 = []

        for i, s in enumerate(self.data.spectra):
            data2 = s.data2
            mgrid, zgrid = np.meshgrid(s.data2[:, 0], s.ztab, indexing='ij')
            mzgrid = np.transpose([np.ravel(mgrid), np.ravel(zgrid), s.mzgrid])

            mztab = ud.make_peaks_mztab(mzgrid, self.pks, self.config.adductmass, index=i)

            ud.make_peaks_mztab_spectrum(mzgrid, self.pks, data2, mztab, index=i)
        for p in self.pks.peaks:
            p.mztab = np.array(p.mztab)
            p.mztab2 = np.array(p.mztab2)

    def peaks_error_replicates(self, pks, spectra, config):
        peakvals = []
        for x in range(0, len(pks.peaks)):
            peakvals.append([])
        for i, pk in enumerate(pks.peaks):
            ints = []
            for spec in spectra:
                chopdat = ud.datachop(spec.massdat, pk.mass - config.peakwindow, pk.mass + config.peakwindow)
                if len(chopdat) > 0:
                    maxind = np.argmax(chopdat[:, 1])
                    peakvals[i].append(chopdat[maxind, 0])
                    ints.append(chopdat[maxind, 1])
                else:
                    peakvals[i].append(0)
                    ints.append(0)
            # print(peakvals[i], ints)
            pk.errorreplicate = ud.weighted_std(peakvals[i], ints)

    def export_params(self, e=None):
        peakparams = []
        for p in self.pks.peaks:
            peakparams.append([str(p.mass), str(p.height), str(p.area), str(p.label)])
        outfile = self.config.outfname + "_peaks.txt"
        try:
            np.savetxt(outfile, np.array(peakparams), delimiter=",", fmt="%s")
        except:
            pass

        peakexts = []
        for p in self.pks.peaks:
            peakexts.append(np.concatenate(([p.mass], p.extracts)))
        outfile = self.config.outfname + "_extracts.txt"
        np.savetxt(outfile, np.array(peakexts))
        print("Peak info saved to:", outfile)

    def export_spectra(self, e=None):
        for s in self.data.spectra:
            outfile = self.config.outfname + "_" + str(s.var1) + ".txt"
            np.savetxt(outfile, s.rawdata)
            print("Writing HDF5 spctrum to text file:", outfile)
            self.config.config_export(self.config.outfname + "_conf.dat")

    def batch_set_config(self, paths):
        for p in paths:
            try:
                self.config.write_hdf5(p)
                print("Assigned Config to:", p)
            except Exception as e:
                print(e)

    def batch_run_unidec(self, paths):
        for p in paths:
            try:
                tstart = time.perf_counter()
                metaunidec_call(self.config, "-all", path=p)
                print("Run:", p, " Time:  %.3gs" % (time.perf_counter() - tstart))
            except Exception as e:
                print(e)

    def batch_extract(self, paths):
        for p in paths:
            try:
                print("Extracting:", p)
                self.open(p)
                self.pick_peaks()
            except Exception as e:
                print(e)

    def fit_data(self, fit="sig"):
        print("Fitting: ", fit)
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
                print("ERROR: Unsupported fit type")
                break
            print(fits)
            self.data.fitgrid.append(fitdat)
            self.data.fits.append(fits)
        self.data.export_fits()

    def import_mzml(self, paths, timestep=1.0, scanstep=None, starttp=None, endtp=None, name=None,
                    startscan=None, endscan=None):
        """
        Tested
        :param paths:
        :param timestep:
        :return:
        """
        errors = []
        self.outpath = None
        if starttp is not None:
            self.outpath = self.parse_multiple_files(paths, starttp=starttp, endtp=endtp, timestep=timestep, name=name)
        elif startscan is not None:
            self.outpath = self.parse_multiple_files(paths, startscan=startscan, endscan=endscan, name=name)
        else:
            for p in paths:
                try:
                    if scanstep is not None:
                        self.outpath = self.parse_file(p, scanstep=scanstep)
                    else:
                        self.outpath = self.parse_file(p, timestep=float(timestep))
                except Exception as e:
                    errors.append(p)
                    print(e)
            if not ud.isempty(errors):
                print("Errors:", errors)

    def parse_file(self, p, timestep=None, scanstep=None, timepoint=None):
        """
        Tested
        :param p:
        :param timestep:
        :return:
        """
        dirname = os.path.dirname(p)
        filename = os.path.basename(p)
        print("Parsing File:", p)
        if scanstep is not None:
            self.outpath = automzml.extract_scans(filename, dirname, scanstep, "hdf5")
        else:
            self.outpath = automzml.extract(filename, dirname, timestep, "hdf5")
        return self.outpath

    def parse_multiple_files(self, paths, starttp=None, endtp=None, timestep=1.0, name=None,
                             startscan=None, endscan=None, scanstep=None):
        dirs = []
        files = []
        print("Parsing multiple files")
        for p in paths:
            dirs.append(os.path.dirname(p))
            files.append(os.path.basename(p))
        if startscan is not None:
            self.outpath = automzml.extract_scans_multiple_files(files, dirs, startscan, endscan, name)
        else:
            self.outpath = automzml.extract_timepoints(files, dirs, starttp, endtp, timestep, outputname=name)
        return self.outpath

    def csv_reader(self, csvpath):
        # Read the CSV file
        print("Reading CSV: ", csvpath)
        seq = np.genfromtxt(csvpath, delimiter=",", skip_header=1, dtype=str)
        keys = seq[0]
        kdict = {keys[i]: i for i in range(len(keys))}
        seq = seq[1:]

        # Sorting the files
        files = seq[:, kdict["File Name"]]
        dirs = seq[:, kdict["Path"]]
        csvdir = os.path.dirname(csvpath)
        # Adding the ability to check if the CSV input path is wrong
        # If it can't find the file, it tries to look in the local directory
        for i in range(len(seq)):
            f = files[i] + ".raw"
            files[i] = f
            path = os.path.join(dirs[i], f)
            if not os.path.isfile(path):
                path2 = os.path.join(csvdir, f)
                if os.path.isfile(path2):
                    dirs[i] = csvdir
                    print("Found file in local directory")
                else:
                    print("Error: Could not find file: ", path)
                    print("Error: Could not find file: ", path2)
                    print("Aborting. Please check csv file.")
                    return None
        print(files, dirs)
        # Creating the outpath
        outname = os.path.split(csvpath)[1]
        outname = os.path.splitext(outname)[0]
        print(outname)

        # Parse the file over all scans
        print("Parsing Files")
        self.outpath = automzml.extract_scans_multiple_files(files, dirs, startscan=0, endscan=-1, outputname=outname,
                                                             vars=seq, keys=kdict)
        print("Completed Parsing. Saved to", self.outpath)
        return self.outpath

    def write_to_imzML(self, outpath):
        print("Writing to imzml:", outpath)
        with ImzMLWriter(outpath) as w:
            for i, s in enumerate(self.data.spectra):
                mzs = s.massdat[:, 0]
                intensities = s.massdat[:, 1]
                if len(mzs) < 2:
                    mzs = [self.config.masslb, self.config.massub]
                    intensities = [0, 0]
                x = s.attrs['xpos']
                y = s.attrs['ypos']
                z = s.attrs['zpos']
                coords = (int(x), int(y), int(z))
                print(i, coords)
                w.addSpectrum(mzs, intensities, coords)
        print("Done Writing to imzml:", outpath)

    def imzml_to_hdf5(self, infile, outfile):
        print("Writing imzML file", infile)
        print("to HDF5 file: ", outfile)
        imzml_to_hdf5(infile, outfile)
        print("Done Writing")

    def generate_image(self, peak_index=0):
        p = self.pks.peaks[peak_index]
        ex = p.extracts
        x = np.array(self.data.var1)
        y = np.array(self.data.var2)
        dat = np.transpose([x, y, ex])
        return dat


if __name__ == '__main__':
    eng = MetaUniDec()
    '''
    testpath = "C:\Python\\unidec\\unidec_src\\unidec\\x64\Release\\test.hdf5"
    eng.data.new_file(testpath)
    data1 = [1, 2, 3]
    data2 = [4, 5, 6]
    data3 = [7, 8, 9]
    eng.data.add_data(data1)
    eng.data.add_data(data2)
    eng.data.add_data(data3)
    eng.data.remove_data([0, 2])
    exit()
    
    testdir = "C:\Python\\UniDec3\\unidec_src\\unidec\\x64\Release"
    testfile = "JAW.hdf5"
    testpath = os.path.join(testdir, testfile)
    eng.open(testpath)
    eng.run_unidec()
    eng.pick_scanpeaks()
    '''
    path = "C:\Data\HTS_Sharon\\20220404-5_sequence_Shortlist.csv"
    eng.csv_reader(path)
