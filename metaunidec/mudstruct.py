import numpy as np
from unidec_modules.hdf5_tools import replace_dataset, get_dataset
import h5py
import unidec_modules.unidectools as ud
import os
from copy import deepcopy
from unidec_modules.peakstructure import Peaks

class MetaDataSet:
    def __init__(self, engine):
        self.names = []
        self.indexes = []
        self.spectra = []
        self.topname = "ms_dataset"
        self.filename = None
        self.massdat = []
        self.data2 = []
        self.massgrid = []
        self.mzgrid = []
        self.mzdat = []
        self.exgrid = []
        self.exvals = []
        self.peaks = []
        self.var1 = []
        self.var2 = []
        self.v1name = "Variable 1"
        self.v2name = "Variable 2"
        self.len = 0
        self.eng = engine
        self.fitgrid = []
        self.fits = []
        pass

    def import_hdf5(self, file=None):
        if file is None:
            file = self.filename
        else:
            self.filename = file
        hdf = h5py.File(file, 'r')
        msdata = hdf.require_group(self.topname)
        keys = list(msdata.keys())
        self.indexes = []
        for k in keys:
            try:
                self.indexes.append(int(k))
            except:
                pass
        self.indexes = np.array(self.indexes)
        self.indexes = sorted(self.indexes)
        self.len = len(self.indexes)
        hdf.close()

        if ud.isempty(self.spectra):
            for i in self.indexes:
                s = Spectrum(self.topname, i, self.eng)
                s.read_hdf5(file)
                self.spectra.append(s)
        else:
            for s in self.spectra:
                s.read_hdf5(file)
        if not ud.isempty(self.spectra):
            self.data2 = self.spectra[0].data2
            self.import_vars()

    def export_hdf5(self, file=None, vars_only=False, delete=False):
        if file is None:
            file = self.filename
        else:
            self.filename = file

        if file is None:
            print("Error: No HDF5 file present. Please create an HDF5 file first")
            return

        # Clear Group
        hdf = h5py.File(file, 'a')
        if delete:
            try:
                del hdf["/" + self.topname]
                pass
            except:
                pass
        group = hdf.require_group(self.topname)
        group.attrs["num"] = len(self.spectra)
        group.attrs["v1name"] = self.v1name
        group.attrs["v2name"] = self.v2name
        config = hdf.require_group("config")
        config.attrs["metamode"] = -1
        hdf.close()

        self.var1 = []
        self.var2 = []
        for i, s in enumerate(self.spectra):
            s.index = str(i)
            s.attrs[self.v1name] = s.var1
            s.attrs[self.v2name] = s.var2
            s.attrs["name"] = s.name
            self.var1.append(s.var1)
            self.var2.append(s.var2)
            s.write_hdf5(self.filename, vars_only=vars_only)
        self.var1 = np.array(self.var1)
        # print("Variable 1:", self.var1)
        self.var2 = np.array(self.var2)
        self.len = len(self.spectra)

    def export_vars(self, file=None):
        for s in self.spectra:
            if s.ignore == 1:
                return
        self.export_hdf5(vars_only=True)

    def import_grids(self):
        hdf = h5py.File(self.filename, 'r')
        msdataset = hdf.require_group(self.topname)
        # Mass Space
        axis = get_dataset(msdataset, "mass_axis")
        sum = get_dataset(msdataset, "mass_sum")
        grid = get_dataset(msdataset, "mass_grid")
        self.massdat = np.transpose([axis, sum])
        try:
            num = int(len(grid) / len(sum))
            grid = grid.reshape((int(num), len(sum))).transpose()
            grid2 = np.transpose([axis for n in range(num)])
            grid3 = [grid2, grid]
            self.massgrid = np.transpose(grid3)
        except:
            print("Mass Grid Warning:", grid.shape, sum.shape)  # , len(grid)/len(sum))
            num = 0
        # MZ Space
        axis = get_dataset(msdataset, "mz_axis")
        sum = get_dataset(msdataset, "mz_sum")
        grid = get_dataset(msdataset, "mz_grid")
        self.mzdat = np.transpose([axis, sum])
        try:
            num = int(len(grid) / len(sum))
            grid = grid.reshape((int(num), len(sum))).transpose()
            grid2 = np.transpose([axis for n in range(num)])
            grid3 = [grid2, grid]
            self.mzgrid = np.transpose(grid3)
        except:
            print("mz grid Warning:", grid.shape, sum.shape)  # , len(grid)/len(sum))
        hdf.close()
        return num

    def import_peaks(self):
        hdf = h5py.File(self.filename, 'r')
        pdataset = hdf.require_group("/peaks")
        self.peaks = get_dataset(pdataset, "peakdata")
        self.exgrid = get_dataset(pdataset, "extracts").transpose()
        self.exvals = self.peaks[:, 0]
        hdf.close()

    def export_fits(self):
        hdf = h5py.File(self.filename, 'r+')
        pdataset = hdf.require_group("/peaks")
        replace_dataset(pdataset, "fit_extracts", self.fitgrid)
        replace_dataset(pdataset, "fits", self.fits)
        hdf.close()

    def import_grids_and_peaks(self):
        if len(self.spectra) > 0:
            self.eng.make_grids()
            num = self.import_grids()
            self.import_peaks()
            if len(self.spectra) != num:
                self.export_hdf5()
                # self.eng.run_unidec()
                self.eng.make_grids()
                self.import_grids()
                self.import_peaks()

    def new_file(self, path):
        self.__init__(self.eng)
        if os.path.isfile(path):
            os.remove(path)
        self.filename = path
        hdf = h5py.File(self.filename, 'w')
        msdata = hdf.require_group(self.topname)
        config = hdf.require_group("config")
        config.attrs["metamode"] = -1
        hdf.close()

    def add_file(self, filename=None, dirname=None, path=None):
        if path is None:
            path = os.path.join(dirname, filename)
        else:
            dirname, filename = os.path.split(path)

        data = ud.load_mz_file(path)
        self.add_data(data, name=filename)

    def add_data(self, data, name="", attrs=None, export=True):
        snew = Spectrum(self.topname, self.len, self.eng)
        snew.rawdata = data
        snew.data2 = deepcopy(data)
        if self.eng.config.datanorm == 1:
            try:
                snew.data2[:, 1] /= np.amax(snew.data2[:, 1])
            except:
                pass

        snew.name = name
        if attrs is not None:
            snew.attrs = attrs

        try:
            if "CID" in name or "SID" in name:
                print(name)
                name = os.path.splitext(name)[0]
                splits = name.split("_")
                for i, s in enumerate(splits):
                    if s == "CID" or s == "SID":
                        try:
                            snew.var1 = float(splits[i + 1])
                        except Exception as e:
                            print(splits, e)

                        self.v1name = "Collision Voltage"
        except Exception as e:
            print("Error in parsing collision voltage:", e)
            print(name)

        self.spectra.append(snew)
        self.len = len(self.spectra)
        if export:
            self.export_hdf5()

    def remove_data(self, indexes):
        for i in sorted(indexes, reverse=True):
            del self.spectra[i]
        self.export_hdf5(delete=True)
        print("Removed")

    def clear(self):
        self.spectra = []
        self.export_hdf5(delete=True)

    def get_spectra(self):
        spectra = []
        for i, s in enumerate(self.spectra):
            if s.ignore == 0:
                spectra.append(s)
        return spectra

    def safemax(self, s):
        try:
            maxval = np.amax(s.data2[:, 1])
        except:
            maxval = 0
        return maxval

    def get_max_data2(self):
        maxes = [self.safemax(s) for s in self.spectra]
        try:
            maxval = np.amax(maxes)
        except:
            maxval = 0
        return maxval


    def get_bool(self):
        bool_array = []
        for i, s in enumerate(self.spectra):
            if s.ignore == 0:
                bool_array.append(True)
            else:
                bool_array.append(False)
        return np.array(bool_array)


    def import_vars(self, get_vnames=True):
        hdf = h5py.File(self.filename, 'r+')
        msdata = hdf.require_group(self.topname)
        if get_vnames:
            try:
                self.v1name = ud.smartdecode(msdata.attrs["v1name"])
            except:
                msdata.attrs["v1name"] = self.v1name
            try:
                self.v2name = ud.smartdecode(msdata.attrs["v2name"])
            except:
                msdata.attrs["v2name"] = self.v2name
        else:
            msdata.attrs["v1name"] = str(self.v1name)
            msdata.attrs["v2name"] = str(self.v2name)
        hdf.close()

        self.update_var_array()


    def update_var_array(self):
        self.var1 = []
        self.var2 = []
        for i in range(0, len(self.spectra)):
            s = self.spectra[i]
            try:
                s.var1 = s.attrs[self.v1name]
            except:
                s.var1 = i
            try:
                s.var2 = s.attrs[self.v2name]
            except:
                s.var2 = 0
            try:
                s.name = s.attrs["name"]
            except:
                s.name = ""
            if s.var1 is None or s.var1 == "None":
                s.var1 = i
            self.var1.append(s.var1)
            self.var2.append(s.var2)
        self.var1 = np.array(self.var1)
        self.var2 = np.array(self.var2)
        try:
            self.var1 = self.var1.astype(float)
        except:
            pass
        try:
            self.var2 = self.var2.astype(float)
        except:
            pass
        # print("Variable 1:", self.var1)


class Spectrum:
    def __init__(self, topname, index, eng):
        # self.fitdat = np.array([])
        # self.baseline = np.array([])
        # self.fitdat2d = np.array([])
        self.rawdata = np.array([])
        self.data2 = np.array([])
        self.massdat = np.array([])
        self.mzgrid = np.array([])
        self.massgrid = np.array([])
        self.ztab = np.array([])
        self.zdata = np.array([])
        self.peaks = np.array([])
        self.pks = Peaks()
        self.index = index
        self.name = ""
        self.topname = topname
        self.filename = ""
        self.attrs = {}
        self.color = "k"
        self.ignore = 0
        self.var1 = 0
        self.var2 = 0
        self.eng = eng

    def write_hdf5(self, file=None, vars_only=False):
        if file is None:
            file = self.filename
        else:
            self.filename = file

        hdf = h5py.File(file, 'a')
        msdata = hdf.require_group(self.topname + "/" + str(self.index))
        if not vars_only:
            replace_dataset(msdata, "raw_data", self.rawdata.astype(self.eng.config.dtype))
            # print(self.eng.config.dtype)
            # replace_dataset(msdata, "fit_data", self.fitdat)
            replace_dataset(msdata, "processed_data", self.data2.astype(self.eng.config.dtype))
            replace_dataset(msdata, "mass_data", self.massdat.astype(self.eng.config.dtype))
            replace_dataset(msdata, "mz_grid", self.mzgrid.astype(self.eng.config.dtype))
            replace_dataset(msdata, "mass_grid", self.massgrid.astype(self.eng.config.dtype))
            # replace_dataset(msdata, "baseline", self.baseline)
            replace_dataset(msdata, "charge_data", self.zdata.astype(self.eng.config.dtype))
        for key, value in list(self.attrs.items()):
            msdata.attrs[key] = value
        hdf.close()

    def read_hdf5(self, file=None):
        if file is None:
            file = self.filename
        else:
            self.filename = file
        hdf = h5py.File(file, 'r')
        msdata = hdf.get(self.topname + "/" + str(self.index))
        self.rawdata = get_dataset(msdata, "raw_data")
        # self.fitdat = get_dataset(msdata, "fit_data")
        self.data2 = get_dataset(msdata, "processed_data")
        if ud.isempty(self.data2) and not ud.isempty(self.rawdata):
            self.data2 = deepcopy(self.rawdata)
        self.massdat = get_dataset(msdata, "mass_data")
        self.zdata = get_dataset(msdata, "charge_data")
        if self.eng.config.datanorm == 1:
            try:
                self.data2[:, 1] /= np.amax(self.data2[:, 1])
            except:
                pass
            try:
                self.massdat[:, 1] /= np.amax(self.massdat[:, 1])
            except:
                pass
            try:
                self.zdata[:, 1] /= np.amax(self.zdata[:, 1])
            except:
                pass
        self.mzgrid = get_dataset(msdata, "mz_grid")
        self.massgrid = get_dataset(msdata, "mass_grid")
        try:
            self.ztab = self.zdata[:, 0]
        except:
            pass
        # self.baseline = get_dataset(msdata, "baseline")
        self.peaks = get_dataset(msdata, "peaks")
        self.setup_peaks()
        self.attrs = dict(list(msdata.attrs.items()))
        hdf.close()

    def setup_peaks(self):
        self.pks.add_peaks(self.peaks, scores_included=True)
        self.pks.default_params()
        #print([[p.mass, p.height, p.dscore] for p in self.pks.peaks])

