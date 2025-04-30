import numpy as np
from scipy.interpolate import interp1d

from unidec.UniDecImporter import ImporterFactory
from unidec.UniDecImporter.Importer import Importer
from unidec.UniDecImporter.Sciex.SciexReader import SciexReader
import matplotlib.pyplot as plt
import matplotlib as mpl
from unidec import tools as ud

mpl.use('WxAgg')


class SciexImporter(Importer):
    def __init__(self, file_path, **kwargs):
        super().__init__(file_path, **kwargs)
        try:
            self.scirun = SciexReader(file_path, selected_sample = None, selected_experiment = None)
        except:
            print("Error reading Sciex file")
            print("Ensure the wiff.scan is in the same path as the wiff file")
        self.times = []
        self.primary_experiment = None
        self.primary_sample = None
        self.scan_ranges = self.get_all_scan_ranges()
        if self.primary_sample is None and self.primary_experiment is None:
            for key in self.scirun.all_scanranges.keys():
                lets = self.scirun.reader.ExperimentInfo(int(key[0]-1), int(key[1]-1))
                val = " ".join([str(i) for i in lets])
                if "MRM" not in val:
                    self.primary_experiment = int(key[1])
                    self.primary_sample = int(key[0])
                    self.scan_ranges = self.scirun.all_scanranges[key]["Scan Range"]
                    break
            print("No TOF scans found, exiting..")
            exit()
        print("Chosen Sample", self.primary_sample)
        print("Chosen Experiment", self.primary_experiment)
        print("Scan Range", [1, self.scan_ranges[1]])
        self.scans = np.arange(self.scan_ranges[0], self.scan_ranges[1] + 1)
        self.ms_level = self.get_ms_level()
        print(self.ms_level)
        print(self.polarity)
        self.imms_support = False


    def get_all_scan_ranges(self):
        if len(self.scirun.all_scanranges) == 1:
            scan_range = next(iter(self.scirun.all_scanranges.values()))["Scan Range"]
            self.primary_sample= self.get_primary_sample()
            self.primary_experiment = self.get_primary_experiment()
            for i in range(scan_range[1]):
                self.times.append(self.scirun.reader.GetRTOfScan(self.primary_sample-1, self.primary_experiment-1, i))
            return scan_range
        else:
            return list(self.scirun.all_scanranges.keys())


    def get_primary_sample(self):
        """Returns the primary sample """
        if len(self.scirun.all_scanranges) == 1:
            return  next(iter(self.scirun.all_scanranges))[0]
        return None

    def get_primary_experiment(self):
        """Returns the primary experiment"""
        if len(self.scirun.all_scanranges) == 1:
            return next(iter(self.scirun.all_scanranges))[1]
        return None

    def get_all_scans(self, threshold=-1):
        self.data = []
        for s in range(len(self.scans) - 1):
            impdat = self.scirun.reader.GetSpectrumIndex(self.primary_sample - 1, self.primary_experiment - 1,
                                                         int(s))
            imp_np = np.array(impdat)
            stacked = np.column_stack((imp_np[0, :], imp_np[1, :]))
            self.data.append(stacked)
        return self.data

    def get_avg_scan(self, scan_range=None, time_range=None):
            first_spectrum = np.array(self.scirun.reader.GetSpectrumIndex(self.primary_sample - 1,
                                                                          self.primary_experiment- 1,
                                                                          0))

            mz_values, intensity_values = first_spectrum[0, :], first_spectrum[1, :]
            common_mz = ud.nonlinear_axis(np.min(mz_values), np.max(mz_values), res=500)
            template = np.column_stack((common_mz, np.zeros_like(common_mz)))

            interp_func = interp1d(mz_values, intensity_values, kind='linear', bounds_error=False, fill_value=0)
            template[:, 1] += interp_func(common_mz)
            total_scans = self.scirun.reader.GetNumCycles(self.primary_sample - 1)


            for scan_idx in range(1, total_scans):
                spectrum = np.array(self.scirun.reader.GetSpectrumIndex(self.primary_sample - 1,
                                                                        self.primary_experiment - 1,
                                                                        scan_idx))
                if spectrum.shape[0] != 2:
                    continue

                mz_values, intensity_values = spectrum[0, :], spectrum[1, :]
                data2 = np.column_stack((mz_values, intensity_values))
                merged_spectrum = ud.mergedata(template, data2)
                template[:, 1] += merged_spectrum[:, 1]

            template[:, 1] /= total_scans
            return template

    def get_single_scan(self, idx):
        spectrum = np.array(self.scirun.reader.GetScanData(self.primary_sample - 1,
                                                                self.primary_experiment - 1,
                                                                int(idx)))

        if spectrum.shape[0] == 2:
            impdat = np.transpose([spectrum[0, :], spectrum[1, :]])
            impdat = impdat[impdat[:, 0] > 10]
            return impdat
            # return self.fill_zeros(impdat)

        return np.array([])

    def get_ms_level(self):
        info = self.scirun.reader.ExperimentInfo(self.primary_sample - 1, self.primary_experiment - 1)
        info = " ".join([str(i) for i in info]).lower()
        if "MRM" in info:
            return None
        if "ms1" in info or ("fullscan" in info and "fragmentation" not in info):
            return "MS1"
        elif "ms2" in info or "fragmentation" in info:
            return "MS2"
        else:
            # I guess here we can assume it is a MS1 scan if neither are determined
            return "MS1"

    def get_tic(self):
        dat = np.array(self.scirun.reader.TicByExp(self.primary_sample - 1, self.primary_experiment - 1))
        mz = dat[0, :]
        intensity = dat[1, :]
        return np.column_stack((mz, intensity))

    def get_polarity(self, scan = None):
        if "Positive" in self.scirun.reader.ExperimentInfo(self.primary_sample - 1, self.primary_experiment - 1):
            self.polarity = "Positive"
        else:
            self.polarity = "Negative"

    def get_mrm_chromatogram(self):
        times = []
        intensities = []
        self.scan_range = self.scirun.all_scanranges[(self.primary_sample, self.primary_experiment)]["Scan Range"]
        self.scans = np.arange(self.scan_range[0], self.scan_range[1] + 1)
        for i in range(1, len(self.scans)):
            curr_dat = np.array(self.scirun.reader.MRMScan(self.primary_sample - 1, self.primary_experiment - 1, int(self.scans[i]-1)))
            times.extend(curr_dat[:, 0])
            intensities.extend(curr_dat[:, 1])
        return np.column_stack((times, intensities))


if __name__ == '__main__':
    path = "Z:\\Group Share\\JGP\\DiverseDataExamples\\Sciex\\20221212_MAG_MPV_RA.wiff"
    importer = ImporterFactory.create_importer(path)
    dat = importer.get_avg_scan()
    for i in dat:
        print(i)


