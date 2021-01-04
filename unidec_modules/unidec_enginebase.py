from unidec_modules import unidecstructure, peakstructure
from copy import deepcopy
from unidec_modules import unidectools as ud
import numpy as np

version = "4.4.0"

class UniDecEngine:
    def __init__(self):
        """
        UniDec Engine Base

        Consists of three main subclasses: Config, DataContiner, Peaks
        
        Establishes some core shared functions

        :return: None
        """
        self.version = version

        print("\nUniDec Engine v." + self.version)
        self.config = None
        self.config_history = []
        self.config_count = 0
        self.data = None
        self.pks = None
        self.initialize()
        pass

    def initialize(self):
        """
        Initialize Config, DataContainer, and Peaks
        :return: None
        """
        self.config = unidecstructure.UniDecConfig()
        self.clear_history()
        self.config.initialize_system_paths()
        self.reset_config()
        self.data = unidecstructure.DataContainer()
        self.pks = peakstructure.Peaks()

    def reset_config(self):
        """
        Resets UniDec config to default. Should not reset paths or filenames.
        :return: None
        """
        self.config.initialize()
        self.update_history()

    def load_config(self, f_name):
        """
        Import UniDec Configuration File
        :param f_name: File name
        :return: None
        """
        if f_name is not None:
            self.config.config_import(f_name)
            self.update_history()
        else:
            print("Load Config Error: No file provided.")

    def export_config(self, f_name=None):
        """
        Export UniDec Configuration File
        :param f_name: File name, Default of None will using config.confname
        :return: None
        """
        self.update_history()
        if f_name is not None:
            self.config.config_export(f_name)
        else:
            self.config.config_export(self.config.confname)

    def save_default(self):
        """
        Saves existing config in default location set at self.config.defaultconfig
        :return: None
        """
        self.config.config_export(self.config.defaultconfig)

    def load_default(self):
        """
        Loads config from default location set at self.config.defaultconfig
        :return: None
        """
        self.config.config_import(self.config.defaultconfig)

    def write_hdf5(self):
        self.update_history()
        self.config.write_hdf5(self.config.hdf_file)
        # self.data.write_hdf5(self.config.hdf_file)

    def read_hdf5(self):
        self.config.read_hdf5(self.config.hdf_file)
        self.update_history()
        # self.data.read_hdf5(self.config.hdf_file)

    def update_history(self):
        # print "Update"
        try:
            if self.config_count > 0 and self.config.check_new(self.config_history[len(self.config_history) - 1]):
                self.config_history.append(self.copy_config(self.config))
                self.config_count = len(self.config_history)
                # print "Updated History", self.config_count
            elif self.config_count == 0:
                self.clear_history()
                # else:
                # print "No changes"
        except:
            self.clear_history()
        # print self.config_count
        pass

    def copy_config(self, config):
        # return deepcopy(config)
        return type("UniDecConfig", (object,), dict(config.__dict__))

    def clear_history(self):
        self.config_history = [self.copy_config(self.config)]
        self.config_count = 1
        pass

    def undo(self):
        if self.config_count > 1:
            self.config_count -= 1
            new = self.config_history[self.config_count - 1]
            old = self.config
            for item in new.__dict__:
                try:
                    old.__dict__[item] = new.__dict__[item]
                except KeyError as e:
                    print(e)
        pass

    def redo(self):
        # print(self.config_count, len(self.config_history))
        if self.config_count < len(self.config_history):
            self.config_count += 1
            new = self.config_history[self.config_count - 1]
            old = self.config
            for item in new.__dict__:
                try:
                    old.__dict__[item] = new.__dict__[item]
                except KeyError as e:
                    print(e)
        pass

    def get_auto_peak_width(self, set=True):
        try:
            fwhm, psfun, mid = ud.auto_peak_width(self.data.data2)
            self.config.automzsig = fwhm
            self.config.autopsfun = psfun
            if set:
                self.config.psfun = psfun
                self.config.mzsig = fwhm
            print("Automatic Peak Width:", fwhm)
        except Exception as e:
            print("Failed Automatic Peak Width:", e)
            print(self.data.data2)

    def check_badness(self):
        """
        Check for problematic variables, such as upper bounds less than lower bounds and raise warning if found.
        :return:
        """
        badness, warning = self.config.check_badness()
        if warning != "":
            print(warning)
        return badness

    def linear_regression_peaks(self):
        print("Starting Linear Regression using a repeating mass of:", self.config.molig)
        fit = [0, 0]
        rsquared = 0
        if self.config.molig != 0:
            x = []
            y = []
            z = []
            for p in self.pks.peaks:
                if p.ignore == 0:
                    y.append(p.mass)
                    z.append(p.height)
                    mnum = np.floor(p.mass / self.config.molig)
                    x.append(mnum)

            x = np.array(x)
            y = np.array(y)
            z = np.array(z)

            fit = np.polyfit(x, y, 1, w=z ** 2)
            slope = fit[0]
            intercept = fit[1]

            fitdat = x * slope + intercept

            sse = np.sum((fitdat - y) ** 2)
            denom = np.sum((y - np.mean(y)) ** 2)
            rsquared = 1 - ud.safedivide1(sse, denom)

            print("Slope:", fit[0], "Intercept:", fit[1], "R-Squared:", rsquared)

            residuals = np.abs(fitdat - y)
            cutoff = self.config.molig * 0.1
            boo1 = residuals > cutoff
            if np.any(boo1):
                print("Removing outliers with residuals greater than:", cutoff)
                print(residuals)
                boo2 = residuals < cutoff
                fit = np.polyfit(x[boo2], y[boo2], 1, w=np.array(z[boo2]) ** 2)
                slope = fit[0]
                intercept = fit[1]
                fitdat = x[boo2] * slope + intercept
                sse = np.sum((fitdat - y[boo2]) ** 2)
                denom = np.sum((y[boo2] - np.mean([boo2])) ** 2)
                rsquared = 1 - ud.safedivide1(sse, denom)
                print("New Slope:", fit[0], "New Intercept:", fit[1], "R-Squared:", rsquared)



        else:
            print("Need to set the mass difference/mass of oligomer")
        return fit, rsquared


    def polydispersity_index(self, e=None):
        pdi = ud.polydispersity_index(self.data.massdat)
        print(pdi)

    def oxidation_analysis(self, e=None):
        data = self.data.massdat

        maxindex = np.argmax(data[:, 1])
        maxmass, maxint = data[maxindex]
        nox = np.arange(0, 4)
        oxmasses = maxmass + nox * 16
        areas = []
        for i, m in enumerate(oxmasses):
            low = m + self.config.integratelb
            high = m + self.config.integrateub
            print("Peak:", m, "Number of Oxidations:", nox[i], "Range:", low, "to", high)
            # area, intdat = ud.integrate(data, low, high) # Peak Area
            area = ud.localmax(data[:, 1], low, high) # Peak Height
            areas.append(area)
        areas = np.array(areas)
        areas /= np.amax(areas)
        print("Relative Areas:", areas)

        totox = np.sum(nox * areas) / np.sum(areas)
        print("Total Oxidations:", totox)

        return np.append(areas, totox)
