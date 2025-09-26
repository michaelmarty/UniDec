import pyteomics.mzxml
import numpy as np
import pandas as pd
import os
from unidec.modules.thermo_reader.ThermoImporter import ThermoDataImporter
import time
import matplotlib.pyplot as plt
import matplotlib.colors as color
from copy import deepcopy
import sklearn.cluster
from scipy import stats, signal

class CDMS():
    def __init__(self, path, mzrange=None, irange=None, injection_time=0.1, resolution=240000, intensity_threshold=200, slope=12.51):
        self.path = path
        self.res = resolution
        self.it = injection_time
        self.times_FWHM = 1
        self.slope = slope
        #super duper important
        #self.slope = 17
        self.read_file()
       #self.maximum = np.maximum(mzrange[0],mzrange[1]self.filtered["Intensity"])
        self.mean = np.mean(self.filtered["Intensity"])
        self.median = np.median(self.filtered["Intensity"])
        #print("Mean: ", self.mean )
        #print("Median: ", self.median)
        self.filter_smaller_int(intensity_threshold)
        self.median_small = np.median(self.filtered_smaller["Intensity"])
        self.signaltonoise_smaller = ((np.max(self.filtered["Intensity"]))/self.median_small)
        #print(len(self.filtered))
        self.filter_int(intensity_threshold)
        #self.maximum = np

        #print(len(self.filtered))
        self.filter_centroids()
        #self.median = np.median(self.filtered["Intensity"])
        self.signaltonoise = np.array((self.filtered["Intensity"]) / self.median)
        #self.filtered["Intensity"] = self.signaltonoise
        #print(len(self.filtered))
        if mzrange is not None:
            self.select_mz(mzrange[0], mzrange[1])
            #print(len(self.filtered))
        if irange is not None:
            self.select_int(irange[0], irange[1])
            #print(len(self.filtered))
        #self.filtered
        self.convert()
        #self.signaltonoise = np.mean((self.filtered["Intensity"]) / self.median)
        self.signaltonosie_max = np.max(self.signaltonoise)
        #self.filtered["Intensity"] = self.signaltonoise
        #print(self.filtered["intensity"])
        #exit()

    def read_file(self):
        extension = os.path.splitext(self.path)[1]
        starttime = time.perf_counter()
        if extension.lower() == ".raw":
            self.TDI = ThermoDataImporter(self.path)
            self.data = self.TDI.grab_data()
            self.scans = self.TDI.scans
            self.it = self.TDI.msrun.injection_time
            self.res = self.TDI.msrun.resolution
            scans = np.concatenate([s * np.ones(len(self.data[i])) for i, s in enumerate(self.scans)])
            mz = np.concatenate([d[:, 0] for d in self.data])
            intensity = np.concatenate([d[:, 1] for d in self.data]) * self.it / 1000.
            self.normalized = pd.DataFrame({"m/z": mz, "Intensity": intensity, "Scan": scans})
            self.filtered = self.normalized
        if extension.lower() == ".mzxml":
            data = pyteomics.mzxml.read(self.path)
            scans = [i for i in data.iterfind("*") if "basePeakIntensity" in i]
            mz = [s["m/z array"] for s in scans]
            intensity = [s["intensity array"] for s in scans]
            scans = [[int(s["num"])] * s["peaksCount"] for s in scans]
            mz = np.concatenate(mz)
            intensity = np.concatenate(intensity) * self.it / 1000.
            scans = np.concatenate(scans)
            self.normalized = pd.DataFrame({"m/z": mz, "Intensity": intensity, "Scan": scans})
            self.filtered = self.normalized
        self.unfiltered = deepcopy(self.filtered)  # copy unfiltered data for later comparison
        print("Load Time:", time.perf_counter() - starttime)

    def select_mz(self, mi, ma):
        '''
        Function filtering centroids based on m/z values.
        '''
        self.filtered = self.filtered[self.filtered["m/z"].between(mi, ma)]

    def select_int(self, mi, ma):
        '''
        Function filtering centroids based on intensity values.
        '''
        self.filtered = self.filtered[self.filtered["Intensity"].between(mi, ma)]

    def filter_int(self, int_threshold):
        self.filtered = self.filtered[self.filtered["Intensity"] >= int_threshold]

    def filter_smaller_int(self, int_threshold):
        self.filtered_smaller = self.filtered[self.filtered["Intensity"] < int_threshold]
        #print(self.filtered_smaller, int_threshold)

    def filter_centroids(self, mz_threshold="auto"):
        select = [self.__remove_adjacent(df["m/z"], mz_threshold) for s, df in self.filtered.groupby(["Scan"])]
        self.filtered = self.filtered[np.concatenate(select)]

    #def filter_centroids_smaller(self, mz_threshold= "auto"):
     #   select = [self.__remove_adjacent(df["m/z"], mz_threshold) for s, df in self.filtered.groupby(["Scan"])]
      #  self.filtered = self.filtered[np.concatenate(select)]

    def __remove_adjacent(self, mz, mz_threshold):
        '''
        Function reads list of m/z centroids and filters for adjacent centroids
        below a certain m/z threshold. The m/z threshold can either be set to a fixed value or
        determined automatically based on the FWHM at the instruments resolution limit.
        General procedure is:
        1. Input MZ_centrois_list if legth n, MZ_threshold
            [1,2,3,3.5,4], 1
        2. Calculate difference between items yielding list with n-1 items
            [1,2,3,3.5,4] -> [1,1,0.5,0.5]
        3. Apply threshold (T=True, F=False)
            [1,1,0.5,0.5] >= 1 -> [T, T, F, F]
        4. In order to get a Boolean list for excluding both peaks which distance is below the threshold we
        generate two Boolean lists with the length n by duplicating the first and the last item, respectively, and
        perform a Boolean conjunction.
            [T, T, F, F] + [F] -> [T, T, F, F, F]
            [T] + [T, T, F, F] -> [T, T, T, F, F]
            [T, T, F, F, F] AND [T, T, T, F, F] -> [T, T, F, F, F]
        5. This Boolean list will be used as selector for the provided list of centroids.
            [T, T, F, F, F] -> [1,2,3,3.5,4] -> [1,2]

        If the list of centroids contains only one item the function will raise an indexing error in step 4
        and the function will return a [False]

        '''
        if mz_threshold == "auto":
            try:
                mz_array = np.array(mz)
                diff_b = np.diff(mz_array) >= (self.__get_FWHM(mz_array[:-1], self.res) * self.times_FWHM)
                diff_b = np.append(diff_b[0], diff_b) & np.append(diff_b, diff_b[-1])
            except:
                diff_b = np.array([False] * len(mz_array))
        else:
            try:
                mz_array = np.array(mz)
                diff_b = np.diff(mz_array) >= mz_threshold
                diff_b = np.append(diff_b[0], diff_b) & np.append(diff_b, diff_b[-1])
            except:
                diff_b = np.array([False] * len(mz_array))
        return diff_b

    def __get_FWHM(self, mz, res):
        '''
        Calculates FWHM of a peaks based on instrument resolution at 400 m/z, assuming R ~ m/z^-0.5
        '''
        A = res / 400 ** -.5
        R = A * mz ** -0.5
        return mz / R
        self.unclustered = deepcopy(self.filtered)

    def cluster(self):
        self.data2 = np.transpose([cd.filtered["m/z"], np.array(cd.filtered["z"])])
        #self.data2 = data2[data2[:, 0].argsort()]
        self.clusters = sklearn.cluster.DBSCAN(eps=2, min_samples=5, metric="euclidean").fit(self.data2)
        self.clabels = np.unique(self.clusters.labels_)
        self.filtered["cz"] = self.filtered["z"]
        for i in self.clabels:
            boo1 = self.clusters.labels_ == i
            if i == -1:
                self.filtered["cz"][boo1] = -1
            else:

                means = np.mean(self.data2[boo1], axis=0)
                means[1] = round(means[1])
                self.filtered["cz"][boo1] = means[1]
        self.filtered = self.filtered[self.filtered["cz"].between(0, np.amax(self.filtered["cz"]))]
        #print(len(self.filtered))

    def convert(self):
        self.filtered["z"] = self.filtered["Intensity"] / self.slope
        self.filtered["mass"] = self.filtered["m/z"] * self.filtered["z"] - self.filtered["z"]

if __name__ == '__main__':
    path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\03112021\\20210311_MK_GroEL_CDMS_10min_NT0.RAW"
    mzrange = [9000, 13000]
    irange = [1, 1200]
    #path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_BSA__CD-MS_Aqu_ACN_10min.RAW"
    #mzrange = [3000, 7000]
    #irange = [100, 600]
    #path = "Z:\\Group Share\\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_AMTB_100V_CD-MS_8min_aqu_15ms.RAW"
    #mzrange = [4000, 12000]
    #irange = [100, 900]
    #path = "Z:\\Group Share\Marius Kostelic\\CD-MS\\10192020 DPPC\\201019_mk_cz_DPPC_24hrs_pos_singleion_25injecttime_pressure05.raw"
    #mzrange = [7500, 15000]
    #irange = [50, 600]
    #path = "Z:\\Group Share\Marius Kostelic\\CD-MS\\03102021\\20210309_MK_GroEL_CDMS_prseeure07_10min_aqu.RAW"
    #mzrange = [9000, 13000]
    #irange = [100, 2000]



    cd = CDMS(path, mzrange=mzrange, irange=irange)

    zrange = np.array(irange) / cd.slope

    custom_map = color.LinearSegmentedColormap.from_list("map", ["white", "red", "darkred"])
    bins = (np.arange(mzrange[0], mzrange[1], 50), np.arange(irange[0], irange[1], 20))
    bins2 = (np.arange(mzrange[0], mzrange[1], 50), np.arange(zrange[0], zrange[1], 1))
    fix, ax = plt.subplots(ncols=2, figsize=(8, 3))
    ax[1].hist2d(cd.filtered["m/z"], cd.filtered["z"], bins=bins2, cmap=custom_map, cmin=1, )
    cd.cluster()
    ax[0].hist2d(cd.filtered["m/z"], cd.filtered["cz"], bins=bins2, cmap=custom_map, cmin=1, )

    if False:
        for i in cd.clabels:
            if i == -1:
                pass
            else:
                boo1 = cd.clusters.labels_ == i
                means = np.mean(cd.data2[boo1], axis=0)
                means[1] = round(means[1])
                print(means, means[0] * means[1] - means[1])
                ax[0].scatter(cd.data2[boo1, 0], cd.data2[boo1, 1])
                ax[0].scatter(cd.filtered["m/z"][boo1], cd.filtered["cz"][boo1])
                ax[0].scatter(means[0], means[1], color="k")

    #unclustered_data = filtered_data.copy()
    #unclustered_data = pd.DataFrame()
    #unclustered_data["z"] = filtered_data.copy["z"]
    #unclustered_data["mass"] = unclustered_data["m/z"] * unclustered_data["z"] - unclustered_data["z"]

    filtered_data=pd.DataFrame()
    filtered_data["z"] = cd.filtered["cz"]
    filtered_data["mass"] = cd.filtered["m/z"] * filtered_data["z"] - filtered_data["z"]

    bins = np.arange(1e4, 3e5, 1e2)
    fix, ax = plt.subplots(figsize=(3, 3))

    #actualdata = np.loadtxt("Z:\\Group Share\Marius Kostelic\\CD-MS\\03082021\\20210308_MK_GroEL_negative_mode_5min_aqu_2_unidecfiles\\20210308_MK_GroEL_negative_mode_5min_aqu_2_mass.txt")
    #actualdata = np.loadtxt("Z:\\Group Share\Marius Kostelic\\CD-MS\\02242021\\02242021_MK_AMTB_100V_1min_qu_unidecfiles\\02242021_MK_AMTB_100V_1min_qu_mass.txt")
    #actualdata = np.loadtxt("Z:\\Group Share\\Marius Kostelic\\CD-MS\\09302020\\20200930_MK_DPPG_pos_10x_dilution_unidecfiles\\20200930_MK_DPPG_pos_10x_dilution_mass.txt")
    #actualdata = np.loadtxt("Z:\\Group Share\\Marius Kostelic\\CD-MS\\02192021 AMTB\\20210219_MMK_CDMS_AMTB_C8E4_1MIN_AQU_unidecfiles\\20210219_MMK_CDMS_AMTB_C8E4_1MIN_AQU_mass.txt")
    actualdata = np.loadtxt("Z:\\Group Share\Marius Kostelic\\CD-MS\\03042021\\210304_MK_ADH_pos_1to100_pos_1min_aqu_unidecfiles\\210304_MK_ADH_pos_1to100_pos_1min_aqu_mass.txt")
    #load actual data
    boo1 = actualdata[:,0]>bins[0]
    plt.plot(actualdata[boo1,0], actualdata[boo1,1]/np.amax(actualdata[boo1,1])*40)

    c, b, f = plt.hist(filtered_data["mass"], color="k", lw=1, bins=bins)
    kde = stats.gaussian_kde(filtered_data["mass"])
    xx = bins
    yy = kde(xx)
    n = c.sum() / yy.sum()

    #c, b, f = plt.hist(unclustered_data["mass"], color="r", lw=1, bins=bins)
    #kde = stats.gaussian_kde(unclustered_data["mass"])
    #xx = bins
    #yy = kde(xx)
    #n = c.sum() / yy.sum()

    peak_idx = signal.find_peaks(yy / yy.max(), height=0.1)[0]

    #ax.plot([xx[peak_idx]] * 2, [0, yy[peak_idx] * n], color="r", label="%d kDa" % (xx[peak_idx] / 1e3))
    #ax.plot(xx, yy * n, "k", lw=1)




    plt.legend(frameon=False)
    plt.title("Mass histogram\n 10 Da bins")
    #ax.set_xticks(np.arange(2.5e6, 3.75e6, 2.5e5))
    #ax.set_xticklabels(np.arange(2.5, 3.75, 0.25))
    plt.xlabel("Mass [kDa]")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("AaLS-neg_MassHistogram.png")
    plt.show()

    plt.show()
