import pyteomics.mzxml
import numpy as np
import pandas as pd


class Preprocess():
    '''
    Class for pre-processing single particle data.
    '''

    def __init__(self, file, it=1, res=220000, times_FWHM=1):
        self.file = file
        self.it = it
        self.res = 240000
        self.times_FWHM = 1
        self.__read_peaks()

    def __read_peaks(self):
        '''
        Function reads mzXML file and normalizes the intensities by their injection time and saves centroids
        with the scan number in a DataFrame.
        '''
        data = pyteomics.mzxml.read(self.file)
        scans = [i for i in data.iterfind("*") if "basePeakIntensity" in i]
        mz = [s["m/z array"] for s in scans]
        intensity = [s["intensity array"] for s in scans]
        scans = [[int(s["num"])] * s["peaksCount"] for s in scans]
        mz = np.concatenate(mz)
        intensity = np.concatenate(intensity) * self.it
        scans = np.concatenate(scans)
        self.normalized = pd.DataFrame({"m/z": mz, "Intensity": intensity, "Scan": scans})
        self.filtered = self.normalized

    def filter_data(self, int_threshold, mz_threshold="auto"):
        '''
        Function filtering centroids data for noise and dephased peaks.
        '''
        self.filtered = self.filtered[self.filtered["Intensity"] >= int_threshold]
        select = [self.__remove_adjacent(df["m/z"], mz_threshold) for s, df in self.filtered.groupby(["Scan"])]
        self.filtered = self.filtered[np.concatenate(select)]

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

    def write_peaks(self, label="_processed"):
        '''
        Writes filtered centroids as CSV list.
        '''
        self.filtered.to_csv(self.file + label + ".csv", index=False)

    def __get_FWHM(self, mz, res):
        '''
        Calculates FWHM of a peaks based on instrument resolution at 400 m/z, assuming R ~ m/z^-0.5
        '''
        A = res / 400 ** -.5
        R = A * mz ** -0.5
        return mz /R

# imports for making the plots in the example
import matplotlib.pyplot as plt
import matplotlib.colors as color
from scipy import stats, signal


custom_map = color.LinearSegmentedColormap.from_list("map", ["white", "red", "darkred"])

AaLS  = Preprocess("Z:\\Group Share\\Marius Kostelic\\CD-MS\\10192020 DPPC\\201019_mk_BSA_pos_singleion_longaquisition.mzXML", it=0.1) # reading mzXML file
#AaLS  = Preprocess("Z:\\Group Share\\Marius Kostelic\\CD-MS\\09302020\\20200930_MK_DPPG_pos_10x_dilution_singleion_25msinject_512transient.mzXML", it=0.025) # reading mzXML file
unfiltered = AaLS.filtered.copy() # copy unfiltered data for later comparison

AaLS.filter_data(300)
mzrange=[2000,8000]
irange=[50,1500]
#set the range and filter threshold^^
AaLS.select_mz(mzrange[0], mzrange[1]) # selects all centroids between 19000 and 21500 m/z
AaLS.select_int(irange[0], irange[1]) # selects all centroids between 1500 and 2500 a.u. intensity
data =np.transpose([AaLS.filtered["m/z"], AaLS.filtered["Intensity"]])

#print(data)
#plt.hist(data[:,1], bins=500)
#plt.show()
#exit()

bins = (np.arange(mzrange[0], mzrange[1], 50), np.arange(irange[0], irange[1], 20))
fix, ax = plt.subplots(ncols=2, figsize=(8,3))
ax[0].hist2d(unfiltered["m/z"], unfiltered["Intensity"], bins=bins,  cmap=custom_map,cmin=1,)
ax[1].hist2d(AaLS.filtered["m/z"], AaLS.filtered["Intensity"], bins=bins, cmap=custom_map,cmin=1)


ax[0].set_title("before filtering")
ax[0].set_ylabel("Intensity")
ax[0].set_xlabel("m/z")
ax[1].set_title("after filtering")
ax[1].set_xlabel("m/z")
plt.tight_layout()
#plt.show()
#plt.savefig"Filtering_effect.png")



#AaLS.write_peaks()

#filtered_data = pd.read_csv("./data/AaLS.mzXML_processed.csv")
slope=12.51
#enter the charge multiplier here^
filtered_data=pd.DataFrame()
filtered_data["z"] = AaLS.filtered["Intensity"] / slope
filtered_data["mass"] = AaLS.filtered["m/z"] * filtered_data["z"] - filtered_data["z"]

bins = np.arange(1e4, 3e5, 1e3)
fix, ax = plt.subplots(figsize=(3, 3))

c, b, f = plt.hist(filtered_data["mass"], color="lightgrey", lw=0.75, bins=bins)

actualdata = np.loadtxt("Z:\\Group Share\\Marius Kostelic\\CD-MS\\10192020 DPPC\\201019_mk_BSA_pos_actual.txt")
#actualdata = np.loadtxt("Z:\\Group Share\\Marius Kostelic\\CD-MS\\09302020\\20200930_MK_DPPG_pos_10x_dilution_unidecfiles\\20200930_MK_DPPG_pos_10x_dilution_mass.txt")
#load actual data
boo1 = actualdata[:,0]>bins[0]
plt.plot(actualdata[boo1,0], actualdata[boo1,1]/np.amax(actualdata[boo1,1])*np.amax(c))


kde = stats.gaussian_kde(filtered_data["mass"])
xx = bins
yy = kde(xx)
n = c.sum() / yy.sum()

peak_idx = signal.find_peaks(yy / yy.max(), height=0.1)[0]

#ax.plot([xx[peak_idx]] * 2, [0, yy[peak_idx] * n], color="r", label="%d kDa" % (xx[peak_idx] / 1e3))
ax.plot(xx, yy * n, "k", lw=0.75)

plt.legend(frameon=False)
plt.title("Mass histogram\n 10 Da bins")
#ax.set_xticks(np.arange(2.5e6, 3.75e6, 2.5e5))
#ax.set_xticklabels(np.arange(2.5, 3.75, 0.25))
plt.xlabel("Mass [kDa]")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("AaLS-neg_MassHistogram.png")
plt.show()



