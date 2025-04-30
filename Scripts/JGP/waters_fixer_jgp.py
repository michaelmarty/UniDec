import unidec.tools as ud
import matplotlib.pyplot as plt
import time
import matplotlib as mpl
from pyopenms import *
import os
import fnmatch
from copy import deepcopy
from unidec.UniDecImporter import ImporterFactory
from unidec.UniDecImporter.ImportTools import *


def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

def datachop(datatop, newmin, newmax):
    """
    Chops the range of the data. The [:,0] column of the data is the column that will be indexed.
    :param datatop: Data array
    :param newmin: Minimum value of chopped data
    :param newmax: Maximum value of chopped data
    :return: New data limited between the two bounds
    """
    boo1 = np.logical_and(datatop[:, 0] <= newmax, datatop[:, 0] >= newmin)
    return datatop[boo1]

def fix_parent_mz(exp, i, p, window=2, purity_thresh=0.5):
    j = deepcopy(i)
    isMS1 = False
    while not isMS1:
        if exp[j].getMSLevel() == 1:
            isMS1 = True
        else:
            j -= 1
    ms1 = exp[j]

    pmass = p.getMZ()
    mz, intensity = ms1.get_peaks()
    data = np.transpose([mz, intensity])
    chopdat = datachop(data, pmass - window, pmass + window)
    if len(chopdat) > 0:
        index = np.argmax(chopdat[:, 1])
        peakmz = chopdat[index, 0]
        peakint = chopdat[index, 1]
        #print(pmass, peakmz, peakint)

        p.setMZ(peakmz)
        p.setIntensity(peakint)
        p.setCharge(1)

        p.setIsolationWindowLowerOffset(window)
        p.setIsolationWindowUpperOffset(window)

        if purity_thresh is not None:
            purity_score = PrecursorPurity().computePrecursorPurity(ms1, p, 10, True)
            prop = purity_score.signal_proportion
            #print("Proportion:", prop)
            if prop < purity_thresh:
                return None

        polarity = ms1.getInstrumentSettings().getPolarity()
        if polarity == IonSource.Polarity.POSITIVE:
            p.setCharge(1)
        elif polarity == IonSource.Polarity.NEGATIVE:
            p.setCharge(-1)

        return p
    else:
        return None

def adjust_mz(s, delta):
    mz, intensity = s.get_peaks()
    mz2 = mz + delta
    s.set_peaks([mz2, intensity])
    return s

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

def datachop(datatop, newmin, newmax):
    """
    Chops the range of the data. The [:,0] column of the data is the column that will be indexed.
    :param datatop: Data array
    :param newmin: Minimum value of chopped data
    :param newmax: Maximum value of chopped data
    :return: New data limited between the two bounds
    """
    boo1 = np.logical_and(datatop[:, 0] <= newmax, datatop[:, 0] >= newmin)
    return datatop[boo1]

def fix_parent_mz(exp, i, p, window=2, purity_thresh=0.5):
    j = deepcopy(i)
    isMS1 = False
    while not isMS1:
        if exp[j].getMSLevel() == 1:
            isMS1 = True
        else:
            j -= 1
    ms1 = exp[j]

    pmass = p.getMZ()
    mz, intensity = ms1.get_peaks()
    data = np.transpose([mz, intensity])
    chopdat = datachop(data, pmass - window, pmass + window)
    if len(chopdat) > 0:
        index = np.argmax(chopdat[:, 1])
        peakmz = chopdat[index, 0]
        peakint = chopdat[index, 1]
        #print(pmass, peakmz, peakint)

        p.setMZ(peakmz)
        p.setIntensity(peakint)
        p.setCharge(1)

        p.setIsolationWindowLowerOffset(window)
        p.setIsolationWindowUpperOffset(window)

        if purity_thresh is not None:
            purity_score = PrecursorPurity().computePrecursorPurity(ms1, p, 10, True)
            prop = purity_score.signal_proportion
            #print("Proportion:", prop)
            if prop < purity_thresh:
                return None

        polarity = ms1.getInstrumentSettings().getPolarity()
        if polarity == IonSource.Polarity.POSITIVE:
            p.setCharge(1)
        elif polarity == IonSource.Polarity.NEGATIVE:
            p.setCharge(-1)

        return p
    else:
        return None

def get_mzs_ppmerrors(indexed_file, calibrants, mz_width=0.3, rt_width=6,
                      min_peak_intensity=1000, scans_to_check=5):
    expected_mzs_ppmerrors = [(0,0)]

    for c in calibrants:
        min_rt = c[1] - (rt_width / 2)
        max_rt = c[1] + (rt_width / 2)
        xic = indexed_file.extract_xic(c[0], mz_width, [min_rt, max_rt])
        max_intensity_index = np.argmax(xic[:,1])
        max_intensity = xic[max_intensity_index, 1]

        exp_mzs = []

        if max_intensity > min_peak_intensity:
            xic_ordered_by_intensity = sorted(xic, key=lambda x:x[1], reverse=True)
            for i in range(scans_to_check):
                if xic_ordered_by_intensity[i][1] > min_peak_intensity:
                    scan_at_rt = indexed_file.get_indexed_spectrum_atRt(xic_ordered_by_intensity[i][0])
                    mz = scan_at_rt.get_exp_mz(c[0], mz_width)
                    exp_mzs.append(mz)
        else:
            print("Calibration molecule at mz: " + str(c[0]) + " not found.")
            plt.plot(xic[:,0], xic[:,1])
            plt.show()

        avg_exp_mz = np.mean(exp_mzs)
        ppm_error = ((c[0]-avg_exp_mz)/c[0])*1e6
        expected_mzs_ppmerrors.append((c[0], ppm_error))

    ordered_mzs_ppmerrors = np.transpose(np.vstack(sorted(expected_mzs_ppmerrors, key=lambda x:x[0])))
    return ordered_mzs_ppmerrors

def calibrate_spectrum_mzs(mzs, indexed_mzcals, min_mz, bin_width):
    bin_indices = np.floor((mzs - min_mz) / bin_width).astype(int)

    mz_cals = indexed_mzcals[bin_indices]

    calibrated_mzs = mzs + mz_cals

    return calibrated_mzs

def calibrate_mzs(mzs_errors, exp):
    filtered = MSExperiment()

    bin_width = 1
    min_mz = 0
    max_mz = 2500
    bin_mzs = np.arange(min_mz, max_mz+bin_width, bin_width)
    theoretical_mzs, ppm_errors = mzs_errors[0,:], mzs_errors[1,:]
    theoretical_mzs = np.array(theoretical_mzs)
    ppm_errors = np.array(ppm_errors)

    indexed_ppmcals = np.interp(bin_mzs, theoretical_mzs, ppm_errors)
    mz_cals = (bin_mzs/1e6)*indexed_ppmcals

    for n,s in enumerate(exp):
        peaks = np.transpose(np.vstack(s.get_peaks()))
        mzs = peaks[:,0]
        intensities = peaks[:,1]

        new_mzs = calibrate_spectrum_mzs(mzs, mz_cals, min_mz, bin_width)
        s.set_peaks([new_mzs, intensities])
        if s.getMSLevel() == 2:
            p = s.getPrecursors()[0]
            p = fix_parent_mz(exp, n, p)
            if p is not None:
                s.setPrecursors([p])
                filtered.addSpectrum(s)
        else:
            filtered.addSpectrum(s)
    return filtered


equisplash_calibrants = [[796.65387, 8.5],
                         [753.54169,7.4],
                         [828.5627,7.4],
                         [811.6199, 9.4],
                         [709.5517,9.6],
                         [740.54669,7.85],
                         [1239.84, 16.1]]

toppath = "Z:\\Group Share\\Annika\\CS Collab Data\\Collaborator Samples (CS)\\Rep2\\mzML\\"
os.chdir(toppath)

mzML_files = match_files(toppath, "*.mzML", exclude="fixed")

for file in mzML_files:
    print(file)
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    filtered = MSExperiment()

    ind_file = IndexedFile()

    for i, s in enumerate(exp):
        if s.getMSLevel() == 1:
            peaks = np.transpose(np.vstack(s.get_peaks()))
            rt = s.getRT()/60
            ind_scan = IndexedScan(peaks, rt, i+1, 0,1)
            ind_file.indexed_scans.append(ind_scan)
    print("Done creating indexed file.")
    mz_errors = get_mzs_ppmerrors(ind_file, equisplash_calibrants)
    filtered = calibrate_mzs(mz_errors, exp)
    MzMLFile().store(file[:-5] + "_fixed.mzML", filtered)

