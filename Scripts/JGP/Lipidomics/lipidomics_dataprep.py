import unidec.IsoDec.datatools as dt
import os
import numpy as np
import fnmatch
from pyopenms import *
import Scripts.JGP.xic_xtractor as xic_xtractor
from copy import deepcopy
import unidec.tools as ud

equisplash_calibrants = [[796.65387, 6.5],
                         [753.54169,5.4],
                         [828.5627,5.4],
                         [811.6199, 7.4],
                         [709.5517,7.6],
                         [740.54669,5.85],
                         [1239.84, 14.1]]

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

def fast_remove_noise(cdata, factor):
    intensities = cdata[:,1]
    median = np.median(intensities)
    threhsold = median * factor
    cdata = [i for i in cdata if i[1] > threhsold]
    return np.vstack(cdata)

class LipidomicsPrepConfig:
    def __init__(self, **kwargs):
        self.remove_noise = True

        self.fix_parent_mz = True
        self.precursor_puritythresh = 0

        self.calibrate_mz = True
        self.calibration_mz_binwidth = 10
        self.xic_rtwidth = 5 #minutes
        self.xic_mzwidth = 20 #ppm
        self.min_peak_intensity = 1000
        self.scans_to_check = 5

def remove_noise_fullfile(exp):
    print("Removing noise...")
    noise_filtered = MSExperiment()
    for i,s in enumerate(exp):
        if i%100 == 0:
            print("Processing scan", i, "of", len(list(exp)), "scans...")
        cdata = np.transpose(np.array(s.get_peaks()))
        cdata_filtered = fast_remove_noise(cdata, 1.5)
        s.set_peaks([cdata_filtered[:,0], cdata_filtered[:,1]])
        noise_filtered.addSpectrum(s)

    return noise_filtered

def get_deltas_expecteds(indexed_file, prep_config, filename=None):
    deltas_ppm_expected_mzs = [[0, 0]]
    delta_rts_expected_rts = [[0, 0]]
    for i in range(len(equisplash_calibrants)):
        cal_mz = equisplash_calibrants[i][0]
        cal_rt = equisplash_calibrants[i][1]
        mass_tol = (cal_mz / 1000000) * prep_config.xic_mzwidth
        xic = indexed_file.extract_xic(cal_mz, mass_tol, cal_rt-prep_config.xic_rtwidth, cal_rt+prep_config.xic_rtwidth)

        intensity_column = xic[:,2]
        ordered_indices = np.argsort(intensity_column)[::-1]

        if xic[ordered_indices[0], 2] < prep_config.min_peak_intensity:
            if filename is not None:
                print("Calibrant ", equisplash_calibrants[i][0], " not found in file ", filename)
            continue

        rts = []
        delta_rts_expected_rts.append([cal_rt - xic[ordered_indices[0], 1]/60, cal_rt])
        current_ind = 0
        while len(rts) < prep_config.scans_to_check:
            rts.append(xic[ordered_indices[current_ind], 1])
            current_ind += 1

        for rt in rts:
            ind_scan = indexed_file.get_indexed_spectrum_atRt(rt)
            exp_mz = ind_scan.get_exp_mz(equisplash_calibrants[i][0], mass_tol)
            delta_ppm = ((cal_mz - exp_mz)/cal_mz)*1e6
            deltas_ppm_expected_mzs.append([delta_ppm, cal_mz])

    return deltas_ppm_expected_mzs, delta_rts_expected_rts

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
    chopdat = ud.datachop(data, pmass - window, pmass + window)
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

def calibrate_mzs(exp, deltas_ppm_experimental_mzs, prep_config):
    print("Calibrating m/zs...")
    filtered = MSExperiment()
    #need to both adjust the mz values in each spectrum and the parent mz values
    #first check for duplicate expected mzs in the deltas_ppm_experimental_mzs list
    #if there are duplicates, average the deltas for those mzs
    unique_mzs = []
    for i in range(len(deltas_ppm_experimental_mzs)):
        if deltas_ppm_experimental_mzs[i][1] not in unique_mzs:
            unique_mzs.append(deltas_ppm_experimental_mzs[i][1])

    unique_deltas = np.zeros(len(unique_mzs))
    for i in range(len(unique_mzs)):
        count = 0
        for j in range(len(deltas_ppm_experimental_mzs)):
            if deltas_ppm_experimental_mzs[j][1] == unique_mzs[i]:
                unique_deltas[i] += deltas_ppm_experimental_mzs[j][0]
                count += 1
        unique_deltas[i] /= count

    sorterd_indices = np.argsort(unique_mzs)
    unique_mzs = np.array(unique_mzs)[sorterd_indices]
    unique_deltas = np.array(unique_deltas)[sorterd_indices]

    min_mz = 100
    max_mz = 2500
    width = 25

    bins = (int)(max_mz - min_mz / width)
    indexed_cals = np.zeros(bins)

    current_mz = min_mz
    current_bin = 0
    while current_mz < max_mz:
        cal_ppm = np.interp(current_mz, unique_mzs, unique_deltas, left=0, right = unique_deltas[len(unique_deltas)-1])
        mass_cal = cal_ppm * current_mz / 1e6
        indexed_cals[current_bin] = mass_cal
        current_mz += width
        current_bin += 1


    for n, s in enumerate(exp):
        mz, intensity = s.get_peaks()
        new_mz = np.zeros(len(mz))
        for i in range(len(mz)):
            bin = (int)((mz[i] - min_mz) / width)
            new_mz[i] = mz[i] + indexed_cals[bin]
        s.set_peaks([new_mz, intensity])
        if s.getMSLevel() == 2:
            p = s.getPrecursors()[0]
            p = fix_parent_mz(exp, n, p, purity_thresh=prep_config.precursor_purity_thresh)
            if p is not None:
                s.setPrecursors([p])
                filtered.addSpectrum(s)
        else:
            filtered.addSpectrum(s)
    return filtered

def prepare_file(file, prep_config):
    print("Preparing file", file)
    exp = MSExperiment()
    MzMLFile().load(file, exp)

    if prep_config.remove_noise:
        exp = remove_noise_fullfile(exp)
    ind_file = xic_xtractor.IndexedFile(exp)

    if prep_config.calibrate_mz:
        delta_ppm_expected, delta_rts_expected = get_deltas_expecteds(ind_file, prep_config, filename=file)
        cal = calibrate_mzs(exp, delta_ppm_expected, prep_config)
    else:
        cal = exp
    MzMLFile().store(file[:-5] + "_fixed.mzML", cal)
    print("Finished preparing file ", file)
    return

if __name__ == "__main__":
    config = LipidomicsPrepConfig()

    file_directory = r"Z:\Group Share\JGP\Ecoli + Equisplash MS1\mzML"
    os.chdir(file_directory)

    mzml_files = match_files(file_directory, "*.mzML", exclude="fixed.mzML")

    for file in mzml_files:
        prepare_file(file, config)