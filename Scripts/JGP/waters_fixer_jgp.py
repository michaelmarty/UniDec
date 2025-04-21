from pyopenms import *
import os
import fnmatch
from copy import deepcopy
import unidec.tools as ud
import xic_xtractor
import matplotlib.pyplot as plt
import time
import matplotlib as mpl

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

def calibrate_mzs(exp, deltas_ppm_experimental_mzs):
    timestart = time.perf_counter()
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
            p = fix_parent_mz(exp, n, p)
            if p is not None:
                s.setPrecursors([p])
                filtered.addSpectrum(s)
        else:
            filtered.addSpectrum(s)
    print("Time for Mass calibration: " + str(time.perf_counter() - timestart) + " seconds")
    return filtered

def calibrate_rts(exp, delta_rts_expected_rts):
    timestart = time.perf_counter()
    cal = MSExperiment()
    #first check for duplicate expected RTs in the delta_rts_expected_rts list
    #if there are duplicates, average the deltas for those RTs
    unique_rts = []

    for i in range(len(delta_rts_expected_rts)):
        if delta_rts_expected_rts[i][1] not in unique_rts:
            unique_rts.append(delta_rts_expected_rts[i][1])

    unique_deltas = np.zeros(len(unique_rts))
    for i in range(len(unique_rts)):
        count = 0
        for j in range(len(delta_rts_expected_rts)):
            if delta_rts_expected_rts[j][1] == unique_rts[i]:
                unique_deltas[i] += delta_rts_expected_rts[j][0]
                count += 1
        unique_deltas[i] /= count

    sorterd_indices = np.argsort(unique_rts)
    unique_rts = np.array(unique_rts)[sorterd_indices]
    unique_deltas = np.array(unique_deltas)[sorterd_indices]
    #now smooth the deltas
    for i in range(1, len(unique_deltas)-1):
        unique_deltas[i] = (unique_deltas[i-1] + unique_deltas[i] + unique_deltas[i+1])/3
    previous_rt = 0
    for n, s in enumerate(exp):
        s.setName(str(n))
        rt = s.getRT()/60
        new_rt = np.interp(rt, unique_rts, unique_deltas + unique_rts, left=0,right=rt + unique_deltas[len(unique_deltas)-1])
        s.setRT(new_rt * 60)
        if new_rt < previous_rt:
            print("RTs not in order, " + str(new_rt) + " is less than " + str(previous_rt))
        previous_rt = new_rt
        cal.addSpectrum(s)
    print("Time for RT calibration: " + str(time.perf_counter() - timestart) + " seconds")
    return cal


# Set backend to Agg
mpl.use('WxAgg')

toppath = "Z:\\Group Share\\Annika\\Whole Membrane\\Controls - 20250311"
os.chdir(toppath)

mzML_files = match_files(toppath, "*.mzML", exclude="fixed")

#now go through each file, find the equisplash-spiked molecules and
#calibrate based on the mz-errors of those.

equisplash_calibrants = [[796.65387, 6.8],
                         [753.64169,7.4],
                         [828.5627,6.4],
                         [811.6199, 7.42],
                         [709.5517,8.9],
                         [740.54669,6.7],
                         [1239.84, 14.1]]

rt_tol = 2
max_scans_to_check = 5
intensity_threshold = 500
mass_tol = 0.05

for file in mzML_files:
    starttime = time.perf_counter()
    print(file)
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    print("Creating indexed file")
    ind_file = xic_xtractor.IndexedFile(exp)
    print("Performing calibration")
    deltas = []
    deltas_ppm_expected_mzs = [[0,0]]
    delta_rts_expected_rts = [[0,0]]
    for i in range(len(equisplash_calibrants)):
        xic_start = time.perf_counter()
        xic = ind_file.extract_xic(equisplash_calibrants[i][0], mass_tol, equisplash_calibrants[i][1]-rt_tol, equisplash_calibrants[i][1]+rt_tol)

        intensity_column = xic[:,2]
        ordered_indices = np.argsort(intensity_column)[::-1]

        if xic[ordered_indices[0], 2] < intensity_threshold:
            print("Calibrant ", equisplash_calibrants[i][0], " not found in file ", file)
            continue

        rts = []
        delta_rts_expected_rts.append([equisplash_calibrants[i][1]  - xic[ordered_indices[0], 1]/60, equisplash_calibrants[i][1]])
        current_ind = 0
        while len(rts) < max_scans_to_check:
            rts.append(xic[ordered_indices[current_ind], 1])
            current_ind += 1

        for rt in rts:
            ind_scan = ind_file.get_indexed_spectrum_atRt(rt)
            deltas.append(equisplash_calibrants[i][0] - ind_scan.get_exp_mz(equisplash_calibrants[i][0], mass_tol))
            deltas_ppm_expected_mzs.append([((equisplash_calibrants[i][0] - ind_scan.get_exp_mz(equisplash_calibrants[i][0], mass_tol))/equisplash_calibrants[i][0]*1e6), equisplash_calibrants[i][0]])

    # cal = calibrate_rts(exp, delta_rts_expected_rts)
    cal = calibrate_mzs(exp, deltas_ppm_expected_mzs)
    del ind_file

    MzMLFile().store(file[:-5] + "_fixed.mzML", cal)
    print("Time taken: ", time.perf_counter() - starttime, " seconds")

