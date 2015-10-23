import os
import platform
import sys
import math
import subprocess
from bisect import bisect_left
from ctypes import *
from copy import deepcopy
import zipfile
import numpy as np
import scipy.ndimage.filters as filt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy import stats
from scipy import signal
from scipy import special
import matplotlib.cm as cm

import mzMLimporter

is_64bits = sys.maxsize > 2 ** 32

"""
Loads initial dll called libmypfunc. Will speed up convolutions and a few functions.

If this isn't present, it will print a warning but use the pure python code later on.
"""

# TODO: Make this import more robust.
dllname = "libmypfunc"

if platform.system() == "Windows":
    if not is_64bits:
        dllname += "32"
    dllname += ".dll"
elif platform.system() == "Darwin":
    dllname += ".dylib"
else:
    dllname += ".so"
testpath = dllname
if os.path.isfile(testpath):
    dllpath = testpath
else:
    # print testpath
    pathtofile = os.path.dirname(os.path.abspath(__file__))
    testpath = os.path.join(pathtofile, dllname)
    if os.path.isfile(testpath):
        dllpath = testpath
    else:
        # print testpath
        testpath = os.path.join(os.path.dirname(pathtofile), dllname)
        if os.path.isfile(testpath):
            dllpath = testpath
        else:
            # print testpath
            testpath = os.path.join(os.path.join(os.path.dirname(pathtofile), "unidec_bin"), dllname)
            if os.path.isfile(testpath):
                dllpath = testpath
            else:
                # print testpath
                print "Unable to find file"

try:
    libs = cdll.LoadLibrary(dllpath)
except WindowsError:
    print dllpath
    print "Failed to load libmypfunc, convolutions in nonlinear mode might be slow"


# ..........................
#
# Utility Functions
#
# ..................................


def isempty(thing):
    """
    Checks wether the thing, a list or array, is empty
    :param thing: Object to check
    :return: Boolean, True if not empty
    """
    try:
        if np.asarray(thing).size == 0 or thing is None:
            out = True
        else:
            out = False
    except TypeError or ValueError or AttributeError:
        print "Error testing emptiness"
        out = False
    return out


def string_to_value(s):
    """
    Try to convert string to float.
    :param s: String
    :return: Float if possible. Otherwise, empty string.
    """
    try:
        v = float(s)
        return v
    except (ValueError, TypeError):
        return ""


def simp_string_to_value(s):
    """
    Try to convert string to float. Raise False if not possible.
    :param s: String
    :return: Float if possible. Otherwise, False
    """
    try:
        v = float(s)
        return v
    except (ValueError, TypeError):
        return False


def string_to_int(s):
    """
    Try to convert string to integer.
    :param s: String
    :return: Integer if possible. Otherwise, False
    """
    try:
        v = int(s)
        return v
    except (ValueError, TypeError):
        return False


def safedivide(a, b):
    """
    Divide numpy arrays b from a. Define a/0=0 to avoid zero division error.
    :param a: Numerator
    :param b: Denominator
    :return: a/b
    """
    c = deepcopy(b)
    c[b != 0] = a[b != 0] / b[b != 0]
    return c


def weighted_std(values, weights):
    """
    Calculate weighted standard deviation.
    :param values: Values
    :param weights: Weights
    :return: Weighted standard deviation.
    """
    average = np.average(values, weights=weights)
    variance = np.average((np.array(values) - average) ** 2, weights=weights)  # Fast and numerically precise
    return np.sqrt(variance)


def interp_pos(array, target):
    """
    On a given array, interpolate the position of the target
    :param array: Array
    :param target: Value
    :return: Interpolated index of value in array.
    """
    indexes = np.arange(0, len(array))
    f = interp1d(array, indexes)
    out = f(target)
    return out


def nearestunsorted(array, target):
    """
    In an unsorted array, find the position of the element closest to the target.

    For a sorted array, the method nearest will be faster.
    :param array: Array
    :param target: Value
    :return: np.argmin(np.abs(array - target))
    """
    return np.argmin(np.abs(array - target))


def nearest(array, target):
    """
    In an sorted array, quickly find the position of the element closest to the target.
    :param array: Array
    :param target: Value
    :return: np.argmin(np.abs(array - target))
    """
    i = bisect_left(array, target)
    if i <= 0:
        return 0
    elif i >= len(array) - 1:
        return len(array) - 1
    if np.abs(array[i] - target) > np.abs(array[i - 1] - target):
        i -= 1
    return i


def get_z_offset(mass, charge):
    """
    For a given mass and charge combination, calculate the charge offset parameter.
    :param mass: Mass in Da
    :param charge: Charge State
    :return: Difference between charge and predicted native charge of the mass
    """
    nativez = predict_charge(mass)
    return charge - nativez


def predict_charge(mass):
    """
    Give predicted native charge state of species of defined mass
    :param mass: Mass in Da
    :return: Float, predicted average native charge state
    """
    m = 0.0467
    s = 0.533
    nativez = m * (mass ** s)
    return nativez


def get_tvalue(dof, ci=0.99):
    """
    Get t-value for a given degrees of freedom and confidence interval
    :param dof: Degrees of freedom
    :param ci: Confidence interval level (default is 0.99, 99% confidence)
    :return: T value
    """
    return stats.t.isf((1 - ci) / 2., dof)


def get_zvalue(ci=0.99):
    """
    Get z-value for area under a Gaussian for a given confidence level
    :param ci: Confidence interval level (default is 0.99, 99% confidence)
    :return: Z-value
    """
    return stats.norm.isf((1 - ci) / 2.)


# ............................
#
# File manipulation
#
# .........................................


def header_test(path):
    """
    A quick test of a file to remove any problems from the header.

    If scans through each line of an input file specificed by path and count the number of lines that cannot be
    completely converted to floats. It breaks the loop on the first hit that is entirely floats, so any header lines
    blow that will not be counted.
    :param path: Path of file to be scanned
    :return: Integer length of the number of header lines
    """
    header = 0
    try:
        with open(path, "r") as f:
            for line in f:
                for sin in line.split():
                    try:
                        float(sin)
                    except ValueError:
                        header += 1
                        break
        if header > 0:
            print "Header Length:", header
    except ImportError or WindowsError or AttributeError or IOError:
        print "Failed header test"
        header = 0
    return header


def load_mz_file(path, config):
    """
    Loads a text or mzml file
    :param path: File path to load
    :param config: UniDecConfig object
    :return: Data array
    """
    if config.extension.lower() == ".txt":
        data = np.loadtxt(path, skiprows=header_test(path))
    elif config.extension.lower() == ".mzml":
        data = mzMLimporter.mzMLimporter(path).get_data()
        txtname = path[:-5] + ".txt"
        np.savetxt(txtname, data)
        print "Saved to:", txtname
    else:
        try:
            data = np.loadtxt(path, skiprows=header_test(path))
        except ImportError or IOError:
            print"Failed to open:", path
            data = None
    return data


def zipdir(path, zip_handle):
    """
    Zips all the files in the path into the zip_handle
    :param path: Directory path
    :param zip_handle: Handle of the zip file that is being created
    :return: None
    """
    files = os.listdir(path)
    for f in files:
        if os.path.isfile(f):
            zip_handle.write(f)


def zip_folder(save_path):
    """
    Zips a directory specified by save_path into a zip file for saving.
    :param save_path: Path to save to zip
    :return: None
    """
    directory = os.getcwd()
    print "Zipping directory:", directory
    zipf = zipfile.ZipFile(save_path, 'w')
    zipdir(directory, zipf)
    zipf.close()
    print "File saved to: " + str(save_path)


def dataexport(datatop, fname):
    np.savetxt(fname, datatop, fmt='%f')
    pass


def mergedata(data1, data2):
    """
    Uses a 1D intepolation to merge the intensity of data2 into the axis of data1
    :param data1: Reference data (N x 2 array)
    :param data2: Data to merge (M x 2 array)
    :return: data2 interpolated into the same axis as data1 (N x 2)
    """
    f = interp1d(data2[:, 0], data2[:, 1], bounds_error=False, fill_value=0)
    intx = data1[:, 0]
    inty = f(intx)
    newdat = np.column_stack((intx, inty))
    return newdat


def mergedata2d(x1, y1, x2, y2, z2):
    """
    Uses a 2D interpolation to merge a 2D grid of data from one set of axes to another.

    Requires rectangular grids of all inputs such as might be created by numpy.meshgrid.
    :param x1: Reference x-value grid (N x M)
    :param y1: Reference y-value grid (N x M)
    :param x2: Data x-value grid (N2 x M2)
    :param y2: Data y-value grid (N2 x M2)
    :param z2: Data intensity grdi (N2 x M2)
    :return: z2 interpolated onto x1 y1 2D grid (N x M)
    """
    oldx = np.ravel(x2)
    oldy = np.ravel(y2)
    newx = np.ravel(x1)
    newy = np.ravel(y1)
    zout = griddata(np.transpose([oldx, oldy]), np.ravel(z2), (newx, newy), method='linear', fill_value=0)
    return zout


# ........................................
#
# Data Processing
#
# ..........................................................


def auto_peak_width(datatop):
    maxpos = np.argmax(datatop[:, 1])
    maxval = datatop[maxpos, 0]

    ac, cpeaks = autocorr(datatop)
    sig = cpeaks[0, 0] / 4.

    boo1 = datatop[:, 0] < maxval + sig
    boo2 = datatop[:, 0] > maxval - sig
    boo3 = np.all([boo1, boo2], axis=0)

    isodat = datatop[boo3]
    fits = np.array([isolated_peak_fit(isodat[:, 0], isodat[:, 1], i) for i in xrange(0, 3)])

    errors = [np.sum(np.array(isodat[:, 1] - f) ** 2.) for f in fits[:, 1]]
    psfun = np.argmin(errors)
    fit, fitdat = fits[psfun]
    fwhm = fit[0, 0]

    return fwhm, psfun, maxval


def auto_noise_level(datatop, buffer=10):
    ndat = np.ravel([datatop[:buffer, 1], datatop[-buffer:, 1]])
    std = np.std(ndat)
    mean = np.mean(ndat)
    return mean + 5 * std


def average_bin_size(datatop):
    diffs = np.diff(datatop[:, 0])
    return np.average(diffs)


def datachop(datatop, newmin, newmax):
    """
    Chops the range of the data. The [:,0] column of the data is the column that will be indexed.
    :param datatop: Data array
    :param newmin: Minimum value of chopped data
    :param newmax: Maximum value of chopped data
    :return: New data limited between the two bounds
    """
    boo1 = np.logical_and(datatop[:, 0] < newmax, datatop[:, 0] > newmin)
    return datatop[boo1]


def datasimpsub(datatop, buff):
    """
    Simple subtraction of a linear background.

    The first n and last n data points are averaged. A line is drawn between the two averages.
    This line is then subtracted from the [:,1] column of the data.
    :param datatop: Data array (N x 2) (m/z, intensity)
    :param buff: Integer, n.
    :return: Subtracted data
    """
    length = len(datatop)
    frontpart = np.mean(datatop[:buff, 1])
    backpart = np.mean(datatop[length - buff - 1:length - 1, 1])
    background = frontpart + (backpart - frontpart) / length * np.arange(length)
    datatop[:, 1] = datatop[:, 1] - background
    return datatop


def datacompsub(datatop, buff):
    """
    Complex background subtraction.

    Taken from Massign Paper

    First creates an array that matches the data but has the minimum value within a window of +/- buff.
    Then, smooths the minimum array with a Gaussian filter of width buff * 2 to form the background array.
    Finally, subtracts the background array from the data intensities.

    :param datatop: Data array
    :param buff: Width parameter
    :return: Subtracted data
    """
    length = len(datatop)
    mins = range(0, length)
    indexes = range(0, length)
    for i in indexes:
        mins[i] = np.amin(datatop[max([0, i - abs(buff)]):min([i + abs(buff), length]), 1])
    background = filt.gaussian_filter(mins, abs(buff) * 2)
    datatop[:, 1] = datatop[:, 1] - background
    return datatop


def intensitythresh(datatop, thresh):
    """
    Sets an intensity threshold. Everything below the threshold is set to 0.
    :param datatop: Data array
    :param thresh: Threshold value
    :return: Thresholded array.
    """
    belowint = datatop[:, 1] < thresh
    datatop[belowint, 1] = 0
    return datatop


def gsmooth(datatop, sig):
    """
    Smooths the data with a Gaussian filter of width sig.
    :param datatop: Data array
    :param sig: Width of Gaussian Array
    :return: Smoothed Data
    """
    datatop[:, 1] = filt.gaussian_filter(datatop[:, 1], sig)
    return datatop


def nonlinear_axis(start, end, res):
    """
    Creates a nonlinear axis with the m/z values spaced with a defined and constant resolution.
    :param start: Minimum m/z value
    :param end: Maximum m/z value
    :param res: Resolution of the axis ( m / delta m)
    :return: One dimensional array of the nonlinear axis.
    """
    axis = []
    i = start
    axis.append(i)
    i += i / float(res)
    while i < end:
        axis.append(i)
        i += i / float(res)
    return np.array(axis)


def linear_interpolation(x1, x2, x):
    """
    :param x1:
    :param x2:
    :param x:
    :return: float(x - x1) / float(x2 - x1)
    """
    return float(x - x1) / float(x2 - x1)


def lintegrate(datatop, intx):
    """
    Linearize x-axis by integration.

    Each intensity value in the old data gets proportionally added into the new x-axis.

    The total sum of the intensity values should be constant.
    :param datatop: Data array
    :param intx: New x-axis for data
    :return: Integration of intensity from original data onto the new x-axis.
        Same shape as the old data but new length.
    """
    length = len(datatop)
    inty = np.zeros_like(intx)
    for i in range(0, length):
        if intx[0] < datatop[i, 0] < intx[len(intx) - 1]:
            index = nearest(intx, datatop[i, 0])
            # inty[index]+=datatop[i,1]
            if intx[index] == datatop[i, 0]:
                inty[index] += datatop[i, 1]
            if intx[index] < datatop[i, 0] and index < length - 1:
                index2 = index + 1
                interpos = linear_interpolation(intx[index], intx[index2], datatop[i, 0])
                inty[index] += (1 - interpos) * datatop[i, 1]
                inty[index2] += interpos * datatop[i, 1]
            if intx[index] > datatop[i, 0] and index > 0:
                index2 = index - 1
                interpos = linear_interpolation(intx[index], intx[index2], datatop[i, 0])
                inty[index] += (1 - interpos) * datatop[i, 1]
                inty[index2] += interpos * datatop[i, 1]
    newdat = np.column_stack((intx, inty))
    return newdat


def linterpolate(datatop, intx):
    """
    Linearize x-axis by interpolation.

    The new x-axis is interpolated on the old data and the corresponding intensities and picked out.

    :param datatop: Data array
    :param intx: New x-axis
    :return: Interpolation of intensity from original data onto the new x-axis.
        Same shape as the old data but new length.
    """
    f = interp1d(datatop[:, 0], datatop[:, 1])
    inty = f(intx)
    newdat = np.column_stack((intx, inty))
    return newdat


def linearize(datatop, binsize, linflag):
    """
    Linearize data by defined x-axis bin size by either interpolation or integration
    :param datatop: Data array (N x 2)
    :param binsize: Spacing between first two m/z points
    :param linflag: Integer flag indicating the type of linearization.

    0=Integration with Linear Axis (constant m/z)
    1=Integration with Nonlinear Axis (constant resolution: mz / delta mz)

    3=Interpolation with Linear Axis
    4=Interpolation with Nonlinear Axis

    :return: Linearized data
    """
    length = len(datatop)
    firstpoint = math.ceil(datatop[0, 0] / binsize) * binsize
    lastpoint = math.floor(datatop[length - 1, 0] / binsize) * binsize
    if linflag == 0 or linflag == 3:
        intx = np.arange(firstpoint, lastpoint, binsize)
    else:
        intx = nonlinear_axis(firstpoint, lastpoint, firstpoint / binsize)

    if linflag < 2:
        newdat = lintegrate(datatop, intx)
        print "Integrating"
    else:
        newdat = linterpolate(datatop, intx)
        print "Interpolating"
    return newdat


def nonlinearize(datatop, num_compressed):
    """
    Compress the data in a simple and nonlinear way.
    :param datatop: Data array (N x 2)
    :param num_compressed:
    :return:
    """
    if num_compressed == 0:
        return datatop
    else:
        num_compressed = int(num_compressed)
        return np.array([np.mean(datatop[index:index + num_compressed], axis=0) for index in
                         xrange(0, len(datatop), num_compressed)])


def removeduplicates(datatop):
    """
    Cleans up data that may have the same x value multiple times.
    Each nonunique x gets merged into a single unique x with the sum of all of the nonunique intensities.
    :param datatop: Data array (N x 2)
    :return: Data with unique x values
    """
    testunique = np.unique(datatop[:, 0])
    if len(testunique) != len(datatop):
        print "Removing Duplicates"
        num, start = np.histogram(datatop[:, 0], bins=testunique)
        means = []
        xvals = []
        index = 0
        for i in range(0, len(testunique) - 1):
            xvals.append(np.mean(datatop[index:index + num[i], 0]))
            means.append(np.sum(datatop[index:index + num[i], 1]))
            index = index + num[i]
        datatop = np.column_stack((xvals, means))
    return datatop


def normalize(datatop):
    """
    Normalize the data so that the maximum intensity is 1.

    maxval = np.amax(datatop[:, 1])
    datatop[:, 1] = datatop[:, 1] / maxval

    :param datatop: Data array (N x 2)
    :return: Normalized data array (N x 2)
    """
    maxval = np.amax(datatop[:, 1])
    datatop[:, 1] = datatop[:, 1] / maxval
    return datatop


def detectoreff(datatop, va):
    """
    Incorporate a model of TOF detector efficiency used by CHAMP to correct the intensity values.
    :param datatop: Data array (N x 2)
    :param va: TOF acceleration voltage in kV (I think...)
    :return: Detector corrected data (N x 2)
    """
    eff = (1 - np.exp(-1620 * (va / datatop[:, 0]) ** 1.75))
    datatop[:, 1] = datatop[:, 1] / eff
    return datatop


def dataprep(datatop, config):
    """
    Main function to process 1D MS data. The order is:

    Crop Data (datacrop)
    Correct for detector efficiency (detectoreff)
    Smooth (gsmooth)
    Remove duplicates (removeduplicates)
    Linearize (linearize or nonlinearize)
    Baseline subtraction (datasimpsub or datacompsub or simply data[:,1]-np.amin(data[:,1])
    Normalization (normalize)
    Intensity Threshold with a threshold of 0 (intthresh)

    :param datatop: Raw data array (N x 2)
    :param config: UniDecConfig object
    :return: Processed data
    """
    newmin = config.minmz
    newmax = config.maxmz
    buff = config.subbuff
    smooth = config.smooth
    binsize = config.mzbins
    # thresh = config.intthresh
    subtype = config.subtype
    va = config.detectoreffva
    linflag = config.linflag
    # Crop Data
    data2 = datachop(datatop, newmin, newmax)

    # correct for detector efficiency
    if va != 0:
        # data2=detectoreff(data2,9.1)
        data2 = detectoreff(data2, va)

    # Smooth Data
    if smooth > 0:
        data2 = gsmooth(data2, smooth)

    # Remove Duplicate Data Points
    data2 = removeduplicates(data2)
    # Linearize Data
    if binsize > 0:
        if linflag != 2:
            data2 = linearize(data2, binsize, linflag)
        else:
            data2 = nonlinearize(data2, binsize)

    # Baseline Subtraction
    buff = abs(buff)
    if subtype == 1 and buff != 0:
        data2 = datasimpsub(data2, buff)
    elif subtype == 2 and buff != 0:
        data2 = datacompsub(data2, buff)
    elif subtype == 0 and buff != 0:
        data2[:, 1] = data2[:, 1] - np.amin(data2[:, 1])

    # Normalization
    data2 = normalize(data2)

    # Intensity Threshold
    data2 = intensitythresh(data2, 0)  # thresh
    # data2=data2[data2[:,1]>0]
    return data2


# ......................................................
#
# UniDec Functions
#
# ............................................................................


def unidec_call(exepath, configfile, silent=False, **kwargs):
    """
    Run the UniDec binary specified by exepath with the configuration file specified by configfile.
    If silent is False (default), the output from exepath will be printed to the standard out.
    If silent is True, the output is suppressed.

    The binary is run as a commandline subprocess with a call of (for example):

    UniDec.exe conf.dat

    :param exepath: Path to UniDec or UniDecIM binary
    :param configfile: Path to configuration file
    :param silent: Whether to print the output of exepath to the standard out
    :param kwargs:
    :return: Standard error of exepath execution
    """
    if "kill" in kwargs:
        killnum = kwargs["kill"]
        # print "Killing Peak:",killnum
        out = subprocess.call([exepath, configfile, "-kill", str(killnum)], stdout=subprocess.PIPE)
    else:
        if silent:
            out = subprocess.call([exepath, configfile], stdout=subprocess.PIPE)
        else:
            out = subprocess.call([exepath, configfile])
    return out


def peakdetect(data, config=None, window=10, threshold=0):
    """
    Simple peak detection algorithm.

    Detects a peak if a given data point is a local maximum within plus or minus config.peakwindow.
    Peaks must also be above a threshold of config.peakthresh * max_data_intensity.

    The mass and intensity of peaks meeting these criteria are output as a P x 2 array.

    :param data: Mass data array (N x 2) (mass intensity)
    :param config: UniDecConfig object
    :return: Array of peaks positions and intensities (P x 2) (mass intensity)
    """
    if config is not None:
        window = config.peakwindow / config.massbins
        threshold = config.peakthresh
    peaks = []
    length = len(data)
    maxval = np.amax(data[:, 1])
    for i in range(1, length):
        if data[i, 1] > maxval * threshold:
            start = i - window
            end = i + window
            if start < 0:
                start = 0
            if end > length:
                end = length
            testmax = np.amax(data[start:end + 1, 1])
            if data[i, 1] == testmax and data[i, 1] != data[i - 1, 1]:
                peaks.append([data[i, 0], data[i, 1]])
    return np.array(peaks)


def mergepeaks(peaks1, peaks2, window):
    """
    Merge a test list of peaks to a reference list.

    Each peak in peaks2 is matched to the closest peak in peak1.
    If it within some defined threshold, the new peak replaces that closes peak in the matched list.
    If no new peaks are matched to a reference peak, that peaks mass is unchanged but intensity is set to 0.

    :param peaks1: Reference peaks (P x 2)
    :param peaks2: Test peaks (R x 2)
    :param window: Tolerance window of the x values
    :return: Intensity of peaks2 matched to values of peaks1 (P x 2)
    """
    newpeaks = deepcopy(peaks1)
    newpeaks[:, 1] = newpeaks[:, 1] * 0
    for peak in peaks2:
        i = np.argmin(abs(peaks1[:, 0] - peak[0]))
        minval = np.amin(abs(peaks1[:, 0] - peak[0]))
        closest = np.sort(np.abs(peaks2[:, 0] - peaks1[i, 0]))[0]
        if minval < window and minval == closest:
            newpeaks[i] = peak
    return newpeaks


def make_peaks_mztab(mzgrid, pks, adductmass):
    """
    For each peak in pks, get the charge state distribution.

    The intensity at each charge state is determined from mzgrid and stored in a list as peak.mztab.
    An array of the results is also output for help speeding up make_peaks_mztab_spectrum

    :param mzgrid: 2D grid of m/z vs charge (N x Z)
    :param pks: Peaks object (length P)
    :param adductmass: Mass of electrospray adduct.
    :return: P x Z array
    """
    xvals = np.unique(mzgrid[:, 0])
    xmin = np.amin(xvals)
    xmax = np.amax(xvals)
    yvals = np.unique(mzgrid[:, 1])
    xlen = len(xvals)
    ylen = len(yvals)
    newgrid = np.reshape(mzgrid[:, 2], (xlen, ylen))
    plen = pks.plen
    ftab = [interp1d(xvals, newgrid[:, k]) for k in xrange(0, ylen)]
    mztab = [[makespecfun(i, k, pks.masses, adductmass, yvals, xvals, ftab, xmax, xmin) for k in xrange(0, ylen)] for i
             in xrange(0, plen)]
    for i in xrange(0, plen):
        pks.peaks[i].mztab = np.array(mztab[i])
    return np.array(mztab)


def makespecfun(i, k, peaks_masses, adductmass, charges, xvals, ftab, xmax, xmin):
    """
    Small helper function to speed things up. May be more better way to do this, but I found this to be fastest.
    :param i: peak mass index
    :param k: charge index
    :param peaks_masses: peak masses list
    :param adductmass: Adduct mass
    :param charges: charge list
    :param xvals: x-values from data
    :param ftab: List of interpolation functions
    :param xmax: Maximum x value
    :param xmin: Minimum x value
    :return:
    """
    intx = np.true_divide((peaks_masses[i] + adductmass * charges[k]), charges[k])
    if xmin < intx < xmax:
        f = ftab[k]
        inty = f(intx)
    else:
        inty = 0
    pos = nearest(xvals, intx)
    return np.array([intx, inty, pos])


def make_peaks_mztab_spectrum(mzgrid, pks, data2, mztab):
    """
    Used for plotting the dots in plot 4.

    Perform the same things as make_peaks_mztab, but instead of using the deconvolved intensities,
    it takes the intensity directly from the spectrum.
    :param mzgrid: m/z grid (N x Z)
    :param pks: Peaks object
    :param data2: Processed data 1D
    :param mztab: Prior mztab from make_peaks_mztab
    :return: mztab but with intensities replaced by value at the spectrum.
    """
    zvals = np.unique(mzgrid[:, 1])
    zlen = len(zvals)
    plen = pks.plen
    mztab2 = deepcopy(mztab)
    mztab2[:, :, 1] = [[data2[pks.peaks[i].mztab[k, 2], 1] for k in range(0, zlen)] for i in range(0, plen)]
    for i in xrange(0, plen):
        pks.peaks[i].mztab2 = np.array(mztab2[i])

    return mztab2


def makeconvspecies(processed_data, pks, config):
    """
    For data and peaks, create a simulated spectrum of each peak.
    Assumes p.mztab has been set for each p in pks.peaks by make_peaks_mztab
    Also makes pks.composite, which is a sum of all the simulated spectra.
    :param processed_data: Processed 2D MS data (data2, N x 2)
    :param pks: Peaks object (P Peak objects in pks.peaks)
    :param config: UniDecConfig object
    :return: Array of simulated spectra (P x N)
    """
    xvals = processed_data[:, 0]
    xlen = len(xvals)

    if config.linflag != 2:
        peakwidth = config.mzsig / config.mzbins
        kernel = conv_peak_shape_kernel(xvals, config.psfun, peakwidth)
        stickdat = [stickconv(p.mztab, kernel) for p in pks.peaks]
    else:
        try:
            stickdat = [cconvolve(xvals, p.mztab, config.mzsig, config.psfun) for p in pks.peaks]
        except WindowsError or TypeError or NameError or AttributeError:
            stickdat = [nonlinstickconv(xvals, p.mztab, config.mzsig, config.psfun) for p in pks.peaks]

    pks.composite = np.zeros(xlen)
    for i in xrange(0, pks.plen):
        pks.peaks[i].stickdat = stickdat[i]
        pks.composite += np.array(stickdat[i])
    pks.convolved = True
    return np.array(stickdat)


def nonlinstickconv(xvals, mztab, fwhm, psfun):
    """
    Python-based Nonlinear convolution. First makes a stick spectrum. Then, convolved with peak shape kernel.
    :param xvals: x-axis
    :param mztab: mztab from make_peaks_mztab
    :param fwhm: Full width half max
    :param psfun: Peak shape function integer
    :return: Convolved output
    """
    if psfun == 0:
        window = 5 * fwhm
    else:
        window = 15 * fwhm
    xlen = len(xvals)
    stick = np.zeros(xlen)
    stick[np.array(mztab[:, 2]).astype(np.int)] = mztab[:, 1]
    bool1 = [np.abs(xvals - xvals[i]) < window for i in xrange(0, xlen)]
    kernels = np.array([make_peak_shape(-xvals[bool1[i]], psfun, fwhm, -xvals[i]) for i in xrange(0, xlen)])
    output = np.array([np.sum(kernels[i] * stick[bool1[i]]) for i in xrange(0, xlen)])
    return output


def stickconv(mztab, kernel):
    """
    Make stick spectrum and then convolves with kernel.
    :param mztab: mztab from make_peaks_mztab
    :param kernel: peak shape kernel
    :return: Convolved output
    """
    xlen = len(kernel)
    temp = np.zeros(xlen)
    temp[np.array(mztab[:, 2]).astype(np.int)] = mztab[:, 1]
    return cconv(temp, kernel)


def make_alpha_cmap(rgb_tuple, alpha):
    """
    Make color map where RGB specified in tup
    :param rgb_tuple: Tuple of RGB vaalues [0,1]
    :param alpha: Maximum Alpha (Transparency) [0,1]
    :return: Color map dictionary with a constant color but varying transprency
    """
    cdict = {'red': ((0.0, rgb_tuple[0], rgb_tuple[0]),
                     (1.0, rgb_tuple[0], rgb_tuple[0])), 'green': ((0.0, rgb_tuple[1], rgb_tuple[1]),
                                                                   (1.0, rgb_tuple[1], rgb_tuple[1])),
             'blue': ((0.0, rgb_tuple[2], rgb_tuple[2]),
                      (1.0, rgb_tuple[2], rgb_tuple[2])), 'alpha': ((0.0, 0, 0),
                                                                    (1.0, alpha, alpha))}
    return cdict


def color_map_array(array, cmap, alpha):
    """
    For a specified array of values, map the intensity to a specified RGB color defined by cmap (output as topcm).
    For each color, create a color map where the color is constant but the transparency changes (output as cmarr).
    :param array: Values
    :param cmap: Color map
    :param alpha: Max Alpha value (transparency) [0,1]
    :return: cmarr, topcm (list of transparent color maps, list of RGB colors for array values)
    """
    rtab = array
    topcm = cm.ScalarMappable(cmap=cmap).to_rgba(rtab)[:, :3]
    cmarr = []
    for i in range(0, len(rtab)):
        cmarr.append(make_alpha_cmap(topcm[i], alpha))
    return cmarr, topcm


# ......................................................
#
# Peak Shape Function
#
# ........................................................................


def ndis_std(x, mid, sig, a=1, norm_area=False):
    """
    Normal Gaussian function normalized to the max of 1.
    :param x: x values
    :param mid: Mean of Gaussian
    :param sig: Standard Deviation
    :param a: Maximum amplitude (default is 1)
    :param norm_area: Boolean, Whether to normalize so that the area under the distribution is 1.
    :return: Gaussian distribution at x values
    """
    if norm_area:
        a *= 1 / (sig * np.sqrt(2 * np.pi))
    return a * np.exp(-(x - mid) * (x - mid) / (2.0 * sig * sig))


def ndis(x, mid, fwhm, **kwargs):
    """
    Gaussian function normalized to a max of 1.

    Note: x and mid are interchangable. At least one should be a single float. The other may be an array.
    :param x: x values
    :param mid: Mean
    :param fwhm: Full Width at Half Max (2.35482 * standard deviation)
    :param kwargs: Allows norm_area flag to be passed
    :return: Gaussian distribution at x values
    """
    sig = fwhm / 2.35482
    return ndis_std(x, mid, sig, **kwargs)


def ndis_fit(x, s, m, a, b):
    """
    Function for fitting normal distribution to peak.
    Adds a background to ndis.
    Prevents negative background, amplitude, and standard deviation.
    :param x: x value
    :param s: full width half max
    :param m: mean
    :param a: amplitude
    :param b: linear background
    :return: peak shape
    """
    if b < 0 or a < 0 or s < 0:
        return x * 0
    return ndis(x, m, s, a=a, norm_area=True) + b


def ldis(x, mid, fwhm, a=1, norm_area=False):
    """
    Lorentzian function normalized to a max of 1.
    Note: x and mid are interchangable. At least one should be a single float. The other may be an array.
    :param x: x values
    :param mid: Mean
    :param fwhm: Full Width at Half Max
    :param a: Amplitude (default is 1)
    :param norm_area: Boolean, Whether to normalize so that the area under the distribution is 1.
    :return: Lorentzian distribution at x values
    """
    if norm_area:
        a *= ((1 / np.pi) * (fwhm / 2.))
    else:
        a *= ((fwhm / 2.0) * (fwhm / 2.0))
    return a / ((x - mid) * (x - mid) + (fwhm / 2.0) * (fwhm / 2.0))


def ldis_fit(x, s, m, a, b):
    """
    Function for fitting Lorentzian distribution to peak.
    Adds a background to ldis.
    Prevents negative background, amplitude, and standard deviation.
    :param x: x value
    :param s: full width half max
    :param m: mean
    :param a: amplitude
    :param b: linear background
    :return: peak shape
    """
    if b < 0 or a < 0 or s < 0:
        return x * 0
    return ldis(x, m, s, a=a, norm_area=True) + b


def splitdis(x, mid, fwhm, a=1, norm_area=False):
    """
    Split Gaussain/Lorentzian function normalized to a max of 1.

    Gaussian < mid. Lorentzian > mid.

    :param mid: Mid point (point of peak intensity)
    :param x: x value or values
    :param fwhm: Full Width at Half Max
    :return: Split Gaussian/Lorentzian distribution at x value
    """
    sig2 = fwhm / (2 * np.sqrt(2 * np.log(2)))
    if norm_area:
        a1 = a * ((1 / np.pi) / (fwhm / 2.)) / 0.83723895067
        a2 = a * 2. / (fwhm * np.pi) / 0.83723895067
    else:
        a1 = a
        a2 = a
    try:
        if mid < x:
            return ldis(x, mid, fwhm, a=a1)
        else:
            return ndis_std(x, mid, sig2, a=a2)
    except ValueError:
        output = np.zeros(len(x))
        output[x > mid] = ldis(x[x > mid], mid, fwhm, a=a1)
        output[x <= mid] = ndis_std(x[x <= mid], mid, sig2, a=a2)
        return output


def splitdis_fit(x, s, m, a, b):
    """
    Function for fitting Split G/L distribution to peak.
    Adds a background to splitdis.
    Prevents negative background, amplitude, and standard deviation.
    :param x: x value
    :param s: full width half max
    :param m: mean
    :param a: amplitude
    :param b: linear background
    :return: peak shape
    """
    if b < 0 or a < 0 or s < 0:
        return x * 0
    return splitdis(x, m, s, a=a, norm_area=True) + b


def voigt(x, mu=0, sigma=1, gamma=1, amp=1, background=0):
    """\
    voigt profile

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))
    """
    if sigma == 0:
        return ldis(x, mu, gamma*2., amp) + background
    elif gamma == 0:
        return ndis_std(x, mu, sigma, amp) + background
    else:
        z = (x - mu + 1j * gamma) / (sigma * np.sqrt(2))
        w = special.wofz(z)
        v = w.real / (sigma * np.sqrt(2 * np.pi))
        v *= (amp / np.amax(v))
        v += background
    return v


# ...........................................................
#
# Peak Shape Tools
#
# ........................................................................


def cconv(a, b):
    """
    Circular convolution of two lists
    :param a: N-length array
    :param b: N-length array
    :return: Convolution
    """
    # return np.fft.ifft(np.fft.fft(a) * np.fft.fft(b)).real
    # return np.convolve(a, np.roll(b, (len(b)) / 2 - 1 + len(b) % 2), mode="same")
    return signal.fftconvolve(a, np.roll(b, (len(b)) / 2 - 1 + len(b) % 2), mode="same")


def single_cwt(a, width, wavelet_type="Ricker"):
    """
    Perform a continuous wavelet transform at a single defined width.
    :param a: 1D numpy array of data (length N)
    :param width: Width of transform defined with respect to the data points in a
    :param wavelet_type: Type of wavelet. Either "Ricker" (Mexican Hat) or "Morlet" (Gabor)
    :return: cwt, wavelet
    Continous wavelet transform and wavelet of choice, both as numpy arrays of length N.
    """
    if wavelet_type is "Morlet":
        wdat = signal.morlet(len(a), width)
    else:
        wdat = signal.ricker(len(a), width)
    return continuous_wavelet_transform(a, [width], wavelet_type=wavelet_type)[0], wdat


def continuous_wavelet_transform(a, widths, wavelet_type="Ricker"):
    """
    Perform a continuous wavelet transform at multiple widths.
    :param a: 1D numpy array of data (length N)
    :param widths: List of widths (length W) of transform defined with respect to the data points in a
    :param wavelet_type: Type of wavelet. Either "Ricker" (Mexican Hat) or "Morlet" (Gabor)
    :return: cwt_matrix (The continuous wavelet transform at the defined widths in a (W x N) array)
    """
    if wavelet_type is "Morlet":
        wavelet = signal.morlet
    else:
        wavelet = signal.ricker
    return signal.cwt(a, wavelet, widths)


def autocorr(datatop, config=None):
    corry = signal.fftconvolve(datatop[:, 1], datatop[:, 1][::-1], mode='same')
    corry /= np.amax(corry)
    xdiff = datatop[1, 0] - datatop[0, 0]
    corrx = np.arange(0.0, len(corry)) * xdiff
    maxpos = np.argmax(corry)
    corrx = corrx - corrx[maxpos]
    autocorr = np.transpose([corrx, corry])
    boo1 = autocorr[:, 0] > 0
    cpeaks = peakdetect(autocorr[boo1], config)
    return autocorr, cpeaks


def cconvolve(xvals, mztab, fwhm, psfun):
    """
    Fast nonlinear convolution algorithm using C dll.
    :param xvals: x-axis
    :param mztab: mztab from make_peaks_mztab()
    :param fwhm: Full width half max of peak
    :param psfun: Peak shape function code (0=Gauss, 1=Lorentzian, 2=Split G/L)
    :return: Convolved output
    """
    stick = np.zeros(len(xvals))
    stick[np.array(mztab[:, 2]).astype(np.int)] = mztab[:, 1]
    cxvals = (c_double * len(xvals))(*xvals)
    cinput = (c_double * len(stick))(*stick)
    cout = (c_double * len(stick))()
    libs.convolve(byref(cxvals), byref(cinput), byref(cout), (c_int)(psfun), (c_double)(fwhm), len(cinput))
    return np.frombuffer(cout)


def conv_peak_shape_kernel(xaxis, psfun, fwhm):
    """
    Creation of an efficient peak shape kernel for circular convolutions.

    Note: peak width must be in units of index not value of x-axis. To get this take fwhm/binsize.
    :param xaxis: Array of x axis values
    :param psfun: Integer

    0=Gaussian
    1=Lorentzian
    2=Split Gaussian/Lorentzian

    :param fwhm: Full width half max of peak (in units matching data indexes not values)
    :return: Peak shape kernel centered at 0.
    """
    limit = int(10 * fwhm)
    length = len(xaxis)
    kernel = np.zeros(length)
    for i in range(-limit, limit):
        if psfun == 0:
            kernel[i % length] = ndis(i, 0, fwhm)
        if psfun == 1:
            kernel[i % length] = ldis(i, 0, fwhm)
        if psfun == 2:
            kernel[i % length] = splitdis(i, 0, fwhm)
    return kernel


def make_peak_shape(xaxis, psfun, fwhm, mid, norm_area=False):
    """
    Make a peak of width fwhm centered at mid for given x axis.

    Note: peak width must be in units of index not value of x-axis. To get this take fwhm/binsize.
    :param xaxis: Array of x axis values
    :param psfun: Integer

    0=Gaussian
    1=Lorentzian
    2=Split Gaussian/Lorentzian

    :param fwhm: Full width half max of peak (in units matching data indexes not values)
    :param mid: Midpoint of peak
    :return: Peak shape centered at midpoint .
    """

    if psfun == 0:
        kernel = ndis(xaxis, mid, fwhm, norm_area=norm_area)
    elif psfun == 1:
        kernel = ldis(xaxis, mid, fwhm, norm_area=norm_area)
    elif psfun == 2:
        kernel = splitdis(xaxis, mid, fwhm, norm_area=norm_area)
    else:
        kernel = xaxis * 0
    return kernel


def gaussfit(xvals, yvals, mguess=None, sguess=0.1, aguess=None, cleanup=True):
    """
    Simple gaussian fitting function.
    :param xvals: X values
    :param yvals: Y values
    :param mguess: Guess for midpoint
    :param sguess: Guess for standard deviation
    :param aguess: Guess for amplitude
    :param cleanup: Boolean Flag, if True, will clean up data to remove background and normalize
    :return: Gaussian fit parameters [mid, sig, amp]
    """
    if cleanup:
        yvals -= np.amin(yvals)
        yvals /= np.amax(yvals)
    if mguess is None:
        mguess = xvals[np.argmax(yvals)]
    if aguess is None:
        aguess = np.amax(yvals)
    guess = [mguess, sguess, aguess]
    fits = curve_fit(ndis_std, xvals, yvals, p0=guess, maxfev=1000000)[0]
    print fits
    return fits


def psfit(x, s, m, a, b, psfun):
    """
    Make peak shape from fit
    :param x: x values
    :param s: fwhm
    :param m: max position
    :param a: amplitude
    :param b: background
    :param psfun: peak shape function integer code
    :return: peak shape fit data
    """
    if psfun == 0:
        return ndis_fit(x, s, m, a, b)
    elif psfun == 1:
        return ldis_fit(x, s, m, a, b)
    elif psfun == 2:
        return splitdis_fit(x, s, m, a, b)


def voigt_fit(xvals, yvals, mguess=0, sguess=0.1, gguess=0, aguess=0, bguess=0):
    """

    """
    guess = [mguess, sguess, gguess, aguess, bguess]
    popt, pcov = curve_fit(voigt, xvals, yvals, p0=guess, maxfev=1000000)
    fitdat = voigt(xvals, popt[0], popt[1], popt[2], popt[3], popt[4])
    return popt, np.sqrt(np.diag(pcov)), fitdat


def fit_peak(xvals, yvals, psfun, midguess, fwhmguess, aguess, bguess):
    """
    Fit peak from xvals and yvals data to defined peak shape function.
    :param xvals: x values of data
    :param yvals: y values of data
    :param psfun: peak shape function integer code
    :param midguess: midpoint guess
    :param fwhmguess: fwhm guess
    :param aguess: amplitude guess
    :param bguess: background guess
    :return: popt, perr, fitdat (optimized parameters [fwhm, mid, a, b], std error of parameters, fit to data)
    """
    guess = [fwhmguess, midguess, aguess, bguess]

    if psfun == 0:
        popt, pcov = curve_fit(ndis_fit, xvals, yvals, p0=guess)
    elif psfun == 1:
        popt, pcov = curve_fit(ldis_fit, xvals, yvals, p0=guess)
    elif psfun == 2:
        popt, pcov = curve_fit(splitdis_fit, xvals, yvals, p0=guess)
    else:
        popt = guess
        pcov = np.ones((len(guess), len(guess)))
        print "Failed"

    fitdat = psfit(xvals, popt[0], popt[1], popt[2], popt[3], psfun)
    return popt, np.sqrt(np.diag(pcov)), fitdat


def isolated_peak_fit(xvals, yvals, psfun, **kwargs):
    """
    Fit an isolated peak to the peak shape model.
    :param xvals: x values of data
    :param yvals: y values of data
    :param psfun: peak shape function integer code
    :param kwargs: keywords (unused)
    :return: fit_array, fit_data (4 x 2 array of (fit, error) for [fwhm,mid,amp,background],fit to data)
    """
    midguess = xvals[np.argmax(yvals)]
    bguess = np.amin(yvals)
    sigguess = weighted_std(xvals, yvals - bguess) * 1
    # Two rounds to guess at area
    if psfun<3:
        testdat = psfit(xvals, sigguess, midguess, 1, bguess, psfun)
        aguess = np.amax(yvals) / np.amax(testdat)
        testdat = psfit(xvals, sigguess, midguess, aguess, bguess, psfun)
        aguess = aguess * np.amax(yvals) / np.amax(testdat)
    else:
        testdat = psfit(xvals, sigguess, midguess, 1, bguess, 0)
        aguess = np.amax(yvals) / np.amax(testdat)
        testdat = psfit(xvals, sigguess, midguess, aguess, bguess, 0)
        aguess = aguess * np.amax(yvals) / np.amax(testdat)
    # Fit it
    if psfun<3:
        fit, err, fitdat = fit_peak(xvals, yvals, psfun, midguess, sigguess, aguess, bguess)
    else:
        fit, err, fitdat = voigt_fit(xvals, yvals, midguess, sigguess, 0, aguess, bguess)
    return np.transpose([fit, err]), fitdat


def correlation_integration(dat1, dat2, alpha=0.01, plot_corr=False, **kwargs):
    """
    Perform MacCoss method (http://pubs.acs.org/doi/pdf/10.1021/ac034790h) of getting peak intensities
    :param dat1: Denominator peak
    :param dat2: Numerator peak
    :param alpha: 1-confidence level
    :param plot_corr: Boolean flag, if True, will plot the results
    :param kwargs: keywords (unused dump)
    :return: slope, cierr, rsquared, slope_std_error * np.sqrt(n), pvalue
    slope of line (ratio of peak areas)
    cierr (confidence interval for error)
    rsquared (R squared value for linear fit)
    slope_std_error * np.sqrt(n) (standard deviation of slope fit)
    pvalue (probability that slope is nonzero)
    """
    n = len(dat1)
    # sigma=0.1
    # r1=np.random.normal(scale=sigma,size=n)
    # r2=np.random.normal(scale=sigma,size=n)

    y1 = dat1[:, 1]  # +r1
    y2 = dat2[:, 1]  # +r2

    outputs = stats.linregress(y1, y2)
    slope, intercept, rvalue, pvalue, slope_std_error = outputs
    rsquared = rvalue ** 2.
    # print rsquared
    df = n - 2  # degrees of freedom
    tval = stats.t.isf(alpha / 2., df)
    # ci=slope + tval*slope_std_error*np.array([-1,1])
    cierr = tval * slope_std_error
    # print ci,cierr

    if slope > 0 and plot_corr:
        print slope, cierr, rsquared, slope_std_error * np.sqrt(n), pvalue
        fitdat = y1 * slope + intercept
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(121)
        plt.scatter(y1, y2)
        plt.plot(y1, fitdat)
        plt.subplot(122)
        plt.plot(y1)
        if pvalue < alpha:
            plt.plot(y2)
        else:
            plt.plot(y2, color="r")
        plt.show()

    return slope, cierr, rsquared, slope_std_error * np.sqrt(n), pvalue


def wide_sse(dat1, dat2, sig):
    """
    Tests the sum of squared errors between two sets of data where the second has been convolved with a Gaussian
    :param dat1: Reference list of intensities.
    :param dat2: Test spectrum list of intensities to broaden.
    :param sig: Sigma for Gaussian filter
    :return: SSE between dat1 and broadedned dat1 (float)
    """
    wdat = filt.gaussian_filter1d(dat2, sig)
    wdat = wdat / np.amax(wdat) * np.amax(dat1)
    return np.sum((wdat - dat1) ** 2)


def broaden(aligned):
    """
    Take a list of aligned peaks and broaden each so that they achieve maximum overlap with the broadest one.
    :param aligned: Array (N x M x 2). N is the number of peaks, M is the number of data points, 2 is (x,y)
    :return: combined, aligned: The combined total of all the aligned and broaded individual peaks.
    """
    # First, normalize everything
    norms = np.array([a[:, 1] / np.amax(a[:, 1]) for a in aligned])
    # Then, find the total area under the normalized peaks
    totals = np.sum(norms, axis=1)

    # Sort the array so that the biggest peak is first
    sortindexes = np.argsort(totals)
    reverse = norms[sortindexes][::-1]

    # For the second peak on, find the minimum deviation as a function of gaussian blur size
    sigs = []
    for i, r in enumerate(reverse):
        if i == 0:
            sigs.append(0)
        else:
            sigrange = np.arange(0, len(r), 1)
            sigdat = [wide_sse(reverse[0], r, s) for s in sigrange]
            minpos = np.argmin(sigdat)
            bestsig = sigrange[minpos]
            sigs.append(bestsig)
    # Reverse everything so that the correct sigmas are restored to their rightful place
    sigs = np.array(sigs)[::-1][sortindexes]

    # Apply the sigmas to the original array
    for i, a in enumerate(aligned):
        a[:, 1] = filt.gaussian_filter1d(a[:, 1], sigs[i])

    # Sum and combine the spectra into a global master
    alignedsum = np.average(aligned[:, :, 1], axis=0)
    combined = np.transpose([aligned[0, :, 0], alignedsum])

    return combined, aligned


if __name__ == "__main__":
    pass
