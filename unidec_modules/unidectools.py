import os
import platform
import sys
import math
import subprocess
import time
import decimal
from bisect import bisect_left
from ctypes import *
from copy import deepcopy
import zipfile
import numpy as np
import scipy.ndimage.filters as filt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy import signal
from scipy import fftpack
import matplotlib.cm as cm
from unidec_modules.mzMLimporter import mzMLimporter
from unidec_modules.fitting import *
from unidec_modules import unidecstructure
import fnmatch

# import unidec_modules.data_reader as data_reader
try:
    import unidec_modules.data_reader as data_reader
except:
    print("Could not import data reader: unidectools")

# from unidec_modules.waters_importer.WatersImporter import WatersDataImporter as WDI
try:
    from unidec_modules.waters_importer.WatersImporter import WatersDataImporter as WDI
except:
    print("Could not import Waters Data Importer")

# from unidec_modules.thermo_reader.ThermoImporter import ThermoDataImporter
try:
    from unidec_modules.thermo_reader.ThermoImporter import ThermoDataImporter
except:
    print("Could not import Thermo Data Importer")

is_64bits = sys.maxsize > 2 ** 32

"""
Loads initial dll called libmypfunc. Will speed up convolutions and a few functions.

If this isn't present, it will print a warning but use the pure python code later on.
"""
dllname = "libmypfunc"
protmass = 1.007276467

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
                print("Unable to find file", testpath)

try:
    libs = cdll.LoadLibrary(dllpath)
except (OSError, NameError):
    print("Failed to load libmypfunc, convolutions in nonlinear mode might be slow")


def get_importer(path):
    if os.path.splitext(path)[1] == ".mzML" or os.path.splitext(path)[1] == ".gz":
        # mzML file
        d = mzMLimporter(path, gzmode=True)
    elif os.path.splitext(path)[1].lower() == ".raw" and not os.path.isdir(path):
        # Thermo Raw File
        try:
            d = ThermoDataImporter(path)
        except:
            d = data_reader.DataImporter(path)
    elif os.path.splitext(path)[1].lower() == ".raw" and os.path.isdir(path):
        # Waters Raw Directory
        d = WDI(path, do_import=False)
    else:
        # Some other file type
        d = data_reader.DataImporter(path)

    return d


# ..........................
#
# Utility Functions
#
# ..................................

def smartdecode(string):
    try:
        string = string.decode()
    except:
        pass
    return string


def commonprefix(args):
    if platform.system() == "Windows":
        sep = "\\"
    else:
        sep = "/"
    return os.path.commonprefix(args).rpartition(sep)[0]


def get_luminance(color, type=2):
    try:
        r = color.Red()
        g = color.Green()
        b = color.Blue()
    except:
        r = color[0]
        g = color[1]
        b = color[2]
    l1 = 0.2126 * r + 0.7152 * g + 0.0722 * b
    l2 = 0.299 * r + 0.587 * g + 0.114 * b
    l3 = np.sqrt(0.299 * r * r + 0.587 * g * g + 0.114 * b * b)
    larray = [l1, l2, l3]
    return larray[type]


def match_files(directory, string):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            files.append(file)
    return np.array(files)


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
    except (TypeError, ValueError, AttributeError):
        print("Error testing emptiness")
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


def string_to_int(s, default=False):
    """
    Try to convert string to integer.
    :param s: String
    :return: Integer if possible. Otherwise, False
    """
    try:
        v = int(s)
        return v
    except (ValueError, TypeError):
        return default


def decimal_formatter(a, template):
    dec = decimal.Decimal(str(template).rstrip("0"))
    dec2 = dec.as_tuple().exponent
    if dec2 == 0:
        label = "{:,}".format(int(a))
    else:
        label = "%g" % round(a, -dec2)
    return label

def round_to_nearest(n, m):
    r = n % m
    return n + m - r if r + r >= m else n - r

def fix_textpos(pos, maxval):
    cutoff = 0.05 * maxval
    posout = []
    for i, p in enumerate(pos):
        newp = p
        if i > 0:
            pastpos = posout[-1]
            if np.abs(p - pastpos) < cutoff:
                if p < pastpos:
                    newp = p + 2 * cutoff
                elif p >= pastpos:
                    newp = p + cutoff
        posout.append(newp)
    return np.array(posout)


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


def safedivide1(a, b):
    if b != 0:
        return a / b
    else:
        return 0


def weighted_std(values, weights):
    """
    Calculate weighted standard deviation.
    :param values: Values
    :param weights: Weights
    :return: Weighted standard deviation.
    """
    average = np.average(values, weights=weights)
    variance = np.average(np.power((np.array(values) - average), 2), weights=weights)  # Fast and numerically precise
    return np.sqrt(variance)


def weighted_avg(values, weights):
    if np.sum(weights) == 0:
        return 0
    return np.average(values, weights=weights)


def mass_weighted_average(value, weights):
    return np.sum(weights * np.power(value, 2)) / np.sum(weights * value)


def polydispersity_index(massdat):
    number_average_molar_mass = weighted_avg(massdat[:, 0], massdat[:, 1])
    mass_average_molar_mass = mass_weighted_average(massdat[:, 0], massdat[:, 1])
    pdi = mass_average_molar_mass / number_average_molar_mass
    return pdi


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


def interp_val(data, target):
    """
    On a given array, interpolate the value of the target
    :param array: data
    :param target: Value
    :return: Interpolated index of value in array.
    """
    f = interp1d(data[:, 0], data[:, 1])
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
    return np.argmin(np.abs(np.array(array) - target))


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
    https://www.pnas.org/content/107/5/2007.long
    :param mass: Mass in Da
    :return: Float, predicted average native charge state
    """
    m = 0.0467
    s = 0.533
    nativez = m * (mass ** s)
    return nativez


def simchargefit(x2):
    fit = [predict_charge(i) for i in x2]
    return np.array(fit)


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


def integrate(data, start, end):
    boo1 = data[:, 0] < end
    boo2 = data[:, 0] > start
    intdat = data[np.all([boo1, boo2], axis=0)]
    integral = np.trapz(intdat[:, 1], x=intdat[:, 0])
    return integral, intdat


def center_of_mass(data, start=None, end=None, relative_cutoff=None, power=1.):
    if start is not None and end is not None:
        boo1 = data[:, 0] < end
        boo2 = data[:, 0] > start
        cutdat = data[np.all([boo1, boo2], axis=0)]
    else:
        cutdat = data

    if relative_cutoff is not None:
        maxval = np.amax(cutdat[:, 1])
        boo3 = cutdat[:, 1] > maxval * relative_cutoff
        cutdat = cutdat[boo3]

    try:
        weightedavg = np.average(cutdat[:, 0], weights=np.power(cutdat[:, 1], float(power)))
        weightstd = weighted_std(cutdat[:, 0], np.power(cutdat[:, 1], float(power)))
        return weightedavg, weightstd
    except ZeroDivisionError:
        # print "Error in center of mass determination. Check the limits:", start, end
        return 0, 0


def stepmax(array, index):
    plus = array[index + 1]
    minus = array[index - 1]
    if plus > minus:
        window = 1
        d = 1
    else:
        window = -1
        d = -1
    val = array[index]
    # newval = array[index + window]
    while array[index + window] > val:
        val = array[index + window]
        window += d
        if index + window >= len(array) or index + window < 0:
            return val
    return val


def localmax(array, start, end):
    start = np.amax([0, start])
    end = np.amin([len(array), end])
    try:
        out = np.amax(array[start:end])
    except (ValueError, TypeError):
        out = 0
    return out


def localmaxpos(data, start, end):
    try:
        boo1 = data[:, 0] < end
        boo2 = data[:, 0] > start
        intdat = data[np.all([boo1, boo2], axis=0)]
        pos = np.argmax(intdat[:, 1])
        return intdat[pos, 0]
    except (ValueError, TypeError):
        return 0


def localmax2(array, index, window):
    start = index - window
    end = index + window
    return localmax(array, start, end)


def simple_extract(data, x, zero_edge=True):
    index = nearest(data[:, 0], x)
    if zero_edge:
        if index != 0 and index != len(data) - 1:
            val = data[index, 1]
        else:
            val = 0
    else:
        val = data[index, 1]
    return val


def limit_data(data, start, end):
    boo1 = data[:, 0] < end
    boo2 = data[:, 0] >= start
    intdat = data[np.all([boo1, boo2], axis=0)]
    return intdat


extractchoices = {0: "Height", 1: "Local Max", 2: "Area", 3: "Center of Mass", 4: "Local Max Position",
                  5: "Center of Mass 50%", 6: "Center of Mass 10%", 7: "Center of Mass Squares",
                  8: "Center of Mass Cubes", 9: "Center of Mass 50% Squares", 10: "Center of Mass 50% Cubes"}


def data_extract(data, x, extract_method, window=None, **kwargs):
    """
    Assumes a sorted array in data[:,0]
    :param data:
    :param x:
    :param extract_method:
    :param window:
    :return:
    """
    if extract_method == 0:
        val = simple_extract(data, x, **kwargs)

    elif extract_method == 1:
        if window is not None:
            if window == 0:
                val = simple_extract(data, x, **kwargs)
            else:
                start = nearest(data[:, 0], (x - window))
                end = nearest(data[:, 0], (x + window))
                val = localmax(data[:, 1], start, end)
        else:
            index = nearest(data[:, 0], x)
            val = stepmax(data[:, 1], index)

    elif extract_method == 2:
        if window is not None:
            start = x - window
            end = x + window
            val, junk = integrate(data, start, end)
        else:
            index = nearest(data[:, 0], x)
            val = data[index, 1]
            print("NEED TO SET INTEGRAL WINDOW!\nUsing Peak Height Instead")

    elif extract_method == 3:
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(data, start, end)
        else:
            val, junk = center_of_mass(data, data[0, 0], data[len(data) - 1, 0])
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 4:
        if window is not None:
            start = x - window
            end = x + window
            val = localmaxpos(data, start, end)
        else:
            val = localmaxpos(data, data[0, 0], data[len(data) - 1, 0])
            print("No window set for local max position!\nUsing entire data range....")

    elif extract_method == 5:
        # Remove data points that fall below 50% threshold
        maxval = np.amax(data[:, 1])
        boo2 = data[:, 1] > maxval * 0.5
        cutdat = data[boo2]
        # Extract from thresholded data
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(cutdat, start, end)
        else:
            val, junk = center_of_mass(cutdat, cutdat[0, 0], cutdat[len(data) - 1, 0])
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 6:
        # Remove data points that fall below 10% threshold
        maxval = np.amax(data[:, 1])
        boo2 = data[:, 1] > maxval * 0.1
        cutdat = data[boo2]
        # Extract from thresholded data
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(cutdat, start, end)
        else:
            val, junk = center_of_mass(cutdat, cutdat[0, 0], cutdat[len(data) - 1, 0])
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 7:
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(data, start, end, power=2)
        else:
            val, junk = center_of_mass(data, data[0, 0], data[len(data) - 1, 0], power=2)
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 8:
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(data, start, end, power=3)
        else:
            val, junk = center_of_mass(data, data[0, 0], data[len(data) - 1, 0], power=3)
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 9:
        # Remove data points that fall below 50% threshold
        maxval = np.amax(data[:, 1])
        boo2 = data[:, 1] > maxval * 0.5
        cutdat = data[boo2]
        # Extract from thresholded data
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(cutdat, start, end, power=2)
        else:
            val, junk = center_of_mass(cutdat, cutdat[0, 0], cutdat[len(data) - 1, 0], power=2)
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 10:
        # Remove data points that fall below 50% threshold
        maxval = np.amax(data[:, 1])
        boo2 = data[:, 1] > maxval * 0.5
        cutdat = data[boo2]
        # Extract from thresholded data
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(cutdat, start, end, power=3)
        else:
            val, junk = center_of_mass(cutdat, cutdat[0, 0], cutdat[len(data) - 1, 0], power=3)
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 11:  # Estimate the peak area
        if window is not None:
            # Calculate height and FWHM
            start = nearest(data[:, 0], (x - window))  # x values should be sorted
            end = nearest(data[:, 0], (x + window))
            data_slice = data[start:end]
            max_index = np.argmax(data_slice[:, 1])
            height = data_slice[max_index, 1]
            # hm1 = data_slice[np.argmin(abs(data_slice[0:max_index, 1] - (height / 2))), 0]
            # hm2 = data_slice[np.argmin(abs(data_slice[max_index:, 1] - (height / 2))) + max_index, 0]
            # fwhm = hm2 - hm1
            try:
                fwhm, psfun, mid = auto_peak_width(data_slice, singlepeak=True)
            except ValueError:
                fwhm = 0
                psfun = 0
            # Calculate estimated area
            gauss_coeff = np.sqrt(np.pi / np.log(2)) / 2
            adjusted_coeff = ((0.5 * gauss_coeff) + (np.pi / 4))
            val = -1  # Error if -1 is returned
            if psfun == 0:  # Gaussian
                val = height * fwhm * gauss_coeff
            elif psfun == 1:  # Lorentzian
                val = height * fwhm * np.pi / 2
            elif psfun == 2:  # Split G/L
                val = height * fwhm * adjusted_coeff

    else:
        val = 0
        print("Undefined extraction choice")
    return val


def data_extract_grid(data, xarray, extract_method=1, window=0):
    igrid = np.zeros_like(xarray)
    dims = igrid.shape
    for j in range(0, dims[0]):
        for k in range(0, dims[1]):
            x = xarray[j, k]
            igrid[j, k] = data_extract(data, x, extract_method=extract_method, window=window)
    return igrid


def normalize_extracts(grid, norm_method=0):
    xlen, ylen = grid.shape
    xarray = range(0, xlen)
    yarray = range(0, ylen)
    # Max
    if norm_method == 1:
        for j in yarray:
            grid[:, j] /= np.amax(grid[:, j])

    # Sum
    if norm_method == 2:
        for j in yarray:
            grid[:, j] /= np.sum(grid[:, j])

    # Peak Max
    if norm_method == 3:
        for j in xarray:
            grid[j] /= np.amax(grid[j])

    # Peak Sum
    if norm_method == 4:
        for j in xarray:
            grid[j] /= np.sum(grid[j])

    return grid


def simple_mass_defect(mass, refmass, centermode=1, normtype=1):
    if refmass == 0:
        print("Error: Reference Mass is 0.")
        return None, None, None, None, None
    kmass = float(mass) / float(refmass)
    if centermode == 1:
        nominalkmass = np.floor(kmass)
    else:
        nominalkmass = np.round(kmass)
    md = kmass - nominalkmass
    if normtype == 1:
        md = md * refmass
    return md


def kendrick_analysis(massdat, kendrickmass, centermode=1, nbins=50, transformmode=1, xaxistype=1):
    # Calculate Defects for Deconvolved Masses
    if kendrickmass == 0:
        print("Error: Kendrick mass is 0.")
        return None, None, None, None, None
    xaxis = massdat[:, 0]
    kmass = np.array(xaxis) / float(kendrickmass)
    if centermode == 1:
        nominalkmass = np.floor(kmass)
    else:
        nominalkmass = np.round(kmass)
    kmdefectexact = kmass - nominalkmass
    # Linearize
    defects = np.linspace(np.amin(kmdefectexact), np.amax(kmdefectexact), int(nbins), endpoint=True)
    nominal = np.unique(nominalkmass)
    m1grid, m2grid = np.meshgrid(nominal, defects, indexing='ij')

    # Get Intensities
    igrid = np.zeros((len(nominal), len(defects)))
    if transformmode == 1:
        # Interpolation
        f = interp1d(massdat[:, 0], massdat[:, 1], bounds_error=False, fill_value=0)
        for j in range(0, len(nominal)):
            for k in range(0, len(defects)):
                nommass = nominal[j]
                defect = defects[k]
                mass = (nommass + defect) * kendrickmass
                intensity = f(mass)
                igrid[j, k] = intensity
    else:
        # Integration
        for j in range(0, len(kmass)):
            nommass = nominalkmass[j]
            defect = kmdefectexact[j]
            pos = nearest(defects, defect)
            pos2 = nearest(nominal, nommass)
            try:
                intensity = massdat[j, 1]
            except IndexError:
                intensity = 0
            igrid[pos2, pos] += intensity
    igrid /= np.amax(igrid)

    # Write Outputs
    if xaxistype == 0:
        factor = 1
    else:
        factor = kendrickmass

    m1grid *= kendrickmass
    m2grid *= factor
    data2 = np.transpose([np.ravel(m1grid), np.ravel(m2grid), np.ravel(igrid)])
    data1 = np.transpose([np.unique(m2grid), np.sum(igrid, axis=0)])

    return data1, data2, m1grid, m2grid, igrid


def solve_for_mass(mz1, mz2, adductmass=1.007276467):
    """
    Simple function to solve for mass from two adjacent charge state peaks in m/z.
    :param mz1: Smaller m/z value (charge = z1)
    :param mz2: Larger m/z value (charge = z2 = z1 + 1)
    :param adductmass: Mass of electrospray adduct (default = Hydrogen)
    :return: mass, z1, z2
    """
    z2 = np.round((mz1 - adductmass) / (mz2 - mz1))
    mass = z2 * mz2
    return mass, z2 + 1, z2


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
                        # print(sin, line)
                        header += 1
                        break
        if header > 0:
            print("Header Length:", header)
    except (ImportError, OSError, AttributeError, IOError) as e:
        print("Failed header test", e)
        header = 0
    return int(header)


def waters_convert(path, config=None, outfile=None):
    if config is None:
        config = unidecstructure.UniDecConfig()
        config.initialize_system_paths()
        print(config.rawreaderpath)
    print(outfile)
    if outfile is None:
        outfile = os.path.join(path, "converted_rawdata.txt")
    call = [config.rawreaderpath, "-i", path, "-o", outfile]
    print(call)
    result = subprocess.call(call)
    print("Conversion Stderr:", result)
    data = np.loadtxt(outfile)
    return data


def waters_convert2(path, config=None, outfile=None):
    data = WDI(path).get_data()

    if outfile is None:
        outfile = os.path.join(path, "converted_rawdata.txt")
    np.savetxt(outfile, data)
    return data


def load_mz_file(path, config=None, time_range=None):
    """
    Loads a text or mzml file
    :param path: File path to load
    :param config: UniDecConfig object
    :return: Data array
    """
    if config is None:
        extension = os.path.splitext(path)[1]
    else:
        extension = config.extension.lower()

    if not os.path.isfile(path):
        if os.path.isdir(path) and os.path.splitext(path)[1].lower() == ".raw":
            try:
                outfile = os.path.splitext(path)[0] + "_rawdata.txt"
                print("Trying to convert Waters File:", outfile)
                data = waters_convert2(path, config, outfile=outfile)
            except:
                print("Attempted to convert Waters Raw file but failed")
                raise IOError
        elif os.path.isdir(path) and os.path.splitext(path)[1].lower() == ".d":
            try:
                print("Trying to convert Agilent File:", path)
                data = data_reader.DataImporter(path).get_data(time_range=time_range)
                txtname = path[:-2] + ".txt"
                np.savetxt(txtname, data)
                print("Saved to:", txtname)
            except:
                print("Attempted to convert Agilent Raw file but failed")
                raise IOError
        else:
            print("Attempted to open:", path)
            print("\t but I couldn't find the file...")
            raise IOError
    else:
        if extension == ".txt":
            data = np.loadtxt(path, skiprows=header_test(path))
            # data = np.loadtxt(path, skiprows=header_test(path, delimiter=","), delimiter=",")
        elif extension == ".csv":
            data = np.loadtxt(path, delimiter=",", skiprows=1, usecols=(0, 1))
        elif extension == ".mzml":
            data = mzMLimporter(path).get_data(time_range=time_range)
            txtname = path[:-5] + ".txt"
            np.savetxt(txtname, data)
            print("Saved to:", txtname)
        elif extension.lower() == ".raw":
            data = ThermoDataImporter(path).get_data(time_range=time_range)
            txtname = path[:-4] + ".txt"
            np.savetxt(txtname, data)
            print("Saved to:", txtname)
        else:
            try:
                data = np.loadtxt(path, skiprows=header_test(path))
            except IOError:
                print("Failed to open:", path)
                data = None
    return data


def zipdir(path, zip_handle):
    """
    Zips all the files in the path into the zip_handle
    :param path: Directory path
    :param zip_handle: Handle of the zip file that is being created
    :return: None
    """
    files = os.scandir(path)
    for f in files:
        if f.is_file():
            zip_handle.write(f.path, arcname=os.path.relpath(f.path, path))  # compress_type=zipfile.ZIP_DEFLATED)


def zip_folder(save_path, directory=None):
    """
    Zips a directory specified by save_path into a zip file for saving.
    :param save_path: Path to save to zip
    :return: None
    """
    if directory is None:
        directory = os.getcwd()
    print("Zipping directory:", directory)
    zipf = zipfile.ZipFile(save_path, 'w')
    zipdir(directory, zipf)
    zipf.close()
    print("File saved to: " + str(save_path))


def dataexport(datatop, fname):
    try:
        np.savetxt(fname, datatop, fmt='%f')
    except:
        path = os.path.join(os.getcwd(), fname)
        path = "\\\\?\\%s" % path
        np.savetxt(path, datatop, fmt='%f')
        print("NOTE: Your path length might exceed the limit for Windows. Please shorten your file name.")
    pass


def savetxt(fname, datatop):
    try:
        np.savetxt(fname, datatop)
    except:
        path = os.path.join(os.getcwd(), fname)
        path = "\\\\?\\%s" % path
        np.savetxt(path, datatop)
        print("NOTE: Your path length might exceed the limit for Windows. Please shorten your file name.")


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


def auto_peak_width(datatop, psfun=None, singlepeak=False):
    maxpos = np.argmax(datatop[:, 1])
    maxval = datatop[maxpos, 0]

    # TODO: This is potentially dangerous if nonlinear!
    ac, cpeaks = autocorr(datatop)
    if singlepeak or not isempty(cpeaks):
        if not singlepeak:
            sig = cpeaks[0, 0] / 2.
            boo1 = datatop[:, 0] < maxval + sig
            boo2 = datatop[:, 0] > maxval - sig
            boo3 = np.all([boo1, boo2], axis=0)
            isodat = datatop[boo3]

            if len(isodat) < 6:
                sig = cpeaks[0, 0]
                boo1 = datatop[:, 0] < maxval + sig
                boo2 = datatop[:, 0] > maxval - sig
                boo3 = np.all([boo1, boo2], axis=0)
                isodat = datatop[boo3]
                if len(isodat) < 6:
                    try:
                        sig = cpeaks[1, 0]
                        boo1 = datatop[:, 0] < maxval + sig
                        boo2 = datatop[:, 0] > maxval - sig
                        boo3 = np.all([boo1, boo2], axis=0)
                        isodat = datatop[boo3]
                    except:
                        pass
                    if len(isodat) < 6:
                        print("Warning: Very small range selected for auto peaks width:", sig, isodat)
        else:
            isodat = datatop

        fits = np.array([isolated_peak_fit(isodat[:, 0], isodat[:, 1], i) for i in range(0, 3)], dtype="object")

        errors = [np.sum(np.array(isodat[:, 1] - f) ** 2.) for f in fits[:, 1]]
        if psfun is None:
            psfun = np.argmin(errors)
        fit, fitdat = fits[psfun]

        fwhm = fit[0, 0]
    else:
        fwhm = 0
    fwhm = np.round(fwhm, 5)
    return fwhm, psfun, maxval


def auto_noise_level(datatop, buffer=10):
    ndat = np.ravel([datatop[:buffer, 1], datatop[-buffer:, 1]])
    std = np.std(ndat)
    mean = np.mean(ndat)
    return mean + 5 * std


def noise_level2(ticdat, percent=0.75, number_stddevs=3):
    sdat = np.sort(ticdat[:, 1])
    index = round(len(sdat) * percent)
    cutoff = sdat[index]
    below = ticdat[:, 1] <= cutoff
    noise = ticdat[below, 1]
    noise = np.std(noise) * number_stddevs + np.mean(noise)
    return noise


def average_bin_size(datatop):
    diffs = np.diff(datatop[:, 0])
    return np.average(diffs)


def cal_data(datatop, poly_coeff=None):
    if poly_coeff is None:
        poly_coeff = [1.54101412e-08, 1.00077531e+00, 1.86397570e-02]
        # This is a default value for when our neg mode went crazy
    datatop[:, 0] = np.polyval(poly_coeff, datatop[:, 0])
    return datatop


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
    buff = int(buff)
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
    mins = list(range(0, length))
    indexes = list(range(0, length))
    for i in indexes:
        mins[i] = np.amin(datatop[int(max([0, i - abs(buff)])):int(min([i + abs(buff), length])), 1])
    background = filt.gaussian_filter(mins, abs(buff) * 2)
    datatop[:, 1] = datatop[:, 1] - background
    return datatop


def calc_local_mins(data, w):
    start = np.amin(data[:, 0])
    stop = np.amax(data[:, 0])
    windows = np.arange(start, stop, step=w)

    localmins = []
    for winstart in windows:
        chopdata = datachop(data, winstart, winstart + w)
        localmin = np.amin(chopdata[:, 1])
        localminpos = chopdata[np.argmin(chopdata[:, 1]), 0]
        localmins.append([localminpos, localmin])
    return np.array(localmins)


def polynomial_background_subtract(datatop, polynomial_order=4, width=20, cutoff_percent=0.25):
    starting_max = np.amax(datatop[:, 1])

    mindat = calc_local_mins(datatop, width)
    coeff = np.polyfit(mindat[:, 0], mindat[:, 1], polynomial_order)
    fitdat = np.polyval(coeff, mindat[:, 0])
    sse = (mindat[:, 1] - fitdat) / fitdat
    good = sse < cutoff_percent
    bad = np.logical_not(good)
    baddat = mindat[bad]
    gooddat = mindat[good]
    while len(baddat) > 1:
        coeff = np.polyfit(gooddat[:, 0], gooddat[:, 1], polynomial_order)
        fitdat = np.polyval(coeff, gooddat[:, 0])
        sse = (gooddat[:, 1] - fitdat) / fitdat
        good = sse < cutoff_percent
        bad = sse >= cutoff_percent
        baddat = gooddat[bad]
        gooddat = gooddat[good]
    coeff = np.polyfit(gooddat[:, 0], gooddat[:, 1], polynomial_order)
    background = np.clip(np.polyval(coeff, datatop[:, 0]), 0, starting_max)
    datatop[:, 1] = np.clip(datatop[:, 1] - background, 0, starting_max)
    return datatop


def savgol(ydata, window=None, order=2):
    # Perform a series of checks to make sure we have enough data points and a sufficient window
    if len(ydata) < 2:
        return ydata
    if len(ydata) < int(order) * 2 + 1:
        order = round((len(ydata) - 1) / 2)
        window = None
    if window is None or window < order:
        window = int(order) * 2 + 1
    else:
        window = round((window - 1) / 2.) * 2 + 1
    # Execute the filter
    return signal.savgol_filter(ydata, window, order)


def savgol_background_subtract(datatop, width, cutoff_percent=0.25):
    starting_max = np.amax(datatop[:, 1])
    sgorder = 2
    buff = sgorder * 2 + 1

    datatop[:, 1] = savgol(datatop[:, 1], window=buff, order=sgorder)

    meddat = calc_local_mins(datatop, width)
    # It is important to keep the end values to allow for interpolation in the end
    firstpoint = [datatop[0]]
    lastpoint = [datatop[len(datatop) - 1]]
    if firstpoint not in meddat:
        meddat = np.concatenate((firstpoint, meddat))
    if lastpoint not in meddat:
        meddat = np.concatenate((meddat, lastpoint))

    def keep_the_extremes(good, bad):
        good[0] = True
        bad[0] = False
        good[len(good) - 1] = True
        bad[len(good) - 1] = False

    sgsmooth = savgol(meddat[:, 1], window=buff, order=sgorder)
    sse = (meddat[:, 1] - sgsmooth) / sgsmooth
    good = sse < cutoff_percent
    bad = np.logical_not(good)
    keep_the_extremes(good, bad)
    baddat = meddat[bad]
    gooddat = meddat[good]

    # import matplotlib.pyplot as plt
    # plt.plot(datatop[:, 0], datatop[:, 1])
    # plt.plot(gooddat[:, 0], gooddat[:, 1], marker="o", linestyle="", color="k")
    while len(baddat) > 1 and len(gooddat) > int(buff):
        sgsmooth = savgol(gooddat[:, 1], window=buff, order=sgorder)
        sse = (gooddat[:, 1] - sgsmooth) / sgsmooth
        good = sse < cutoff_percent
        bad = sse >= cutoff_percent
        keep_the_extremes(good, bad)
        baddat = gooddat[bad]
        gooddat = gooddat[good]
        # plt.plot(gooddat[:,0],gooddat[:,1],marker="o",linestyle="")
    sgsmooth = savgol(gooddat[:, 1], window=buff, order=sgorder)
    f = interp1d(gooddat[:, 0], sgsmooth, kind="linear", fill_value=0, bounds_error=True)
    background = f(datatop[:, 0])
    background = savgol(background, window=buff, order=sgorder)
    background = np.clip(background, 0, starting_max)
    datatop[:, 1] = np.clip(datatop[:, 1] - background, 0, starting_max)
    # plt.plot(datatop[:, 0], background, color="r")
    # plt.plot(gooddat[:, 0], sgsmooth, color="r", marker="o", linestyle="")
    # plt.plot(datatop[:, 0], datatop[:, 1],color = "k")
    # plt.show()

    return datatop


def gaussian_backgroud_subtract(datatop, sig):
    background = deepcopy(datatop)
    # background[:,1]=np.log(background[:,1])
    gsmooth(background, sig)
    datatop[:, 1] -= background[:, 1]
    normalize(datatop)
    datatop[datatop[:, 1] < 0, 1] = 0
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


def intensitythresh_del(datatop, thresh):
    """
    Sets an intensity threshold. Everything below the threshold is set to 0.
    :param datatop: Data array
    :param thresh: Threshold value
    :return: Thresholded array.
    """
    above = datatop[:, 1] > thresh
    datatop = datatop[above]
    return datatop


def intensitythresh_sub(datatop, thresh):
    """
    Sets an intensity threshold. Everything below the threshold is set to 0. Everything else is subtracted by threshold.
    :param datatop: Data array
    :param thresh: Threshold value
    :return: Thresholded array.
    """
    datatop[:, 1] -= thresh * np.amax(datatop[:, 1])
    belowint = datatop[:, 1] < 0
    datatop[belowint, 1] = 0
    return datatop


def remove_noise(datatop, percent=None):
    l1 = len(datatop)
    if percent is None:
        cutoff = np.average(datatop[:, 1]) * 2
        print("Removing points below the mean intensity times 2:", cutoff)
        # datatop = intensitythresh_del(datatop, cutoff)
    else:
        sdat = np.sort(datatop[:, 1])
        index = round(l1 * percent / 100.)
        cutoff = sdat[index]
    datatop = intensitythresh(datatop, cutoff)
    # datatop = remove_middle_zeros(datatop)
    print(l1, len(datatop))
    return datatop


def gsmooth(datatop, sig):
    """
    Smooths the data with a Gaussian filter of width sig.
    :param datatop: Data array
    :param sig: Width of Gaussian Array
    :return: Smoothed Data
    """
    print(len(datatop), sig)
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
        # print "Integrating"
    else:
        newdat = linterpolate(datatop, intx)
        # print "Interpolating"
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
                         range(0, len(datatop), num_compressed)])


def removeduplicates(datatop):
    """
    Cleans up data that may have the same x value multiple times.
    Each nonunique x gets merged into a single unique x with the sum of all of the nonunique intensities.
    :param datatop: Data array (N x 2)
    :return: Data with unique x values
    """
    testunique = np.unique(datatop[:, 0])
    if len(testunique) != len(datatop):
        print("Removing Duplicates")
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


def fake_log(data):
    non_zero_min = np.amin(data[data > 0])
    return np.log10(np.clip(data, non_zero_min, np.amax(data)))


def remove_middle_zeros(data):
    boo1 = data[1:len(data) - 1, 1] != 0
    boo2 = data[:len(data) - 2, 1] != 0
    boo3 = data[2:len(data), 1] != 0
    boo4 = np.logical_or(boo1, boo2)
    boo5 = np.logical_or(boo3, boo4)
    boo6 = np.concatenate(([True], boo5, [True]))
    return data[boo6]


def dataprep(datatop, config, peaks=True, intthresh=True):
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
    thresh = config.intthresh
    subtype = config.subtype
    va = config.detectoreffva
    linflag = config.linflag
    redper = config.reductionpercent
    # Crop Data
    data2 = datachop(deepcopy(datatop), newmin, newmax)
    if len(data2) == 0:
        print("Error: m/z range is too small. No data fits the range.")
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
    elif subtype == 4 and buff != 0:
        data2 = polynomial_background_subtract(data2, buff)
    elif subtype == 5 and buff != 0:
        data2 = savgol_background_subtract(data2, buff)
    elif subtype == 3 and buff != 0:
        data2 = intensitythresh_sub(data2, buff)
        pass
    elif buff == 0:
        pass
    else:
        print("Background subtraction code unsupported", subtype, buff)

    # Scale Adjustment
    print(config.intscale, config.intscale == "Square Root")
    if config.intscale == "Square Root":
        data2[:, 1] = np.sqrt(data2[:, 1])
        print("Square Root Scale")
    elif config.intscale == "Logarithmic":
        data2[:, 1] = fake_log(data2[:, 1])
        data2[:, 1] -= np.amin(data2[:, 1])
        print("Log Scale")

    # Data Reduction
    if redper > 0:
        data2 = remove_noise(data2, redper)

    if config.datanorm == 1:
        # Normalization
        data2 = normalize(data2)

    # Intensity Threshold
    if thresh > 0 and intthresh:
        print(len(data2))
        data2 = intensitythresh(data2, thresh)  # thresh
        print("Intensity Threshold Applied:", thresh, len(data2), np.amax(data2[:, 1]))
        # data2=data2[data2[:,1]>0]
    else:
        data2 = intensitythresh(data2, 0)  # thresh
        # data2=data2[data2[:,1]>0]

    if redper < 0 and peaks:
        # widths = [1, 2, 4, 8, 16]
        # peakind = signal.find_peaks_cwt(data2[:,1], widths)
        # data2 = data2[peakind]
        data2 = peakdetect(data2, window=-redper, threshold=thresh)
        print(data2)

    if linflag == 2:
        try:
            data2 = remove_middle_zeros(data2)
        except:
            pass
        pass

    return data2


# ......................................................
#
# UniDec Functions
#
# ............................................................................


def unidec_call(config, silent=False, **kwargs):
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

    call = [config.UniDecPath, str(config.confname)]

    if config.autotune:
        call.append("-autotune")

    if silent:
        out = subprocess.call(call, stdout=subprocess.PIPE)
    else:
        out = subprocess.call(call)
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
    for i in range(0, length):
        if data[i, 1] > maxval * threshold:
            start = i - window
            end = i + window
            if start < 0:
                start = 0
            if end > length:
                end = length
            start = int(start)
            end = int(end) + 1
            testmax = np.amax(data[start:end, 1])
            if data[i, 1] == testmax and np.all(data[i, 1] != data[start:i, 1]):
                peaks.append([data[i, 0], data[i, 1]])

    return np.array(peaks)


def peakdetect_nonlinear(data, config=None, window=1, threshold=0):
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
        window = config.peakwindow
        threshold = config.peakthresh
    peaks = []
    length = len(data)
    maxval = np.amax(data[:, 1])
    for i in range(0, length):
        if data[i, 1] > maxval * threshold:
            start = data[i, 0] - window
            end = data[i, 0] + window
            isodat = datachop(data, start, end)
            testmax = np.amax(isodat[:, 1])
            index = nearest(isodat[:, 0], data[i, 0])
            if data[i, 1] == testmax and np.all(data[i, 1] != isodat[:index, 1]):
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


def make_peaks_mztab(mzgrid, pks, adductmass, index=None):
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
    ftab = [interp1d(xvals, newgrid[:, k]) for k in range(0, ylen)]
    mztab = [[makespecfun(i, k, pks.masses, adductmass, yvals, xvals, ftab, xmax, xmin) for k in range(0, ylen)] for i
             in range(0, plen)]
    if index is None:
        for i in range(0, plen):
            pks.peaks[i].mztab = np.array(mztab[i])
    else:
        for i in range(0, plen):
            pks.peaks[i].mztab.append(np.array(mztab[i]))
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


def make_peaks_mztab_spectrum(mzgrid, pks, data2, mztab, index=None):
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

    if index is None:
        mztab2[:, :, 1] = [[data2[int(pks.peaks[i].mztab[k, 2]), 1] for k in range(0, zlen)] for i in range(0, plen)]
        for i in range(0, plen):
            pks.peaks[i].mztab2 = np.array(mztab2[i])
    else:
        mztab2[:, :, 1] = [[data2[int(pks.peaks[i].mztab[index][k, 2]), 1] for k in range(0, zlen)] for i in
                           range(0, plen)]
        for i in range(0, plen):
            pks.peaks[i].mztab2.append(np.array(mztab2[i]))

    return mztab2


def fix_double_zeros(processed_data, new_data):
    for i in range(len(processed_data) - 1):
        if processed_data[i, 1] == 0 and processed_data[i + 1, 1] == 0:
            new_data[i] = 0
            new_data[i + 1] = 0
    return new_data


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

    if config.mzsig == 0:
        stickdat = [sticks_only(p.mztab, xvals) for p in pks.peaks]
    else:
        if config.linflag != 2:
            peakwidth = config.mzsig / config.mzbins
            kernel = conv_peak_shape_kernel(xvals, config.psfun, peakwidth)
            stickdat = [stickconv(p.mztab, kernel) for p in pks.peaks]
        else:
            try:
                stickdat = [cconvolve(xvals, p.mztab, config.mzsig, config.psfun) for p in pks.peaks]
            except (OSError, TypeError, NameError, AttributeError):
                stickdat = [nonlinstickconv(xvals, p.mztab, config.mzsig, config.psfun) for p in pks.peaks]

    pks.composite = np.zeros(xlen)
    for i in range(0, pks.plen):
        try:
            sd = fix_double_zeros(processed_data, stickdat[i])
        except:
            sd = stickdat[i]
        pks.peaks[i].stickdat = sd
        pks.composite += np.array(sd)
    pks.convolved = True
    return np.array(stickdat)  # Note: this is not processed with fix_double_zeros and may cause issues...


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
    bool1 = [np.abs(xvals - xvals[i]) < window for i in range(0, xlen)]
    kernels = np.array([make_peak_shape(-xvals[bool1[i]], psfun, fwhm, -xvals[i]) for i in range(0, xlen)])
    output = np.array([np.sum(kernels[i] * stick[bool1[i]]) for i in range(0, xlen)])
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


def sticks_only(mztab, kernel):
    """
    Make stick spectrum and then convolves with kernel.
    :param mztab: mztab from make_peaks_mztab
    :param kernel: peak shape kernel
    :return: Convolved output
    """
    xlen = len(kernel)
    temp = np.zeros(xlen)
    temp[np.array(mztab[:, 2]).astype(np.int)] = mztab[:, 1]
    return temp


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


# ..............................................
#
# Matching Functions
#
# ..................................................

def lengths(array):
    top = []
    num = len(array)
    for i in range(0, num):
        start = int(array[i, 2])
        end = int(array[i, 3])
        top.append(end + 1 - start)
    return top


# TODO: Merge into make_isolated_matches
def combine(array2):
    lens = lengths(array2)
    tup = tuple(lens)
    startindex = array2[:, 2]
    startindex = startindex.astype(np.int)
    basemass = array2[:, 0]
    basemass = basemass.astype(np.float)
    omass = array2[:, 1]
    omass = omass.astype(np.float)
    finlist = []
    for index in np.ndindex(tup):
        total = np.sum((index + startindex) * omass + basemass)
        if total > 0:
            finlist.append(total)
    return finlist


def combine_all(array2):
    lens = lengths(array2)
    tup = tuple(lens)
    startindex = array2[:, 2]
    startindex = startindex.astype(np.int)
    basemass = array2[:, 0]
    basemass = basemass.astype(np.float)
    omass = array2[:, 1]
    omass = omass.astype(np.float)
    names = array2[:, 4]
    finlist = []
    namelist = []
    for index in np.ndindex(tup):
        name = ""
        for i in range(0, len(index)):
            val = index[i] + startindex[i]
            if val > 0:
                if names[i] == "":
                    name = name + str(val)
                else:
                    name = name + str(val) + "[" + names[i] + "] "
            else:
                pass
        total = np.sum((index + startindex) * omass + basemass)
        if total > 0:
            finlist.append(total)
            namelist.append(name)
            # print index,name
    return finlist, namelist


def make_isolated_match(oligos):
    oligomasslist = []
    oligonames = []
    for i in range(0, len(oligos)):
        start = int(oligos[i][2])
        end = int(oligos[i][3])
        for j in range(start, end + 1, 1):
            newmass = float(oligos[i][0]) + j * float(oligos[i][1])
            if newmass > 0:
                oligomasslist.append(newmass)
                if j > 0 or oligos[i][4] == "":
                    if oligos[i][4] == "":
                        oligonames.append(str(j))
                    else:
                        oligonames.append(str(j) + "[" + oligos[i][4] + "]")
                else:
                    oligonames.append("")
                    # self.oligonames.append(str(j)+""+oligos[i][4])
    oligomasslist = np.array(oligomasslist)
    return oligomasslist, oligonames


def make_all_matches(oligos):
    if len(oligos) > 1:
        oligos = np.array(oligos)
        oligomasslist, oligonames = combine_all(oligos)
        oligomasslist = np.array(oligomasslist)
    else:
        oligomasslist, oligonames = make_isolated_match(oligos)
    return oligomasslist, oligonames


def match(pks, oligomasslist, oligonames, tolerance=None):
    matches = []
    errors = []
    peaks = []
    names = []
    for i in range(0, pks.plen):
        p = pks.peaks[i]
        target = p.mass
        nearpt = nearestunsorted(oligomasslist, target)
        match = oligomasslist[nearpt]
        error = target - match
        if tolerance is None or np.abs(error) < tolerance:
            name = oligonames[nearpt]
        else:
            name = ""

        p.label = name
        p.match = match
        p.matcherror = error
        matches.append(match)
        errors.append(error)
        peaks.append(target)
        names.append(name)
    matchlist = [peaks, matches, errors, names]
    return matchlist


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
    return signal.fftconvolve(a, np.roll(b, (len(b)) // 2 - 1 + len(b) % 2), mode="same")


def FD_gauss_wavelet(length, width):
    xvals = np.arange(0, length).astype(np.float)
    mu = np.round(length / 2)
    sigma = width
    return [stats.norm.pdf(x, mu, sigma) * (mu - x) / sigma ** 2 for x in xvals]


def single_cwt(a, width, wavelet_type="Ricker"):
    """
    Perform a continuous wavelet transform at a single defined width.
    :param a: 1D numpy array of data (length N)
    :param width: Width of transform defined with respect to the data points in a
    :param wavelet_type: Type of wavelet. Either "Ricker" (Mexican Hat) or "Morlet" (Gabor)
    :return: cwt, wavelet
    Continous wavelet transform and wavelet of choice, both as numpy arrays of length N.
    """
    if wavelet_type == "Morlet":
        wdat = signal.morlet(len(a), width)
    elif wavelet_type == "1DG":
        wdat = FD_gauss_wavelet(len(a), width)
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
    if wavelet_type == "Morlet":
        wavelet = signal.morlet
    elif wavelet_type == "1DG":
        wavelet = FD_gauss_wavelet
    else:
        wavelet = signal.ricker
    return signal.cwt(a, wavelet, widths)


def autocorr(datatop, config=None):
    """
    Note: ASSUMES LINEARIZED DATA
    :param datatop: 1D data
    :param config: Config file (optional file)
    :return: Autocorr spectrum, peaks in autocorrelation.
    """
    corry = signal.fftconvolve(datatop[:, 1], datatop[:, 1][::-1], mode='same')
    corry /= np.amax(corry)
    maxpos1 = np.argmax(datatop[:, 1])
    start = np.amax([maxpos1 - len(datatop) / 10, 0])
    end = np.amin([len(datatop) - 1, maxpos1 + len(datatop) / 10])
    cutdat = datatop[int(start):int(end)]
    if len(cutdat) < 20:
        cutdat = datatop
    # cutdat=datatop # Other old
    xdiff = np.mean(cutdat[1:, 0] - cutdat[:len(cutdat) - 1, 0])  # Less dangerous but still dangerous when non-linear
    # xdiff = datatop[1, 0] - datatop[0, 0] #OLD
    corrx = np.arange(0.0, len(corry)) * xdiff
    maxpos = np.argmax(corry)
    corrx = corrx - corrx[maxpos]
    autocorr = np.transpose([corrx, corry])
    boo1 = autocorr[:, 0] > xdiff
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
    elif psfun == 3:
        kernel = ndis(xaxis, mid, fwhm, norm_area=norm_area)
    else:
        kernel = xaxis * 0
    return kernel


def fft(data, phases=False):
    # X-axis
    massdif = data[1, 0] - data[0, 0]
    # fvals = np.fft.fftfreq(len(data), d=massdif)
    fvals = fftpack.fftfreq(len(data), d=massdif)
    # Y-axis
    # fft = np.fft.fft(data[:, 1])
    fft = fftpack.fft(data[:, 1])
    # fft = pyfftw.interfaces.numpy_fft.fft(data[:,1])
    # Cleanup
    aftdat = np.transpose([fvals, np.abs(fft)])
    aftdat[:, 1] /= np.amax(aftdat[:, 1])
    if not phases:
        return aftdat
    if phases:
        phases = np.arctan2(fft.imag, fft.real)
        return aftdat, phases


def fft_diff(data, diffrange=[500., 1000.]):
    fftdat = fft(data)
    ftrange = [1. / diffrange[1], 1. / diffrange[0]]
    boo1 = fftdat[:, 0] < ftrange[1]
    boo2 = fftdat[:, 0] > ftrange[0]
    boo3 = np.all([boo1, boo2], axis=0)
    ftext = fftdat[boo3]
    maxpos = localmaxpos(ftext, ftrange[0], ftrange[1])

    # fit, err, fitdat = voigt_fit(ftext[:, 0], ftext[:, 1], np.average(ftrange), np.average(ftrange) / 10., 0, 1, 0)
    return 1. / maxpos, ftext


def pad_data(linear_mzdata, pad_until=50000):
    mzdiff = linear_mzdata[1, 0] - linear_mzdata[0, 0]
    maxmz = np.amax(linear_mzdata[:, 0])
    paddedmz = np.arange(maxmz + mzdiff, pad_until, mzdiff)
    paddedint = np.zeros(len(paddedmz))
    paddat = np.transpose([paddedmz, paddedint])
    return np.concatenate((linear_mzdata, paddat))


def pad_data_length(linear_mzdata, pad_until_length=50000):
    mzdiff = linear_mzdata[1, 0] - linear_mzdata[0, 0]
    maxmz = np.amax(linear_mzdata[:, 0])
    paddedmz = np.arange(maxmz + mzdiff, maxmz + mzdiff * (pad_until_length - len(linear_mzdata) + 1), mzdiff)
    paddedint = np.zeros(len(paddedmz))
    paddat = np.transpose([paddedmz, paddedint])
    return np.concatenate((linear_mzdata, paddat))


def pad_two_power(data):
    mzdiff = data[1, 0] - data[0, 0]
    maxmz = np.amax(data[:, 0])
    n_original = len(data)
    n_power_of_2 = 2 ** int(math.ceil(math.log(n_original, 2)))
    n_pad = n_power_of_2 - n_original
    z = np.zeros(n_pad)
    x = np.arange(n_pad) * mzdiff + maxmz
    padding = np.transpose([x, z])
    return np.concatenate((data, padding))


def double_fft_diff(mzdata, diffrange=None, binsize=0.1, pad=None, preprocessed=False):
    tstart = time.perf_counter()
    if diffrange is None:
        diffrange = [600, 900]
    if not preprocessed:
        mzdata = linearize(mzdata, binsize, 3)
        mzdata = pad_two_power(mzdata)
        if pad is not None:
            mzdata = pad_data(mzdata, pad_until=pad)

    fftdat = fft(mzdata)
    fft2 = fft(fftdat)

    maxpos = localmaxpos(fft2, diffrange[0], diffrange[1])
    boo1 = fft2[:, 0] < diffrange[1]
    boo2 = fft2[:, 0] > diffrange[0]
    boo3 = np.all([boo1, boo2], axis=0)
    ftext2 = fft2[boo3]

    return maxpos, ftext2


def fft_process(mzdata, diffrange=None, binsize=0.1, pad=None, preprocessed=False):
    tstart = time.perf_counter()
    if diffrange is None:
        diffrange = [600, 900]
    if not preprocessed:
        mzdata = linearize(mzdata, binsize, 3)
        mzdata = pad_two_power(mzdata)
        if pad is not None:
            mzdata = pad_data(mzdata, pad_until=pad)

    fftdat = fft(mzdata)
    fft2 = fft(fftdat)

    maxpos = localmaxpos(fft2, diffrange[0], diffrange[1])
    boo1 = fft2[:, 0] < diffrange[1]
    boo2 = fft2[:, 0] > diffrange[0]
    boo3 = np.all([boo1, boo2], axis=0)
    ftext2 = fft2[boo3]

    return maxpos, ftext2, fftdat, fft2


def windowed_fft(data, mean, sigma, diffrange=None, norm=True):
    if diffrange is None:
        diffrange = [740, 770]
    window = ndis_std(data[:, 0], mean, sigma)
    newdata = deepcopy(data)
    newdata[:, 1] = newdata[:, 1] * window
    maxpos, fft2 = double_fft_diff(newdata, diffrange=diffrange, preprocessed=True)
    if norm:
        factor = np.sum(newdata[:, 1])
        fft2[:, 1] *= factor
    fft2[:, 1] -= np.amin(fft2[:, 1])
    return maxpos, fft2


def windowed_fft_single(data, mean, sigma, diffrange=None, norm=True):
    if diffrange is None:
        diffrange = [0, 10]
    window = ndis_std(data[:, 0], mean, sigma)
    newdata = deepcopy(data)
    newdata[:, 1] = newdata[:, 1] * window
    maxpos, fft2 = fft_diff(newdata, diffrange=diffrange)
    if norm:
        factor = np.sum(newdata[:, 1])
        fft2[:, 1] *= factor
    fft2[:, 1] -= np.amin(fft2[:, 1])
    return maxpos, fft2


def win_fft_grid(rawdata, binsize, wbin, window_fwhm, diffrange, norm=True):
    # Prepare data
    mindat = np.amin(rawdata[:, 0])
    maxdat = np.amax(rawdata[:, 0])

    mzdata = linearize(rawdata, binsize, 3)
    mzdata = pad_two_power(mzdata)

    xvals = np.arange(mindat, maxdat, wbin)

    results = np.array([windowed_fft(mzdata, x, window_fwhm, diffrange=diffrange, norm=norm)[1] for x in xvals])

    intdat = results[:, :, 1]
    yvals = np.unique(results[:, :, 0])
    xgrid, ygrid = np.meshgrid(xvals, yvals, indexing="ij")
    out = np.transpose([np.ravel(xgrid), np.ravel(ygrid), np.ravel(intdat)])
    return out


def win_fft_grid_single(rawdata, binsize, wbin, window_fwhm, diffrange, norm=True):
    # Prepare data
    mindat = np.amin(rawdata[:, 0])
    maxdat = np.amax(rawdata[:, 0])

    mzdata = linearize(rawdata, binsize, 3)
    mzdata = pad_two_power(mzdata)

    xvals = np.arange(mindat, maxdat, wbin)

    results = np.array([windowed_fft_single(mzdata, x, window_fwhm, diffrange=diffrange, norm=norm)[1] for x in xvals])

    intdat = results[:, :, 1]
    yvals = np.unique(results[:, :, 0])
    xgrid, ygrid = np.meshgrid(xvals, yvals, indexing="ij")
    out = np.transpose([np.ravel(xgrid), np.ravel(ygrid), np.ravel(intdat)])
    return out


def win_fft_diff(rawdata, binsize=0.05, sigma=1000, diffrange=None, norm=True):
    smoothdata = deepcopy(rawdata)
    smoothdata = gsmooth(smoothdata, 1000)
    mzdata = linearize(rawdata, binsize, 3)
    mzdata = pad_two_power(mzdata)
    maxpos = smoothdata[np.argmax(smoothdata[:, 1]), 0]
    maxdiff, fftdat = windowed_fft(mzdata, maxpos, sigma, diffrange=diffrange, norm=norm)
    print("Difference:", maxdiff)
    return maxdiff, fftdat


def windowed_autocorr(data, mean, sigma, diffrange=None):
    if diffrange is None:
        diffrange = [0, 500]
    window = ndis_std(data[:, 0], mean, sigma)
    newdata = deepcopy(data)
    newdata[:, 1] = newdata[:, 1] * window
    acdat, acpeaks = autocorr(newdata)
    acdat = limit_data(acdat, diffrange[0], diffrange[1])
    return acpeaks, acdat


def win_autocorr_grid(rawdata, binsize, wbin, window_fwhm, diffrange):
    # Prepare data
    mindat = np.amin(rawdata[:, 0])
    maxdat = np.amax(rawdata[:, 0])

    mzdata = linearize(rawdata, binsize, 3)
    mzdata = pad_two_power(mzdata)

    xvals = np.arange(mindat, maxdat, wbin)

    results = np.array([windowed_autocorr(mzdata, x, x * window_fwhm, diffrange=diffrange)[1] for x in xvals])

    intdat = results[:, :, 1]
    yvals = np.unique(results[:, :, 0])
    xgrid, ygrid = np.meshgrid(xvals, yvals, indexing="ij")
    out = np.transpose([np.ravel(xgrid), np.ravel(ygrid), np.ravel(intdat)])
    return out


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
        print(slope, cierr, rsquared, slope_std_error * np.sqrt(n), pvalue)
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
        pass

    # Sum and combine the spectra into a global master
    alignedsum = np.average(aligned[:, :, 1], axis=0)
    combined = np.transpose([aligned[0, :, 0], alignedsum])

    return combined, aligned


def calc_FWHM(peak, data):
    index = nearest(data[:, 0], peak)
    int = data[index, 1]
    leftwidth = 0
    rightwidth = 0
    counter = 1
    leftfound = False
    rightfound = False
    while rightfound is False or leftfound is False:
        if leftfound is False:
            if index - counter < 0:
                leftfound = True
            elif data[index - counter, 1] <= (int) / 2.:
                leftfound = True
                leftwidth += 1
            else:
                leftwidth += 1
        if rightfound is False:
            if index + counter >= len(data):
                rightfound = True
            elif data[index + counter, 1] <= (int) / 2.:
                rightfound = True
                rightwidth += 1
            else:
                rightwidth += 1
        counter += 1

    indexstart = index - leftwidth
    indexend = index + rightwidth
    if indexstart < 0:
        indexstart = 0
    if indexend >= len(data):
        indexend = len(data) - 1

    FWHM = data[indexend, 0] - data[indexstart, 0]
    return FWHM, [data[indexstart, 0], data[indexend, 0]]


def peaks_error_FWHM(pks, data):
    """
    Calculates the error of each peak in pks using FWHM.
    Looks for the left and right point of the peak that is 1/2 the peaks max intensity, rightmass - leftmass = error
    :param pks:
    :param data: self.data.massdat
    :return:
    """
    pmax = np.amax([p.height for p in pks.peaks])
    try:
        datamax = np.amax(np.asarray(data)[:, 1])
    except:
        datamax = 0
    div = datamax / pmax
    for p in pks.peaks:
        int = p.height
        index = nearest(data[:, 0], p.mass)
        leftwidth = 0
        rightwidth = 0
        counter = 1
        leftfound = False
        rightfound = False
        val = (int * div) / 2.
        while rightfound is False or leftfound is False:
            if leftfound is False and index - counter >= 0:
                if data[index - counter, 1] <= val:
                    leftfound = True
                    leftwidth += 1
                else:
                    leftwidth += 1
            else:
                leftfound = True
            if rightfound is False and index + counter < len(data):
                if data[index + counter, 1] <= val:
                    rightfound = True
                    rightwidth += 1
                else:
                    rightwidth += 1
            else:
                rightfound = True
            counter += 1

        indexstart = index - leftwidth
        indexend = index + rightwidth
        if indexstart < 0:
            indexstart = 0
        if indexend >= len(data):
            indexend = len(data) - 1

        mlow = data[indexstart, 0]
        mhigh = data[indexend, 0]

        p.errorFWHM = mhigh - mlow
        p.intervalFWHM = [mlow, mhigh]

        diff = np.abs(np.array(p.intervalFWHM) - p.mass)
        r = safedivide1(diff, p.errorFWHM)
        # Check to make sure that the FWHM is symmetric. If not, flag it and bring things back.
        cutoff = 0.75
        if r[0] > cutoff:
            mlow = p.mass - diff[1] * cutoff / (1 - cutoff)
            p.errorFWHM = mhigh - mlow
            p.intervalFWHM = [mlow, mhigh]
            p.badFWHM = True
            pass
        elif r[1] > cutoff:
            mhigh = p.mass + diff[0] * cutoff / (1 - cutoff)
            p.errorFWHM = mhigh - mlow
            p.intervalFWHM = [mlow, mhigh]
            p.badFWHM = True
            pass
        else:
            p.badFWHM = False

        start = p.intervalFWHM[0]
        end = p.intervalFWHM[1]
        p.centroid = center_of_mass(data, start, end)[0]
        # print("Apex:", p.mass, "Centroid:", p.centroid, "FWHM Range:", p.intervalFWHM)


def peaks_error_mean(pks, data, ztab, massdat, config):
    """
    Calculates error using the masses at different charge states.
    For each peak, finds the local max of the peak at each charge state, and does a weighted mean and weighted std. dev.
    :param pks:
    :param data: self.data.massgrid
    :param ztab: self.data.ztab
    :param massdat: self.data.massdat
    :param config: self.config
    :return:
    """
    # Reshape the data
    length = len(data) / len(ztab)
    data2 = data.reshape(int(length), len(ztab))
    # Loop over each peak
    for pk in pks.peaks:
        # Set the window as 2 x FWHM if possible
        window = config.peakwindow
        if pk.errorFWHM > 0:
            window = 2 * pk.errorFWHM

        # Grab a slice of the data around the peak mass
        index = nearest(massdat[:, 0], pk.mass)
        startindmass = nearest(massdat[:, 0], massdat[index, 0] - window)
        endindmass = nearest(massdat[:, 0], massdat[index, 0] + window)
        data3 = data2[startindmass:endindmass]

        # Get the peak mass and intensity at each charge state
        masses = []
        ints = []
        for z in range(0, len(ztab)):
            tmparr = data3[:, z]
            ints.append(np.amax(tmparr))
            masses.append(massdat[startindmass + np.argmax(tmparr), 0])

        # Convert to Arrays
        ints = np.array(ints)
        masses = np.array(masses)

        # Set a 1% Threshold
        b1 = ints > 0.01 * np.amax(ints)
        # Calculate weighted standard deviation
        std = weighted_std(masses[b1], ints[b1])

        # Set the parameters
        pk.errormean = std
        pk.errormasses = masses
        pk.errorints = ints
    pass


def subtract_and_divide(pks, basemass, refguess):
    avgmass = []
    ints = []
    for p in pks.peaks:
        if p.ignore == 0:
            mass = p.mass - basemass
            num = np.round(mass / float(refguess))
            if num != 0:
                avg = mass / num
                avgmass.append(avg)
                ints.append(p.height)

    return np.average(avgmass, weights=ints)


if __name__ == "__main__":
    testfile = "C:\Python\\UniDec\TestSpectra\\test_imms.raw"
    # data = waters_convert(testfile)
    # print(np.amax(data))
    data = waters_convert2(testfile)
    print(np.amax(data))

    exit()

    x = [0., 1., 2., 3., 4.]
    y = [1, 0.7, 0.5, 0.4, 0.3]
    import matplotlib.pyplot as plt

    y = logistic(np.array(x), 2, -10, 1, 1)
    plt.plot(x, y)
    plt.show()
    # fit=sig_fit(x,y)
    # print fit
    pass
