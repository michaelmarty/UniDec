import os
import platform
import sys
import math
import subprocess
import time
import re
import decimal
from bisect import bisect_left
from copy import deepcopy
import zipfile
# noinspection PyUnresolvedReferences
import numpy as np
import scipy.fft
import scipy.ndimage.filters as filt
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy import signal
from scipy import fftpack
import matplotlib.cm as cm
import matplotlib.colors as colors
from unidec.modules.fitting import *
from itertools import cycle
#from numba import jit

try:
    from unidec.modules.mzMLimporter import mzMLimporter
except:
    pass
from unidec.modules import unidecstructure
import fnmatch
from unidec.modules.mzXML_importer import mzXMLimporter

try:
    from unidec.modules import data_reader
except:
    print("Could not import data reader: unidectools")

# from modules.waters_importer.WatersImporter import WatersDataImporter as WDI
try:
    from unidec.modules.waters_importer.WatersImporter import WatersDataImporter as WDI
except:
    print("Could not import Waters Data Importer")

# from modules.thermo_reader.ThermoImporter import ThermoDataImporter
try:
    from unidec.modules.thermo_reader.ThermoImporter import ThermoDataImporter
except:
    print("Could not import Thermo Data Importer")

is_64bits = sys.maxsize > 2 ** 32
protmass = 1.007276467
oxmass = 15.994914

known_extensions = [".raw", ".d", ".mzml", ".gz", ".mzml.gz", ".mzxml", ".mzML", ".mzML.gz", ".mzXML",
                    ".txt", ".csv", ".dat", ".npz"]

HTseqDict = {'2': '101', '3': '1110100', '4': '000100110101111', '5': '0000100101100111110001101110101',
             '6': '000001000011000101001111010001110010010110111011001101010111111',
             '7': '0000001000001100001010001111001000101100111010100111110100001110001001001101101011011110110001101001011101110011001010101111111',
             '8': '000000010111000111011110001011001101100001111001110000101011111111001011110100101000011011101101111101011101000001100101010100011010110001100000100101101101010011010011111101110011001111011001000010000001110010010011000100111010101101000100010100100011111',
             '-2': '010', '-3': '0001011', '-4': '111011001010000', '-5': '1111011010011000001110010001010',
             '-103': '0010111', '-105': '1010111011000111110011010010000', '5S4': '1001011001111100011011101010000',
             '5S13': '1111100011011101010000100101100', '5S15': '1110001101110101000010010110011'
             }


def get_importer(path):
    extension = os.path.splitext(path)[1]
    extension = extension.lower()
    if os.path.isfile(path) or os.path.isdir(path):
        if extension == ".mzml" or extension == ".gz":
            # mzML file
            d = mzMLimporter(path, gzmode=True)
        elif extension == ".mzxml":
            d = mzXMLimporter(path)
        elif extension == ".raw" and not os.path.isdir(path):
            # Thermo Raw File
            try:
                d = ThermoDataImporter(path)
            except:
                d = data_reader.DataImporter(path)
        elif extension == ".raw" and os.path.isdir(path):
            # Waters Raw Directory
            d = WDI(path, do_import=False)
        elif extension == ".txt" or extension == ".dat" or extension == ".csv":
            # print("Text Files Not Supported for This Operation")
            return None
        else:
            # Some other file type
            d = data_reader.DataImporter(path)
    else:
        print("Get Importer Failed. File Not Found:", path)
        return None

    return d


def get_max_time(path):
    if os.path.splitext(path)[1] == ".txt":
        return 0

    importer = get_importer(path)
    maxtime = importer.get_max_time()
    return maxtime


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


luminance_cutoff = 135


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


def create_color_cycle(seq='bgrcmk'):
    cycol = cycle(seq)
    return cycol


def get_color_from_index(index, seq="custom1"):
    if seq == "tab":
        carray = colors.TABLEAU_COLORS.keys()
    elif seq == "custom1":
        carray = ["purple", "blue", "dodgerblue", "cyan", "green", "lime", "gold", "orange", "coral", "red", "magenta"]
    else:
        carray = cm.get_cmap(seq).colors
    color = list(carray)[index % len(carray)]
    return color




def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)


def match_dirs_recursive(topdir, ending=".raw"):
    found_dirs = []
    for root, dirs, files in os.walk(topdir):
        for d in dirs:
            if d.endswith(ending):
                found_dirs.append(os.path.join(root, d))
    return np.array(found_dirs)


def match_files_recursive(topdir, ending=".raw"):
    found_files = []
    for root, dirs, files in os.walk(topdir):
        for f in files:
            if f.endswith(ending):
                found_files.append(os.path.join(root, f))
    return np.array(found_files)


def find_kernel_file(kernel_path):
    kernel_name = os.path.basename(kernel_path)

    if kernel_name.split('_')[-1] != "mass.txt":
        # If the filename doesn't end in "_mass.txt", it's not a mass file
        is_mass = False
    else:
        # Otherwise, it might be a mass file. Check if linear
        kernel_dat = np.loadtxt(kernel_path)
        kernel_diff = np.diff(kernel_dat[:, 0])
        is_mass = np.all(kernel_diff == kernel_diff[0])

    if is_mass:
        return kernel_path
    else:
        # If not a mass file, find the mass file
        bare_name = os.path.splitext(kernel_name)[0]
        kernel_path2 = os.path.dirname(kernel_path) + "\\" + bare_name + "_unidecfiles\\" + \
                       bare_name + "_mass.txt"

        try:
            with open(kernel_path2, "r") as f:
                return kernel_path2
        except (IOError, FileNotFoundError) as err:
            # print("Please deconvolve the m/z file [" + kernel_name + "] with UniDec first.")
            return None


def isempty(thing):
    """
    Checks wether the thing, a list or array, is empty
    :param thing: Object to check
    :return: Boolean, True if not empty
    """
    try:
        if np.asarray(thing, dtype=object).size == 0 or thing is None:
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


def strip_char_from_string(s, char):
    """
    Strip a character from a string
    :param s: String
    :param char: Character to strip
    :return: String without character
    """
    return s.replace(char, "")


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
        return a * 0


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
    :param data: data array
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
    In a sorted array, quickly find the position of the element closest to the target.
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
        # Center of Mass
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(data, start, end)
        else:
            val, junk = center_of_mass(data, data[0, 0], data[len(data) - 1, 0])
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 4:
        # Local Max Position
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
        # Center of Mass for intensity squared
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(data, start, end, power=2)
        else:
            val, junk = center_of_mass(data, data[0, 0], data[len(data) - 1, 0], power=2)
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 8:
        # Center of Mass for intensity cubed
        if window is not None:
            start = x - window
            end = x + window
            val, junk = center_of_mass(data, start, end, power=3)
        else:
            val, junk = center_of_mass(data, data[0, 0], data[len(data) - 1, 0], power=3)
            print("No window set for center of mass!\nUsing entire data range....")

    elif extract_method == 9:
        # Center of mass Remove data points that fall below 50% threshold squared
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
        # Center of Mass Remove data points that fall below 50% threshold cubed
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

    elif extract_method == 12:
        # Peak Area with 10% relative threshold
        start = nearest(data[:, 0], (x - window))  # x values should be sorted
        end = nearest(data[:, 0], (x + window))
        data_slice = data[start:end]
        max_index = np.argmax(data_slice[:, 1])
        local_height = data_slice[max_index, 1]

        boo2 = data[:, 1] > local_height * 0.1
        cutdat = data[boo2]
        if window is not None:
            start = x - window
            end = x + window
            val, junk = integrate(cutdat, start, end)
        else:
            index = nearest(data[:, 0], x)
            val = data[index, 1]
            print("NEED TO SET INTEGRAL WINDOW!\nUsing Peak Height Instead")

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


def extract_from_data_matrix(xaxis, matrix, midpoint, extract_method=0, window=0):
    if extract_method == 0:
        index = nearest(xaxis, midpoint)
        data = matrix[:, index]
    elif extract_method == 1:
        startindex = nearest(xaxis, midpoint - window)
        endindex = nearest(xaxis, midpoint + window)
        data = np.amax(matrix[:, startindex:endindex], axis=1)
    elif extract_method == 2:
        startindex = nearest(xaxis, midpoint - window)
        endindex = nearest(xaxis, midpoint + window)
        data = np.sum(matrix[:, startindex:endindex], axis=1)
    else:
        print("Grid Extraction Method Not Recognized")
        data = []
    return data


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


def kendrick_analysis(massdat, kendrickmass, centermode=1, nbins=50, transformmode=1, xaxistype=1, massrange=None):
    # Calculate Defects for Deconvolved Masses
    if kendrickmass == 0:
        print("Error: Kendrick mass is 0.")
        return None, None, None, None, None
    if massrange is not None:
        massdat = datachop(massdat, massrange[0], massrange[1])
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


def header_test(path, deletechars=None, delimiter=" |\t", strip_end_space=True):
    """
    A quick test of a file to remove any problems from the header.

    If scans through each line of an input file specificed by path and count the number of lines that cannot be
    completely converted to floats. It breaks the loop on the first hit that is entirely floats, so any header lines
    blow that will not be counted.
    :param path: Path of file to be scanned
    :param deletechars: Characters to be removed from the file before scanning
    :param delimiter: Delimiter to be used in the file
    :param strip_end_space: Boolean to strip the end space from the line
    :return: Integer length of the number of header lines
    """
    header = 0
    try:
        with open(path, "r") as f:
            for line in f:
                # If the line is empty, skip it and add to the header
                if line == "\n":
                    header += 1
                    continue
                # Remove any characters that are defined in deletechars
                if deletechars is not None:
                    for c in deletechars:
                        line = line.replace(c, "")

                # Strip the end space if strip_end_space is True
                if strip_end_space:
                    if line[-1] == "\n" and len(line) > 1:
                        line = line[:-1]
                    if line[-1] == " " and len(line) > 1:
                        line = line[:-1]

                # Split the line by the delimiter and try to convert each element to a float
                for sin in re.split(delimiter, line):  # line.split(delimiter):
                    try:
                        float(sin)
                    except ValueError:
                        # print(sin, line)
                        header += 1
                        break
        # If the header is greater than 0, print the header length
        if header > 0:
            print("Header Length:", header)
    except (ImportError, OSError, AttributeError, IOError) as e:
        print("Failed header test", e)
        header = 0
    return int(header)


'''
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
    return data'''


def waters_convert2(path, config=None, outfile=None, time_range=None):
    data = WDI(path).get_data(time_range=time_range)

    if outfile is None:
        outfile = os.path.join(path, "converted_rawdata.txt")
    np.savetxt(outfile, data)
    return data


def get_polarity(path, importer=None):
    if importer is None:
        d = get_importer(path)
    else:
        d = importer
    if d is not None:
        polarity = d.get_polarity()
        return polarity
    else:
        return None


def load_mz_file(path, config=None, time_range=None, imflag=0):
    """
    Loads a text or mzml file
    :param path: File path to load
    :param config: UniDecConfig object
    :return: Data array
    """
    if config is None:
        extension = os.path.splitext(path)[1]
        extension = extension.lower()
    else:
        extension = config.extension.lower()

    if not os.path.isfile(path) and not os.path.isdir(path):
        print("Load MZ File Failed. File Not Found:", path)
        raise IOError

    if not os.path.isfile(path):
        if os.path.isdir(path) and extension == ".raw":
            try:
                outfile = os.path.splitext(path)[0] + "_rawdata.txt"
                print("Trying to convert Waters File:", outfile)
                data = waters_convert2(path, config, outfile=outfile, time_range=time_range)
            except Exception as e:
                print("Attempted to convert Waters Raw file but failed", e, path)
                raise IOError
        elif os.path.isdir(path) and extension == ".d":
            try:
                print("Trying to convert Agilent File:", path)
                data = data_reader.DataImporter(path).get_data(time_range=time_range)
                txtname = path[:-2] + ".txt"
                np.savetxt(txtname, data)
                print("Saved to:", txtname)
            except Exception as e:
                print("Attempted to convert Agilent Raw file but failed", e, path)
                raise IOError
        else:
            print("Attempted to open:", path)
            print("\t but failed. File type may not be supported...")
            raise IOError
    else:
        if extension == ".txt" or extension == ".dat":
            try:
                data = np.loadtxt(path, skiprows=header_test(path))
            except:
                c = lambda s: float(str(s.decode("UTF-8")).replace(",", ""))  # .replace("'", "").replace("b", "")
                data = np.genfromtxt(path, skip_header=header_test(path, deletechars=","), converters={0: c, 1: c})
            # data = np.loadtxt(path, skiprows=header_test(path, delimiter=","), delimiter=",")
        elif extension == ".csv":
            try:
                data = np.loadtxt(path, delimiter=",", skiprows=1, usecols=(0, 1))
            except:
                data = np.loadtxt(path, delimiter=",", skiprows=8, usecols=(0, 1))
        elif extension == ".mzml" or extension == ".gz":
            if extension == ".mzml":
                num = -5
            if extension == '.gz':
                num = -8
            txtname = path[:num] + ".txt"
            if imflag == 1:
                if config.compressflag == 1:
                    mzbins = config.mzbins
                else:
                    mzbins = None
                data = mzMLimporter(path).get_im_data(time_range=time_range, mzbins=mzbins)

                np.savetxt(txtname, sparse(data))
                print("Saved to:", txtname)
            else:
                try:
                    data = mzMLimporter(path).get_data(time_range=time_range)
                except Exception as e:
                    print("Error with loading mzML file:", e)
                    data = mzMLimporter(path, nogz=True).get_data(time_range=time_range)
                np.savetxt(txtname, data)
                print("Saved to:", txtname)
        elif extension.lower() == ".raw":
            data = ThermoDataImporter(path).get_data(time_range=time_range)
            txtname = path[:-4] + ".txt"
            np.savetxt(txtname, data)
            print("Saved to:", txtname)
        elif extension.lower() == ".mzxml":
            data = mzXMLimporter(path).get_data(time_range=time_range)
            txtname = path[:-6] + ".txt"
            np.savetxt(txtname, data)
            print("Saved to:", txtname)
        elif extension.lower() == ".wiff":
            print("Trying to convert Sciex File (no promises...):", path)
            data = data_reader.DataImporter(path).get_data(time_range=time_range)
            txtname = path[:-4] + ".txt"
            np.savetxt(txtname, data)
            print("Saved to:", txtname)
        elif extension.lower() == ".npz":
            data = np.load(path, allow_pickle=True)['data']
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


def dataexport(datatop, fname, header=""):
    try:
        np.savetxt(fname, datatop, fmt='%f', header=header)
    except:
        path = os.path.join(os.getcwd(), fname)
        path = "\\\\?\\%s" % path
        np.savetxt(path, datatop, fmt='%f', header=header)
        print("NOTE: Your path length might exceed the limit for Windows. Please shorten your file name.")
    pass


def dataexportbin(datatop, fname):
    try:
        datatop.astype(np.single).tofile(fname, sep="")
    except:
        path = os.path.join(os.getcwd(), fname)
        path = "\\\\?\\%s" % path
        datatop.astype(np.single).tofile(path, sep="")
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


def mergedata2d(x1, y1, x2, y2, z2, method="linear"):
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
    zout = griddata(np.transpose([oldx, oldy]), np.ravel(z2), (newx, newy), method=method, fill_value=0)
    return zout


# ........................................
#
# Data Processing
#
# ..........................................................


def auto_peak_width(datatop, psfun=None, singlepeak=False):
    maxpos = np.argmax(datatop[:, 1])
    maxval = datatop[maxpos, 0]

    if maxval == 0:
        return 1, 0, 0

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


def remove_noise(datatop, percent=None, silent=True):
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
    if not silent:
        print("Removed", l1 - len(datatop), "points.")
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


def linear_axis(raw_data):
    # Get the sparse data
    sparse_data = np.unique(raw_data)
    # Find the sample rate
    sample_rate = np.amin(np.diff(sparse_data))
    # Setup the new axis
    min = np.amin(sparse_data)
    max = np.amax(sparse_data)
    newaxis = np.arange(min, max + sample_rate, sample_rate)
    # Round the old axis to match the new. Otherwise, small differences will be deadly
    oldraw = np.round((raw_data - min) / sample_rate) * sample_rate + min
    oldaxis = np.unique(oldraw)
    return newaxis, oldaxis, oldraw


def mergedata2dexact(x1, y1, x2, y2, z2):
    p1 = np.vectorize(complex)(x1, y1)
    p2 = np.vectorize(complex)(x2, y2)
    boo3 = np.isin(p1, p2)

    z1 = np.zeros_like(x1)
    z1[boo3] = z2
    return z1


def unsparse(rawdata):
    # Create new axes
    mzaxis2, mzaxis, rawdata[:, 0] = linear_axis(rawdata[:, 0])
    dtaxis2, dtaxis, rawdata[:, 1] = linear_axis(rawdata[:, 1])

    # Create new grids
    mzaxis2grid, dtaxis2grid = np.meshgrid(mzaxis2, dtaxis2, sparse=False, indexing='ij')

    # Unsparse
    ints2 = mergedata2dexact(np.ravel(mzaxis2grid), np.ravel(dtaxis2grid), rawdata[:, 0], rawdata[:, 1], rawdata[:, 2])

    # Wrap everything together
    rawdata3 = np.transpose([np.ravel(mzaxis2grid), np.ravel(dtaxis2grid), ints2])

    # Create 1D sum
    intgrid = ints2.reshape((len(mzaxis2), len(dtaxis2)))
    rawdata = np.transpose([mzaxis2, np.sum(intgrid, axis=1)])

    return rawdata3, rawdata


def sparse(rawdata):
    boo1 = rawdata[:, 2] != 0
    return rawdata[boo1]


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


def lintegrate(datatop, intx, fastmode=False):
    """
    Linearize x-axis by integration.

    Each intensity value in the old data gets proportionally added into the new x-axis.

    The total sum of the intensity values should be constant.
    :param datatop: Data array
    :param intx: New x-axis for data
    :param fastmode: If True, uses a faster but less accurate method of integration that just picks the nearest point.
    :return: Integration of intensity from original data onto the new x-axis.
    Same shape as the old data but new length.
    """
    length = len(datatop)
    l2 = len(intx)
    inty = np.zeros_like(intx)
    for i in range(0, length):
        x = datatop[i, 0]
        y = datatop[i, 1]
        if intx[0] < x < intx[len(intx) - 1]:
            index = nearest(intx, x)
            if fastmode:
                inty[index] += y
            else:
                if intx[index] == x:
                    inty[index] += y
                elif intx[index] < x and index < l2 - 1:
                    index2 = index + 1
                    interpos = linear_interpolation(intx[index], intx[index2], x)
                    inty[index] += (1 - interpos) * y
                    inty[index2] += interpos * y
                elif intx[index] > x and index > 0:
                    index2 = index - 1
                    interpos = linear_interpolation(intx[index], intx[index2], x)
                    inty[index] += (1 - interpos) * y
                    inty[index2] += interpos * y
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
    f = interp1d(datatop[:, 0], datatop[:, 1], fill_value=0, bounds_error=False)
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
    config = unidecstructure.UniDecConfig()
    floatdata = np.array(datatop[:, 0], dtype=config.dtype)
    testunique = np.unique(floatdata)
    if len(testunique) != len(floatdata):
        print("Removing Duplicates")
        num, start = np.histogram(floatdata, bins=testunique)
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
    try:
        maxval = np.amax(datatop[:, 1])
        datatop[:, 1] = datatop[:, 1] / maxval
    except Exception as e:
        pass
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


def fake_log(data, percent=50):
    nonzerodata = data[data > 0]
    l1 = len(nonzerodata)
    sdat = np.sort(np.ravel(nonzerodata))
    index = round(l1 * percent / 100.)
    non_zero_min = sdat[index]
    newdata = np.log10(np.clip(data, non_zero_min, np.amax(data)))
    newdata -= np.amin(newdata)
    return newdata


def remove_middle_zeros(data):
    boo1 = data[1:len(data) - 1, 1] != 0
    boo2 = data[:len(data) - 2, 1] != 0
    boo3 = data[2:len(data), 1] != 0
    boo4 = np.logical_or(boo1, boo2)
    boo5 = np.logical_or(boo3, boo4)
    boo6 = np.concatenate(([True], boo5, [True]))
    return data[boo6]


def smash(data, midpoint, window):
    # Set all the data from midpoint-window to midpoint+window to zero
    b1 = data[:, 0] > midpoint - window
    b2 = data[:, 0] < midpoint + window
    b3 = np.logical_and(b1, b2)
    data[b3, 1] = 0
    return data


def smash2d(data, midpoint, window, ymidpoint, ywindow):
    # Set all the data from midpoint-window to midpoint+window to zero
    b1 = data[:, 0] > midpoint - window
    b2 = data[:, 0] < midpoint + window
    b3 = np.logical_and(b1, b2)

    b1 = data[:, 1] > ymidpoint - ywindow
    b2 = data[:, 1] < ymidpoint + ywindow
    b4 = np.logical_and(b1, b2)

    b5 = np.logical_and(b3, b4)

    data[b5, 2] = 0
    return data


def dataprep(datatop, config, peaks=True, intthresh=True, silent=False):
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
    :param peaks: Boolean to perform peak detection
    :param intthresh: Boolean to perform intensity thresholding
    :param silent: Boolean to suppress output
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
    if type(newmin) != str and type(newmax) != str:
        # Crop Data
        data2 = datachop(deepcopy(datatop), newmin, newmax)
    else:
        data2 = deepcopy(datatop)
    print("Test")
    if config.smashflag == 1:
        print("Smashing!:", config.smashlist)
        for i in range(0, len(config.smashlist)):
            data2 = smash(data2, config.smashlist[i][0], config.smashlist[i][1])

    if len(data2) == 0:
        print("Error: m/z range is too small. No data fits the range.")
    # correct for detector efficiency
    if va != 0:
        # data2=detectoreff(data2,9.1)
        data2 = detectoreff(data2, va)

    # Smooth Data
    if smooth > 0:
        data2 = gsmooth(data2, smooth)

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

    # Data Reduction
    if redper > 0:
        data2 = remove_noise(data2, redper)

    if config.datanorm == 1:
        # Normalization
        data2 = normalize(data2)

    # Intensity Threshold
    if thresh > 0 and intthresh:
        # print(len(data2))
        data2 = intensitythresh(data2, thresh)  # thresh
        if not silent:
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
        if not silent:
            print(data2)

    if linflag == 2:
        try:
            data2 = remove_middle_zeros(data2)
        except:
            pass
        pass

    # Remove Duplicate Data Points
    data2 = removeduplicates(data2)

    return data2


# ......................................................
#
# unidec Functions
#
# ............................................................................


def exe_call(call, silent=False):
    """
    Run the  binary specified in the call.
    If silent is False (default), the output from exepath will be printed to the standard out.
    If silent is True, the output is suppressed.

    :param call: Call arguments to be passed ot the shell
    :param silent: Whether to print the output of exepath to the standard out
    :return: Standard error of exepath execution
    """
    print("System Call:", call)
    result = subprocess.run(call, shell=False, capture_output=True, text=True)
    out = result.returncode
    if not silent:
        print(result.stdout)
    return out


def unidec_call(config, silent=False, conv=False, **kwargs):
    """
    Run the unidec binary specified by exepath with the configuration file specified by configfile.
    If silent is False (default), the output from exepath will be printed to the standard out.
    If silent is True, the output is suppressed.

    The binary is run as a commandline subprocess with a call of (for example):

    unidec.exe conf.dat

    :param config: Config object
    :param silent: Whether to print the output of exepath to the standard out
    :param conv: Whether to call the convolution function only rather than the standard deconvolution
    :param kwargs:
    :return: Standard error of exepath execution
    """

    call = [config.UniDecPath, str(config.confname)]

    if config.autotune:
        call.append("-autotune")

    if conv:
        call.append("-conv")

    out = exe_call(call, silent=silent)
    return out

#@jit(nopython=True)
def peakdetect(data, config=None, window=10, threshold=0, ppm=None, norm=True):
    """
    Simple peak detection algorithm.

    Detects a peak if a given data point is a local maximum within plus or minus config.peakwindow.
    Peaks must also be above a threshold of config.peakthresh * max_data_intensity.

    The mass and intensity of peaks meeting these criteria are output as a P x 2 array.

    :param data: Mass data array (N x 2) (mass intensity)
    :param config: UniDecConfig object
    :param window: Tolerance window of the x values
    :param threshold: Threshold of the y values
    :param ppm: Tolerance window in ppm
    :param norm: Whether to normalize the data before peak detection
    :return: Array of peaks positions and intensities (P x 2) (mass intensity)
    """
    if config is not None:
        window = config.peakwindow / config.massbins
        threshold = config.peakthresh
        norm = config.normthresh

    peaks = []
    length = len(data)
    shape = np.shape(data)
    if length == 0 or shape[1] != 2:
        return np.array(peaks)

    if norm:
        maxval = np.amax(data[:, 1])
    else:
        maxval = 1
    for i in range(0, length):
        if data[i, 1] > maxval * threshold:
            if ppm is not None:
                ptmass = data[i, 0]
                newwin = ppm * 1e-6 * ptmass
                start = nearest(data[:, 0], ptmass - newwin)
                end = nearest(data[:, 0], ptmass + newwin)
            else:
                start = i - window
                end = i + window

                start = int(start)
                end = int(end) + 1

                if start < 0:
                    start = 0
                if end > length:
                    end = length

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
    stick[np.array(mztab[:, 2]).astype(int)] = mztab[:, 1]
    bool1 = [np.abs(xvals - xvals[i]) < window for i in range(0, xlen)]
    kernels = np.array([make_peak_shape(-xvals[bool1[i]], psfun, fwhm, -xvals[i]) for i in range(0, xlen)],
                       dtype=object)
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
    temp[np.array(mztab[:, 2]).astype(int)] = mztab[:, 1]
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
    temp[np.array(mztab[:, 2]).astype(int)] = mztab[:, 1]
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
    cmap = colors.LinearSegmentedColormap("newcmap", cdict)
    return cmap


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
    startindex = startindex.astype(int)
    basemass = array2[:, 0]
    basemass = basemass.astype(float)
    omass = array2[:, 1]
    omass = omass.astype(float)
    finlist = []
    for index in np.ndindex(tup):
        total = np.sum((index + startindex) * omass + basemass)
        if total > 0:
            finlist.append(total)
    return finlist


def index_to_oname(index, startindex, names):
    # print(index, startindex, names)
    name = ""
    for i in range(0, len(index)):
        val = index[i] + startindex[i]
        if val > 0:
            if names[i] == "":
                name = name + str(val)
            else:
                name = name + str(val) + "[" + names[i] + "] "
    return name


def combine_all(array2):
    lens = lengths(array2)
    tup = tuple(lens)
    print("Starting combining all: ", np.prod(lens))
    startindex = array2[:, 2]
    startindex = startindex.astype(int)
    basemass = array2[:, 0]
    basemass = basemass.astype(float)
    omass = array2[:, 1]
    omass = omass.astype(float)
    # names = array2[:, 4]

    # namelist = np.array([index for index in np.ndindex(tup)])
    namelist = np.array(list(np.ndindex(tup)), dtype=np.dtype((int, len(lens))))
    # finlist = np.array([np.sum((index + startindex) * omass + basemass) for index in namelist])
    finlist = np.sum((namelist + startindex) * omass + basemass, axis=1)

    """
    for index in np.ndindex(tup):
        vals = index + startindex
        b1 = vals > 0
        vstr = vals[b1].astype(str)
        nstr = "[" + np.char.array(names[b1]) + "] "
        nstr[nstr == "[] "] = ""
        fstr = vstr + nstr
        name = ""
        for i in range(0, len(vstr)):
            name += fstr[i]
        namelist.append(name)"""
    '''
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
            # print index,name'''

    b1 = finlist != 0
    finlist = finlist[b1]
    namelist = namelist[b1]
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
                index = np.zeros(len(oligos)).astype(int)
                index[i] = j
                oligonames.append(index - start)
                '''
                if j > 0 or oligos[i][4] == "":
                    if oligos[i][4] == "":
                        oligonames.append(str(j))
                    else:
                        oligonames.append(str(j) + "[" + oligos[i][4] + "]")
                else:
                    oligonames.append("")
                    # self.oligonames.append(str(j)+""+oligos[i][4])'''
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


def get_glyco_indexes(oligomerlist, printoutput=False):
    names = oligomerlist[:, 4]
    sname = ''
    hname = ''
    gname = ''
    fname = ''
    hindex = -1
    gindex = -1
    sindex = -1
    findex = -1
    for i, n in enumerate(names):
        if "SA" in n:
            sname = n
            sindex = i
        if "Hex" in n:
            hname = n
            hindex = i
        if "GlcNAc" in n:
            gname = n
            gindex = i
        if "Fuc" in n:
            fname = n
            findex = i
    if printoutput:
        print("Sialic Acid Name:", sname, sindex)
        print("Hexose Name:", hname, hindex)
        print("GlcNAc Name:", gname, gindex)
        print("Fucose Name:", fname, findex)
    return sindex, hindex, gindex, findex


def pair_glyco_matches(oligomasslist, oligonames, oligomerlist):
    oligomerlist = np.array(oligomerlist)
    startindex = oligomerlist[:, 2].astype(int)

    sindex, hindex, gindex, findex = get_glyco_indexes(oligomerlist)

    nmatches = len(oligonames)
    ns = np.zeros(nmatches)
    nh = np.zeros(nmatches)
    ng = np.zeros(nmatches)

    if sindex != -1:
        ns = oligonames[:, sindex] + startindex[sindex]

    if hindex != -1:
        nh = oligonames[:, hindex] + startindex[hindex]

    if gindex != -1:
        ng = oligonames[:, gindex] + startindex[gindex]
    # Number of Sialic Acid is less than or equal to number of hexoses
    b1 = ns <= nh
    # Number of Sialic Acids is also less than or equal to number of glcnacs
    b2 = ns <= ng
    b3 = np.logical_and(b1, b2)
    print(nmatches, len(oligomasslist[b3]))
    return oligomasslist[b3], oligonames[b3]


def match(pks, oligomasslist, oligonames, oligomerlist, tolerance=None, return_numbers=False):
    print("Starting Match")
    starttime = time.perf_counter()
    matches = []
    errors = []
    peaks = []
    names = []
    numbers = []

    oligomerlist = np.array(oligomerlist)
    startindex = oligomerlist[:, 2].astype(int)
    onames = oligomerlist[:, 4]

    for i in range(0, pks.plen):
        p = pks.peaks[i]
        target = p.mass
        nearpt = nearestunsorted(oligomasslist, target)
        match = oligomasslist[nearpt]
        error = target - match
        number = np.zeros(len(startindex))
        if tolerance is None or np.abs(error) < tolerance:
            name = index_to_oname(oligonames[nearpt], startindex, onames)
            if return_numbers:
                number = oligonames[nearpt] + startindex
        else:
            name = ""

        p.label = name
        p.match = match
        p.matcherror = error
        matches.append(match)
        errors.append(error)
        peaks.append(target)
        names.append(name)
        numbers.append(number)
    matchlist = [peaks, matches, errors, names]
    endtime = time.perf_counter()
    print("Matched in: ", endtime - starttime, "s")
    if return_numbers:
        return np.array(matchlist), np.array(numbers)
    else:
        return np.array(matchlist)


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
    xvals = np.arange(0, length).astype(float)
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


def autocorr(datatop, config=None, window=10):
    """
    Note: ASSUMES LINEARIZED DATA
    :param datatop: 1D data
    :param config: Config file (optional file)
    :param window: Window for peak detection, default 10 data points
    :return: Autocorr spectrum, peaks in autocorrelation.
    """
    corry = signal.fftconvolve(datatop[:, 1], datatop[:, 1][::-1], mode='same')
    if np.amax(corry) != 0:
        corry /= np.amax(corry)
        maxpos1 = np.argmax(datatop[:, 1])
        start = np.amax([maxpos1 - len(datatop) / 10, 0])
        end = np.amin([len(datatop) - 1, maxpos1 + len(datatop) / 10])
        cutdat = datatop[int(start):int(end)]
        if len(cutdat) < 20:
            cutdat = datatop
        # cutdat=datatop # Other old
        xdiff = np.mean(
            cutdat[1:, 0] - cutdat[:len(cutdat) - 1, 0])  # Less dangerous but still dangerous when non-linear
        # xdiff = datatop[1, 0] - datatop[0, 0] #OLD
        corrx = np.arange(0.0, len(corry)) * xdiff
        maxpos = np.argmax(corry)
        corrx = corrx - corrx[maxpos]
        autocorr = np.transpose([corrx, corry])
        boo1 = autocorr[:, 0] > xdiff
        cpeaks = peakdetect(autocorr[boo1], config, window=window)
        return autocorr, cpeaks

    else:
        return [[]], [[]]


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


def make_peak_shape(xaxis, psfun, fwhm, mid, norm_area=False, speedy=False):
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
    # If Speedy, limit the range to just within a few fwhm around the mid
    if speedy:
        # b1 = xaxis < mid+2*fwhm
        # b2 = xaxis > mid-2*fwhm
        # b1 = np.logical_and(b1, b2)
        b1 = np.abs(xaxis - mid) < 2 * fwhm
        fullx = xaxis
        xaxis = xaxis[b1]
        # print(len(xaxis), len(fullx))

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

    # Merge back in the speedy results
    if speedy:
        fullkernel = fullx * 0
        fullkernel[b1] = kernel
        kernel = fullkernel

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


def rfft(data, nmult=1):
    # X-axis
    massdif = data[1, 0] - data[0, 0]
    # fvals = np.fft.fftfreq(len(data), d=massdif)
    fvals = scipy.fft.rfftfreq(len(data) * nmult, d=massdif)
    # Y-axis
    # fft = np.fft.fft(data[:, 1])
    fft = scipy.fft.rfft(data[:, 1], n=len(data) * nmult)
    # fft = pyfftw.interfaces.numpy_fft.fft(data[:,1])
    # Cleanup
    aftdat = np.transpose([fvals, np.abs(fft)])
    aftdat[:, 1] /= np.amax(aftdat[:, 1])
    return aftdat


def fft_diff(data, diffrange=None):
    if diffrange is None:
        diffrange = [500., 1000.]
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
    Perform MacCoss method (https://pubs.acs.org/doi/pdf/10.1021/ac034790h) of getting peak intensities
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
            elif data[index - counter, 1] <= int / 2.:
                leftfound = True
                leftwidth += 1
            else:
                leftwidth += 1
        if rightfound is False:
            if index + counter >= len(data):
                rightfound = True
            elif data[index + counter, 1] <= int / 2.:
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


def peaks_error_FWHM(pks, data, level=0.5):
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
    try:
        div = datamax / pmax
    except:
        div = 1
    for p in pks.peaks:
        int = p.height
        index = nearest(data[:, 0], p.mass)
        leftwidth = 0
        rightwidth = 0
        counter = 1
        leftfound = False
        rightfound = False
        val = (int * div) * level
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
    pks.centroids = np.array([p.centroid for p in pks.peaks])
    pks.fwhms = np.array([p.errorFWHM for p in pks.peaks])


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


def subtract_and_divide(pks, basemass, refguess, outputall=False):
    avgmass = []
    ints = []
    nums = []
    masses = []
    for p in pks.peaks:
        if p.ignore == 0:
            mass = p.mass - basemass
            num = np.round(mass / float(refguess))
            if num != 0:
                avg = mass / num
                avgmass.append(avg)
                ints.append(p.height)
                nums.append(num)
                masses.append(p.mass)
                p.sdnum = num
                p.sdval = avg

    if outputall:
        return np.average(avgmass, weights=ints), avgmass, ints, nums, masses
    else:
        return np.average(avgmass, weights=ints)

def calc_swoop(sarray, adduct_mass=1):
    mz_mid, z_mid, z_spread, z_width = sarray
    z_mid = np.round(z_mid)
    z_width = np.round(z_width)
    mass = mz_mid * z_mid - adduct_mass * z_mid
    minz = z_mid - z_width
    maxz = z_mid + z_width
    z = np.arange(minz, maxz + 1)
    mz = (mass + adduct_mass * z) / z

    zup = np.round(z + z_spread)
    zdown = np.round(z - z_spread)

    return mz, z, zup, zdown

def get_swoop_mz_minmax(mz, i):
    if i == 0:
        mzmax = mz[i]
    else:
        mzmax = (mz[i] + mz[i - 1]) / 2
    if i == len(mz) - 1:
        mzmin = mz[i]
    else:
        mzmin = (mz[i] + mz[i + 1]) / 2
    return mzmin, mzmax

def within_ppm(theo, exp, ppmtol):
    return np.abs(((theo - exp) / theo) * 1e6) <= ppmtol



if __name__ == "__main__":
    testdir = "C:\\Data\\AgilentData"
    testfile = "LYZ-F319-2-11-22-P1-A1.D"
    # testfile = "2019_05_15_bsa_ccs_02.d"
    path = os.path.join(testdir, testfile)

    data = load_mz_file(path)

    exit()
    testfile = "C:\\Python\\UniDec3\\TestSpectra\\test_imms.raw"
    # data = waters_convert(testfile)
    # print(np.amax(data))
    data = waters_convert2(testfile)  # , time_range=[0,1])
    # print(np.amax(data))
    print(data)
    exit()
    testfile = "C:\\Python\\UniDec3\\TestSpectra\\test_imms.raw\\test_imms_imraw.txt"
    print(header_test(testfile))

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
