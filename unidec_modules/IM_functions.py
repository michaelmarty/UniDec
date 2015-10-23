import math

import numpy as np
import scipy.ndimage.filters as filt

# import time
from unidec_modules import unidectools as ud

__author__ = 'Michael.Marty'


def smooth_2d(igrid, smooth, smoothdt):
    """
    Gaussian smooth of 2D grid
    :param igrid: 2D numpy array
    :param smooth: width in x dimension (in data points)
    :param smoothdt: width in y dimension (in data points)
    :return: smoothed 2D numpy array
    """
    igrid = filt.gaussian_filter(igrid, [smooth, smoothdt])
    return igrid


def min_array(igrid, xwidth, ywidth):
    """
    Creates an array of local minimum values in a neighborhood of +/- width
    :param igrid: 2D numpy array
    :param xwidth: width in x direction
    :param ywidth: width in y direction
    :return: array of local minimum values in the same dimensions as igrid
    """
    mins = np.zeros_like(igrid)
    shape = igrid.shape
    for i in range(0, shape[0]):
        xstart = max([0, i - abs(xwidth)])
        xend = min([shape[0], i + abs(xwidth) + 1])
        for j in range(0, shape[1]):
            ystart = max([0, j - abs(ywidth)])
            yend = min([shape[1], j + abs(ywidth) + 1])
            mins[i, j] = np.amin(igrid[xstart:xend, ystart:yend])
    return mins


def intensitythresh(igrid, thresh):
    """
    Intensity threshold for 2D grid
    :param igrid: 2D numpy array
    :param thresh: Threshold
    :return: igrid with all elements below thresh set to 0
    """
    belowint = igrid < thresh
    igrid[belowint] = 0
    return igrid


def subtract_complex_2d(igrid, config):
    """
    Subtract a complex curved baseline.
    Creates a baseline using min_array, smooths it using smooth_2d, and subtracts it.
    :param igrid: 2D numpy array
    :param config: UniDecConfig object specifying width parameters in x and y dimensions
    :return: igrid with background subtraction applied
    """
    mins = min_array(igrid, config.subbuff, config.subbufdt)
    mins = smooth_2d(mins, config.subbuff * 2, config.subbufdt * 2)
    igrid = intensitythresh(igrid - mins, 0)
    return igrid


def compress_2d(xgrid, ygrid, zgrid, num):
    """
    Compresses a 2D grid in to smaller grids by merging n data points in the x dimension together as their average.
    :param xgrid: 2D numpy array
    :param ygrid: 2D numpy array
    :param zgrid: 2D numpy array
    :param num: number of consecutive data points in the x-dimension to average together
    :return: x, y, and z grids compressed along the x dimension to a new size
    """
    # Assumes all are in grid form
    x2 = np.array([np.average(xgrid[i:i + num], axis=0) for i in xrange(0, len(xgrid), num)])
    y2 = np.array([np.average(ygrid[i:i + num], axis=0) for i in xrange(0, len(ygrid), num)])
    z2 = np.array([np.average(zgrid[i:i + num], axis=0) for i in xrange(0, len(zgrid), num)])
    return x2, y2, z2


def linearize_2d(xvals, yvals, igrid, binsize):
    """
    Linearize the x-axis of a 2D grid
    :param xvals: 2D numpy array (x values)
    :param yvals: 2D numpy array (y values)
    :param igrid: 2D numpy array (intensity values)
    :param binsize: difference between consecutive linear x points
    :return: x, y, and z grid linearized along the x-axis
    """
    firstpoint = math.ceil(xvals[0] / binsize) * binsize
    lastpoint = math.floor(xvals[len(xvals) - 1] / binsize) * binsize
    intx = np.arange(firstpoint, lastpoint, binsize)
    iout = np.zeros((len(intx), len(yvals)))
    shape = igrid.shape
    for i in xrange(0, shape[0]):
        if intx[0] < xvals[i] < intx[len(intx) - 1]:
            index = ud.nearest(intx, xvals[i])
            # iout[index]+=C[i]
            if intx[index] == xvals[i]:
                iout[index] += igrid[i]
            if intx[index] < xvals[i] and index < shape[0] - 1:
                index2 = index + 1
                interpos = ud.linear_interpolation(intx[index], intx[index2], xvals[i])
                iout[index] += (1 - interpos) * igrid[i]
                iout[index2] += interpos * igrid[i]
            if intx[index] > xvals[i] and index > 0:
                index2 = index - 1
                interpos = ud.linear_interpolation(intx[index], intx[index2], xvals[i])
                iout[index] += (1 - interpos) * igrid[i]
                iout[index2] += interpos * igrid[i]

    xout, yout = np.meshgrid(intx, yvals, indexing='ij')
    return xout, yout, iout


def detectoreff_2d(igrid, xgrid, va):
    """
    Corrects for Q-TOF detector efficiency
    :param igrid: Intensity 2D numpy array
    :param xgrid: 2D numpy array of m/z values
    :param va: TOF acceleration voltage
    :return: igrid corrected for detector efficiency
    """
    eff = (1 - np.exp(-1620 * (va / xgrid) ** 1.75))
    igrid /= eff
    return igrid


def process_data_2d(xgrid, ygrid, igrid, config):
    """
    Process IM-MS raw data.

    1. Chop data to defined limits
    2. Linearize
    3. Smooth
    4. Subtract background
    5. Detector efficiency correction
    6. Normalize
    :param xgrid: 2D numpy array (x values)
    :param ygrid: 2D numpy array (y values)
    :param igrid: 2D numpy array (intensity values values)
    :param config: UniDecConfig object carrying parameters for processing
    :return: x, y, and z grids of processed data
    """
    # tstart = time.clock()
    if config.pusher != 0:
        ygrid = np.array(ygrid) * config.pusher * 0.001
    boo1 = xgrid > config.minmz
    boo2 = xgrid < config.maxmz
    boo3 = ygrid > config.mindt
    boo4 = ygrid < config.maxdt
    boofin = np.all([boo1, boo2, boo3, boo4], axis=0)
    xgrid = xgrid[boofin]
    ygrid = ygrid[boofin]
    igrid = igrid[boofin]
    xvals = np.unique(xgrid)
    yvals = np.unique(ygrid)
    xlen = len(xvals)
    ylen = len(yvals)
    xgrid = xgrid.reshape((xlen, ylen))
    ygrid = ygrid.reshape((xlen, ylen))
    igrid = igrid.reshape((xlen, ylen))
    if config.mzbins != 0:
        xgrid, ygrid, igrid = linearize_2d(xvals, yvals, igrid, config.mzbins)
    # tend = time.clock()
    # print "Time2: %.2gs" % (tend - tstart)
    if config.smooth != 0 or config.smoothdt != 0:
        igrid = smooth_2d(igrid, config.smooth, config.smoothdt)
    if config.subbuff != 0 or config.subbufdt != 0:
        igrid = subtract_complex_2d(igrid, config)
    # tend = time.clock()
    # print "Time3: %.2gs" % (tend - tstart)
    if config.detectoreffva != 0:
        igrid = detectoreff_2d(igrid, xgrid, config.detectoreffva)
    igrid /= np.amax(igrid)
    print "Shape of Processed Data: ", igrid.shape
    return xgrid, ygrid, igrid


def calc_twave_dt(mass, z, ccs, config):
    """
    Calculate drift time in a T-wave cell
    :param mass: Mass in Da
    :param z: Charge
    :param ccs: CCS in Angstrom^2
    :param config: UniDecConfig object (carries the key parameters of calibration)
    :return: Drift time (ms)
    """
    rmass = (mass * config.gasmass) / (mass + config.gasmass)
    rccs = ccs / (z * np.sqrt(1. / rmass))
    rdt = np.exp((np.log(rccs) - config.tcal2) / config.tcal1)
    out = rdt + (config.edc * np.sqrt(mass / z) / 1000.)
    return out


def calc_linear_dt(mass, z, ccs, config):
    """
    Calculate drift time in a linear cell
    :param mass: Mass in Da
    :param z: Charge
    :param ccs: CCS in Angstrom^2
    :param config: UniDecConfig object (carries the key parameters of temp, pressure, etc)
    :return: Drift time (ms)
    """
    e = 1.60217657E-19
    kb = 1.3806488E-23
    n = 2.6867774E25
    po = 760
    tempo = 273.15
    temp = config.temp + tempo
    ccsconst = (np.sqrt(18.0 * np.pi) / 16.0) * (e / np.sqrt(kb * temp)) / n * (
        config.volt / pow(config.driftlength, 2.0)) * (po / config.pressure) * (temp / tempo) * 1E20
    ac = 1.6605389E-27
    rmass = ac * (mass * config.gasmass) / (mass + config.gasmass)
    td = ccs / (ccsconst * z * np.sqrt(1 / rmass))
    dt = (td * 1000) + config.to
    return dt


def calc_native_ccs(mass, gasmass):
    """
    Predict native CCS value for a particular mass in a given mass
    :param mass: Mass of protein in Da
    :param gasmass: Mass of IM background gas in Da
    :return: Predicted CCS in Angstroms^2
    """
    if gasmass < 10:
        a = 4.06739
        b = 0.629424
    else:
        a = 5.33111
        b = 0.613072
    return a * pow(mass, b)
