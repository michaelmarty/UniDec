

import numpy as np
from unidec.modules import unidectools as ud

__author__ = 'Michael.Marty'


def make_mass_spectrum(array, zrange=(10, 50), mzrange=(2000, 10000), mz_bin_size=1, adductmass=1.00727647, psfun=0,
                       noise=0, baseline=0, **kwargs):
    """
    Create a new mass spectrum.

    First, create the mzaxis from mzrange and the possible charge states from zrange.
    Then, add in peaks at the appropriate mz values for the parameters defined in array.
    Then, normalize and add optional background and noise.
    Finally, renormalize.

    :param array: P x 5 array of parameters [mass, mass fwhm, z avg, z std dev, intensity]
    :param zrange: Tuple of list define range of allowed charge states
    :param mzrange: Tuple or list defining range of m/z values
    :param mz_bin_size: delta mz of mzaxis
    :param adductmass: Mass of electrospray adduct species
    :param psfun: Peak shape function integer code
    :param noise: Std deviation of random Gaussian noise to be added to normalized spectrum
    :param baseline: Peak intensity of sinusoidal baseline to be added to normalized spectrum
    :param kwargs: Keyword catcher
    :return: Spectrum, ztab (new spectrum in N x 2 (m/z, intensity) and charge states allowed)
    """
    mzaxis = np.arange(mzrange[0], mzrange[1], mz_bin_size)
    ztab = np.arange(zrange[0], zrange[1], 1)

    array = np.array(array)
    num = int(len(array.flatten()) / 5)
    array = np.reshape(array, (num, 5))
    array = np.abs(array)
    output = np.zeros(len(mzaxis))
    for i in range(0, len(array)):
        mmid = array[i, 0]
        msig = array[i, 1]
        zmid = array[i, 2]
        zsig = array[i, 3]
        inten = array[i, 4]

        zint = np.array(ud.ndis_std(ztab, zmid, zsig, norm_area=True))

        mzvals = (mmid + adductmass * ztab) / ztab
        mzsigs = msig / ztab

        for j in range(0, len(ztab)):
            output += ud.make_peak_shape(mzaxis, psfun, mzsigs[j], mzvals[j], norm_area=True) * zint[j] * inten
    output /= np.amax(output)

    if baseline > 0:
        background = np.sin((mzaxis - mzrange[0]) / (mzrange[1] - mzrange[0]) * np.pi) * baseline
        output = output + background
    if baseline < 0:
        output += np.abs(baseline)
    output /= np.amax(output)

    if noise > 0:
        noisedat = np.random.normal(0, noise, size=len(output))
        output = np.abs(output + noisedat)

    output /= np.amax(output)

    if "scramble" in kwargs:
        if kwargs['scramble']:
            np.random.shuffle(output)

    return np.transpose([mzaxis, output]), ztab


def get_zrange(params):
    """
    From paramaters defining mass spectrum, pick a reasonable charge state range.
    :param params: P x 5 array of parameters [mass, mass fwhm, z avg, z std dev, intensity]
    :return: charge state range [min z , max z]
    """
    zmax = np.amax(params[:, 2]) + np.amax(params[:, 3]) * 5
    zmin = max(np.amin(params[:, 2]) - np.amax(params[:, 3]) * 3, 1)
    zrange = [int(zmin), int(zmax)]
    return zrange


def get_mzrange(params):
    """
    From paramaters defining mass spectrum, pick a reasonable m/z range.
    :param params: P x 5 array of parameters [mass, mass fwhm, z avg, z std dev, intensity]
    :return: m/z range [min m/z , max m/z]
    """
    zrange = get_zrange(params)
    mmax = np.amax(params[:, 0])
    mmin = np.amin(params[:, 0])
    mzrange = [mmin / zrange[1], mmax / zrange[0]]
    return mzrange


def simple_params(masslist, intlist=None, resolution=1000, zwidth=2, rlist=None, **kwargs):
    """
    Set simple parameters [mass, mass fwhm, z avg, z std dev, intensity] for a mass list and optional arguments.
    Average charge is defined as the native charge state. Other parameters are defaults.

    :param masslist: list of masses in Da
    :param intlist: list of intensities (default is None which gives all as 1)
    :param resolution: Mass resolution (delta m/m), Overridden if rlist is not None
    :param zwidth: Standard deviation of charge state distribution
    :param rlist: List of resolution values for each corresponding mass.
    If None, will use a single resolution for all defined by resolution.
    :param kwargs: Extra keywords
    :return: P x 5 array of parameters [mass, mass fwhm, z avg, z std dev, intensity]
    """
    if intlist is None:
        intlist = np.ones(len(masslist))
    if rlist is None:
        rlist = np.ones(len(masslist)) * resolution
    params = np.array(
            [[m, m / float(rlist[i]), ud.predict_charge(m), zwidth, intlist[i]] for i, m in enumerate(masslist)])
    return params


def simple_spectrum(masslist, **kwargs):
    """
    For a mass list, create the parameters and ranges for a simple default spectrum.
    :param masslist: Mass values
    :param kwargs: Keywords of parameters to pass to simple_params and make_mass_spectrum
    :return: N x 2 spectrum (m/z, intensity)
    """
    params = simple_params(masslist, **kwargs)
    zrange = get_zrange(params)
    mzrange = get_mzrange(params)
    spec = make_mass_spectrum(params, zrange=zrange, mzrange=mzrange, **kwargs)
    return spec


def simple_spectrum2(masslist, **kwargs):
    """
    For a mass list, create the parameters and ranges for a simple default spectrum.
    :param masslist: Mass values
    :param kwargs: Keywords of parameters to pass to simple_params and make_mass_spectrum
    :return: N x 2 spectrum (m/z, intensity)
    """
    params = simple_params(masslist, **kwargs)
    zrange = get_zrange(params)
    mzrange = get_mzrange(params)
    spec = make_mass_spectrum(params, zrange=zrange, mzrange=mzrange, **kwargs)
    return spec, params



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    dat, ztab = simple_spectrum([200000, 205000], intlist=[1, 2], resolution=500, psfun=2)
    plt.plot(dat[:, 0], dat[:, 1])
    plt.show()
