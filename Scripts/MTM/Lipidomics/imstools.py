from pyopenms import *
import numpy as np
from copy import deepcopy


def get_im_bins(exp, max=200):
    dtall = []
    for s in exp:
        fda = s.getFloatDataArrays()[0]
        dtdata = np.unique(fda.get_data())
        dtall.append(dtdata)
        unique = np.unique(dtall)
        if len(unique) == max:
            return unique
    return unique


def slice(mzdata, dtdata, intdat, mzslice=None, dtslice=None):
    if mzslice is None and dtslice is None:
        # print("Warning: No values for either mzslice or dtslice")
        pass

    if mzslice is not None:
        b1 = mzdata < mzslice[1]
        b2 = mzdata > mzslice[0]
        b3 = b1 * b2
        mzdata = mzdata[b3]
        dtdata = dtdata[b3]
        intdat = intdat[b3]

    if dtslice is not None:
        b1 = dtdata < dtslice[1]
        b2 = dtdata > dtslice[0]
        b3 = b1 * b2
        mzdata = mzdata[b3]
        dtdata = dtdata[b3]
        intdat = intdat[b3]

    return mzdata, dtdata, intdat


def merge_spec(mzdata, dtdata, intdat, mz=False, dt=False, dttemplate=None):
    if mz and dt:
        return np.sum(intdat)

    if dttemplate is None:
        dttemplate = np.sort(np.unique(dtdata))
    if len(dttemplate) == 0:
        raise Exception

    if mz:
        merger = SpectraMerger()
        params = merger.getDefaults()
        # params.setValue("mz_binning_width", float(mz_binning_width))
        # params.setValue("mz_binning_width_unit", b"ppm" if is_ppm else b"Da")
        params.setValue("block_method:rt_block_size", 201)
        # params.setValue("block_method:ms_levels", [ms_level])
        merger.setParameters(params)

        exp = MSExperiment()
        for d in dttemplate:
            s = MSSpectrum()
            b1 = dtdata == d
            mzspec = mzdata[b1]
            intspec = intdat[b1]
            s.set_peaks((mzspec, intspec))
            exp.addSpectrum(s)
        exp.updateRanges()
        merger.mergeSpectraBlockWise(exp)
        mzout, intout = exp[0].get_peaks()
        return mzout, intout

    elif dt:
        newint = [np.sum(intdat[dtdata == d]) for d in dttemplate]
        return dttemplate, np.array(newint)

    elif not dt and not mz:
        print("Not slicing. Need to set either mz, dt, or both")
        return False


def merge_rt(exp):
    merger = SpectraMerger()
    params = merger.getDefaults()
    # params.setValue("mz_binning_width", float(mz_binning_width))
    # params.setValue("mz_binning_width_unit", b"ppm" if is_ppm else b"Da")
    params.setValue("block_method:rt_block_size", len(list(exp)) + 1)
    # params.setValue("block_method:ms_levels", [ms_level])
    merger.setParameters(params)

    exp.updateRanges()
    merger.mergeSpectraBlockWise(exp)
    mzout, intout = exp[0].get_peaks()
    return mzout, intout


def get_slice(s, mzslice=None, dtslice=None):
    fda = s.getFloatDataArrays()[0]
    dtdata = np.array(fda.get_data())
    mzdata, intdata = np.array(s.get_peaks())
    if mzslice is not None or dtslice is not None:
        mzdata, dtdata, intdata = slice(mzdata, dtdata, intdata, mzslice=mzslice, dtslice=dtslice)
    return mzdata, dtdata, intdata


def merge_to_mz(s, dtslice=None):
    mzdata, dtdata, intdata = get_slice(s, dtslice=dtslice)
    mz, it = merge_spec(mzdata, dtdata, intdata, mz=True)
    s2 = deepcopy(s)
    s2.set_peaks((mz, it))
    return s2


def merge_to_dt(s, mzslice=None):
    mzdata, dtdata, intdata = get_slice(s, mzslice=mzslice)
    dt, it = merge_spec(mzdata, dtdata, intdata, dt=True)
    s2 = MSSpectrum()
    s2.setRT(s.getRT())
    s2.set_peaks((dt, it))
    return s2


def centroid(s):
    picker = PeakPickerHiRes()
    news = MSSpectrum()
    picker.pick(s, news)
    return news
