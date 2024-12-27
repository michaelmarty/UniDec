import re
import numpy as np
import unidec.tools as ud
from copy import deepcopy

def header_test(path, deletechars=None, delimiter=" |\t|,", strip_end_space=True):
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


def get_resolution(testdata):
    """
    Get the median resolution of 1D MS data.
    :param testdata: N x 2 data (mz, intensity)
    :return: Median resolution (float)
    """
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    resolutions = ud.safedivide(diffs[:, 0], diffs[:, 1])
    return np.median(resolutions)


def fit_line(x, a, b):
    return a * x ** b


def get_longest_index(datalist):
    lengths = [len(x) for x in datalist]
    return np.argmax(lengths)

def merge_spectra(datalist, mzbins=None, type="Interpolate"):
    """
    Merge together a list of data.
    Interpolates each data set in the list to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    # Filter out junk spectra that are empty
    datalist = [x for x in datalist if len(x) > 0]
    # Find which spectrum in this list is the largest. This will likely have the highest resolution.
    maxlenpos = get_longest_index(datalist)

    # Concatenate everything for finding the min/max m/z values in all scans
    concat = np.concatenate(datalist)
    # xvals = concat[:, 0]
    # print "Median Resolution:", resolution
    # axis = nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    # for d in datalist:
    #    print(d)
    # If no m/z bin size is specified, find the average resolution of the largest scan
    # Then, create a dummy axis with the average resolution.
    # Otherwise, create a dummy axis with the specified m/z bin size.

    if mzbins is None or float(mzbins) == 0:
        resolution = get_resolution(datalist[maxlenpos])
        if resolution < 0:
            print("ERROR with auto resolution:", resolution, maxlenpos, datalist[maxlenpos])
            print("Using ABS")
            resolution = np.abs(resolution)
        elif resolution == 0:
            print("ERROR, resolution is 0, using 20000.", maxlenpos, datalist[maxlenpos])
            resolution = 20000

        axis = ud.nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        axis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]), float(mzbins))
    template = np.transpose([axis, np.zeros_like(axis)])

    print("Length merge axis:", len(template))

    # Loop through the data and resample it to match the template, either by integration or interpolation
    # Sum the resampled data into the template.
    for d in datalist:
        if len(d) > 2:
            if type == "Interpolate":
                newdat = ud.mergedata(template, d)
            elif type == "Integrate":
                newdat = ud.lintegrate(d, axis)
            else:
                print("ERROR: unrecognized trdtrmerge spectra type:", type)
                continue

            template[:, 1] += newdat[:, 1]

    # Trying to catch if the data is backwards
    try:
        if template[1, 0] < template[0, 0]:
            template = template[::-1]
    except:
        pass
    return template

def get_resolution_im(data):
    """
    Get the median resolution of 1D MS data.
    :param testdata: N x 2 data (mz, intensity)
    :return: Median resolution (float)
    """
    testdata = deepcopy(data)
    testdata = testdata[testdata[:, 2] > 0]
    diffs = np.transpose([testdata[1:, 0], np.diff(testdata[:, 0])])
    b1 = diffs[:,1] > 0
    diffs = diffs[b1]
    resolutions = ud.safedivide(diffs[:, 0], diffs[:, 1])
    return np.median(resolutions)

def merge_im_spectra(datalist, mzbins=None, type="Integrate"):
    """
    Merge together a list of ion mobility data.
    Interpolates each data set in the list to a new nonlinear axis with the median resolution of the first element.
    Optionally, allows mzbins to create a linear axis with each point spaced by mzbins.
    Then, adds the interpolated data together to get the merged data.
    :param datalist: M x N x 2 list of data sets
    :return: Merged N x 2 data set
    """
    # Find which spectrum in this list is the largest. This will likely have the highest resolution.
    maxlenpos = get_longest_index(datalist)

    # Concatenate everything for finding the min/max m/z values in all scans
    if len(datalist) > 1:
        concat = np.concatenate(datalist)
    else:
        concat = np.array(datalist[0])
    # If no m/z bin size is specified, find the average resolution of the largest scan
    # Then, create a dummy axis with the average resolution.
    # Otherwise, create a dummy axis with the specified m/z bin size.
    if mzbins is None or float(mzbins) == 0:
        resolution = get_resolution_im(datalist[maxlenpos])
        mzaxis = ud.nonlinear_axis(np.amin(concat[:, 0]), np.amax(concat[:, 0]), resolution)
    else:
        mzaxis = np.arange(np.amin(concat[:, 0]), np.amax(concat[:, 0]), float(mzbins))

    # For drift time, use just unique drift time values. May need to make this fancier.
    dtaxis = np.sort(np.unique(concat[:, 1]))

    # Create the mesh grid from the new axes
    X, Y = np.meshgrid(mzaxis, dtaxis, indexing="ij")

    template = np.transpose([np.ravel(X), np.ravel(Y), np.ravel(np.zeros_like(X))])
    print("Shape merge axis:", X.shape)
    xbins = deepcopy(mzaxis)
    xbins[1:] -= np.diff(xbins)
    xbins = np.append(xbins, xbins[-1] + np.diff(xbins)[-1])
    ybins = deepcopy(dtaxis)
    ybins[1:] -= np.diff(ybins) / 2.
    ybins = np.append(ybins, ybins[-1] + np.diff(ybins)[-1])
    # Loop through the data and resample it to match the template, either by integration or interpolation
    # Sum the resampled data into the template.
    for d in datalist:
        if len(d) > 2:
            if type == "Interpolate":
                newdat = ud.mergedata2d(template[:, 0], template[:, 1], d[:, 0], d[:, 1], d[:, 2])
            elif type == "Integrate":
                newdat, xedges, yedges = np.histogram2d(d[:, 0], d[:, 1], bins=[xbins, ybins], weights=d[:, 2])
            else:
                print("ERROR: unrecognized merge spectra type:", type)
                continue
            template[:, 2] += np.ravel(newdat)
    return template


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
    i += i / fit_line(i, res[0], res[1])
    while i < end:
        axis.append(i)
        i += i / fit_line(i, res[0], res[1])
    return np.array(axis)

#
# def waters_convert2(path, config=None, outfile=None, time_range=None):
#     data = WDI(path).get_data(time_range=time_range)
#
#     if outfile is None:
#         outfile = os.path.join(path, "converted_rawdata.txt")
#     np.savetxt(outfile, data)
#     return data
