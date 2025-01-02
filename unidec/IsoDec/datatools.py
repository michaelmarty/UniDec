import numpy as np
from numba import njit
import scipy.ndimage.filters as filt

# from bisect import bisect_left

# Make a python code from this C code
@njit(fastmath=True)
def fastnearest(array, target):
    """
    In a sorted array, quickly find the position of the element closest to the target.
    :param array: Array
    :param target: Value
    :return: np.argmin(np.abs(array - target))
    """
    if array.size == 0:
        return 0
    start = 0
    length = len(array) - 1
    end = 0
    diff = length - start
    while diff > 1:
        if target < array[start + (length - start) // 2]:
            length = start + (length - start) // 2
        elif target == array[start + (length - start) // 2]:
            end = start + (length - start) // 2
            length = start + (length - start) // 2
            start = start + (length - start) // 2
            return end
        elif target > array[start + (length - start) // 2]:
            start = start + (length - start) // 2
        diff = length - start
    if np.abs(target - array[start]) >= np.abs(target - array[length]):
        end = length
    else:
        end = start
    return end

@njit(fastmath=True)
def fastwithin_abstol(array, target, tol):
    result = []

    if len(array) == 0:
        return result

    nearest_idx = fastnearest(array, target)

    if np.abs(array[nearest_idx] - target) <= tol:
        result.append(nearest_idx)
    else:
        return result

    adding_upper = True
    adding_lower = True

    current_upper = nearest_idx + 1
    current_lower = nearest_idx - 1

    while adding_upper or adding_lower:
        if adding_upper:
            if current_upper >= len(array):
                adding_upper = False
            elif np.abs(array[current_upper] - target) <= tol:
                result.append(current_upper)
                current_upper += 1
            else:
                adding_upper = False

        if adding_lower:
            if current_lower < 0:
                adding_lower = False
            elif np.abs(array[current_lower] - target) <= tol:
                result.append(current_lower)
                current_lower -= 1
            else:
                adding_lower = False

    return result

@njit(fastmath=True)
def fastwithin_abstol_withnearest(array, target, tol):
    result = []

    if len(array) == 0:
        return -1, result

    nearest_idx = fastnearest(array, target)

    if np.abs(array[nearest_idx] - target) <= tol:
        result.append(nearest_idx)
    else:
        return nearest_idx,result

    adding_upper = True
    adding_lower = True

    current_upper = nearest_idx + 1
    current_lower = nearest_idx - 1

    while adding_upper or adding_lower:
        if adding_upper:
            if current_upper >= len(array):
                adding_upper = False
            elif np.abs(array[current_upper] - target) <= tol:
                result.append(current_upper)
                current_upper += 1
            else:
                adding_upper = False

        if adding_lower:
            if current_lower < 0:
                adding_lower = False
            elif np.abs(array[current_lower] - target) <= tol:
                result.append(current_lower)
                current_lower -= 1
            else:
                adding_lower = False

    return nearest_idx, result

@njit(fastmath=True)
def fastpeakdetect(data, window: int = 10, threshold=0.0, ppm=None, norm=True):
    peaks = np.empty((0, 2))
    if data is None or data.size == 0 or data.ndim != 2 or data.shape[1] != 2:
        return peaks

    length = data.shape[0]
    maxval = np.amax(data[:, 1]) if norm else 1

    if maxval == 0:
        return peaks

    for i in range(length):
        if data[i, 1] > maxval * threshold:
            if ppm is not None:
                ptmass = data[i, 0]
                newwin = ppm * 1e-6 * ptmass
                start = fastnearest(data[:, 0], ptmass - newwin)
                end = fastnearest(data[:, 0], ptmass + newwin)
            else:
                start = max(0, i - window)
                end = min(length, i + window + 1)

            if start < end:
                testmax = np.amax(data[start:end, 1])
                if data[i, 1] == testmax and np.all(data[i, 1] != data[start:i, 1]):
                    peaks = np.append(peaks, np.array([[data[i, 0], data[i, 1]]]), axis=0)

    return peaks


@njit(fastmath=True)
def fastcalc_FWHM(peak, data):
    index = fastnearest(data[:, 0], peak)
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
    #print("FWHM:", [data[indexstart, 0], data[indexend, 0]])
    return FWHM, [data[indexstart, 0], data[indexend, 0]]


def get_noise(data, n=20):
    """
    Get the noise level of the data.
    :param data: 2D numpy array of data
    :return: float, noise level
    """
    ydat = data[:, 1]
    # take n random data points from the data
    noise = np.concatenate([ydat[:n], ydat[-n:]])
    # return the standard deviation of the noise
    return np.std(noise)


def remove_noise_cdata(data, localmin=100, factor=1.5, mode="median"):
    """
    Remove noise from the data.
    :param data: 2D numpy array of data
    :param localmin: int, number of data points local width to take for min calcs
    :return: data with noise removed
    """
    if len(data) < localmin:
        return data

    ydat = data[:, 1]

    # find local minima
    if mode == "median":
        localmins = [np.median(ydat[i:i + localmin]) for i in range(len(ydat) - localmin)]
    elif mode == "min":
        localmins = [np.amin(ydat[i:i + localmin]) for i in range(len(ydat) - localmin)]
    elif mode == "mean":
        localmins = [np.mean(ydat[i:i + localmin]) for i in range(len(ydat) - localmin)]
    else:
        localmins = [np.median(ydat[i:i + localmin]) for i in range(len(ydat) - localmin)]
    # Extend the last bit at the end value to be the same length as data
    localmins = np.concatenate([localmins, np.full(len(data) - len(localmins), localmins[-1])])
    # smooth
    noiselevel = np.convolve(localmins, np.ones(localmin * 1) / localmin * 1, mode="same")
    # Subtract the noise level from the data
    newdata = data[:, 1] - noiselevel * factor

    data = data[newdata > 0]
    # return the standard deviation of the noise
    return data


@njit(fastmath=True)
def get_fwhm_peak(data, peakmz):
    """
    Get the full width half max of the peak.
    :param data: 2D numpy array of data
    :param peakmz: float, m/z value of the peak
    :return: float, full width half max of the peak
    """
    fwhm, interval = fastcalc_FWHM(peakmz, data)
    return fwhm




#@njit(fastmath=True)
def get_centroid(data, peakmz, fwhm=1):
    """
    Get the centroid of the peak.
    :param data: 2D numpy array of data
    :param peakmz: float, m/z value of the peak
    :return: float, centroid of the peak
    """
    b1 = data[:, 0] > peakmz - fwhm
    b2 = data[:, 0] < peakmz + fwhm
    b = b1 & b2
    if np.sum(b) == 0:
        return peakmz

    d = data[b]
    return np.sum(d[:, 0] * d[:, 1]) / np.sum(d[:, 1])


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

#@njit(fastmath=True)
def get_all_centroids(data, window=5, threshold=0.0001, background=100, moving_average_smoothing=3):
    if len(data) < 3:
        return np.empty((0, 2))

    if background > 0:
        data = datacompsub(data, 100)

    if moving_average_smoothing > 1:
        data[:,1] = np.convolve(data[:,1], np.ones(moving_average_smoothing)/moving_average_smoothing, mode='same')

    fwhm = get_fwhm_peak(data, data[np.argmax(data[:, 1]), 0])
    peaks = fastpeakdetect(data, window=window, threshold=threshold)

    p0 = [get_centroid(data, p[0], fwhm) for p in peaks]
    outpeaks = np.empty((len(p0), 2))
    outpeaks[:, 0] = p0
    outpeaks[:, 1] = peaks[:, 1]
    # print("Final length:", len(outpeaks))
    return outpeaks


# @njit(fastmath=True)
def get_centroids(data, peakmz, mzwindow=None):
    if mzwindow is None:
        mzwindow = [-1.5, 3.5]
    chopped = data[(data[:, 0] > peakmz + mzwindow[0]) & (data[:, 0] < peakmz + mzwindow[1])]

    if len(chopped) < 3:
        return np.array([]), np.array([])

    fwhm = get_fwhm_peak(chopped, peakmz)

    peaks = fastpeakdetect(chopped, window=5, threshold=0.01)
    # print(peaks)

    for p in peaks:
        d = chopped[(chopped[:, 0] > p[0] - fwhm) & (chopped[:, 0] < p[0] + fwhm)]
        if len(d) == 0:
            continue
        p[0] = np.sum(d[:, 0] * d[:, 1]) / np.sum(d[:, 1])
    # print(peaks)
    return peaks, chopped


def isotope_finder(data, mzwindow=1.5):
    nl = get_noise(data)
    # Chop data below noise level
    peaks = fastpeakdetect(data, window=50, threshold=nl / np.amax(data[:, 1]))
    # sort peaks
    peaks = np.array(sorted(peaks, key=lambda x: x[1], reverse=True))
    return peaks

# @njit(fastmath=True)
def check_spacings(spectrum: np.ndarray):
    spacings = np.diff(spectrum[:, 0])
    return np.median(spacings)

def simp_charge(centroids, silent=False):
    """
    Simple charge prediction based on the spacing between peaks. Picks the largest peak and looks at the spacing
    between the two nearest peaks. Takes 1/avg spacing as the charge.
    :param centroids: Centroid data, 2D numpy array as [m/z, intensity]
    :param silent: Whether to print the charge state, default False
    :return: Charge state as int
    """
    diffs = np.diff(centroids[:, 0])
    maxindex = np.argmax(centroids[:, 1])
    start = maxindex - 1
    end = maxindex + 2
    if start < 0:
        start = 0
    if end > len(diffs):
        end = len(diffs)
    centraldiffs = diffs[start:end]
    if len(centraldiffs) == 0:
        return 0
    avg = np.mean(centraldiffs)
    if avg == 0:
        return 0
    charge = 1 / avg
    charge = round(charge)
    if not silent:
        print("Simple Prediction:", charge)
    return charge


def subtract_matched_centroid_range(profile_data, matched_theoretical, centroids, noise_threshold=1000, tolerance = 0.001, window_size=5):
    mz_values = profile_data[:, 0].copy()
    intensity_values = profile_data[:, 1].copy()

    stored_widths = {}
    peak_ranges = []

    for centroid in centroids:
        centroid_x = centroid[0]

        fwhm = get_fwhm_peak(profile_data, centroid_x)
        stored_widths[centroid_x] = fwhm
        #Old
        #left_range = max((centroid_x - fwhm / 2), 0)
        left_range = max((centroid_x - fwhm),0)
        right_range = max((centroid_x + fwhm),0)
        left_range -= window_size
        right_range += window_size

        peak_ranges.append((left_range, right_range, centroid_x))
    # 0 is left
    # 1 is right
    # 2 is the peak itself
    matched_ranges = []
    theo_intensites = matched_theoretical[:, 1]
    theo_mz = matched_theoretical[:, 0]
    for peak_range in peak_ranges:
        centroid_x = peak_range[2]

        for t_mz, t_intent in zip(theo_mz, theo_intensites):
            if np.abs(t_mz - centroid_x) < tolerance:
                matched_ranges.append((peak_range[0], peak_range[1], centroid_x))
                break


    matched_ranges.sort(key=lambda x: x[2])
    for left_range, right_range, peak in matched_ranges:
        indices_to_zero = np.logical_and(mz_values >= left_range, mz_values <= right_range)
        filtered = indices_to_zero[indices_to_zero< noise_threshold]
        intensity_values[filtered] = 0


    profile_data[:, 1] = intensity_values
    return profile_data


if __name__ == "__main__":
    a = 14501.69584
    blist = [8031.29980469, 9462.05664062, 9820.18164062, 9935.2109375,
             10980.79394531, 11164.91015625, 12485.61914062, 12727.8125,
             12973.95507812, 14500.69335938, 14719.77832031]

    index = fastnearest(np.array(blist), a)

    print(index, blist[index])
