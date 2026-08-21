"""Peak-width utilities used by IsoDec.

These focused helpers replace the dependency on UniDec's much larger
``fwhmtools`` module.
"""

import numpy as np

from .tools import safedivide


def _interpolate_x(point1, point2, y_value):
    x1, y1 = point1
    x2, y2 = point2
    if y1 == y2:
        return (x1 + x2) / 2
    return x1 + (y_value - y1) * (x2 - x1) / (y2 - y1)


def _sigma(data):
    if len(data) < 2 or np.max(data[:, 1]) <= 0:
        return -1.0
    area = np.trapz(data[:, 1], data[:, 0])
    return area / (np.sqrt(2 * np.pi) * np.max(data[:, 1]))


def _find_fwhm(data, index, wfactor=4, maxppmtol=1000, maxasymmetry=4):
    default = np.full(9, -1.0)
    mz, intensity = data[index]
    if intensity <= 0:
        return default

    def within_tolerance(value):
        return maxppmtol is None or abs(value - mz) <= abs(mz) * maxppmtol * 1e-6

    half_max = intensity / 2
    left = index
    while left > 0 and data[left, 1] >= half_max and within_tolerance(data[left, 0]):
        left -= 1
    right = index
    while right < len(data) - 1 and data[right, 1] >= half_max and within_tolerance(data[right, 0]):
        right += 1
    if left == index or right == index:
        return default

    left_mz = _interpolate_x(data[left], data[left + 1], half_max)
    right_mz = _interpolate_x(data[right - 1], data[right], half_max)
    left_width = mz - left_mz
    right_width = right_mz - mz
    if left_width <= 0 or right_width <= 0:
        return default
    if maxasymmetry is not None and max(left_width / right_width, right_width / left_width) > maxasymmetry:
        return default

    centroid_slice = data[left:right + 1]
    centroid = np.average(centroid_slice[:, 0], weights=centroid_slice[:, 1])
    broad_left = max(0, index - int(wfactor) * (index - left))
    broad_right = min(len(data) - 1, index + int(wfactor) * (right - index))
    sigma = _sigma(data[broad_left:broad_right + 1])
    return np.array([
        left_width + right_width,
        left_width,
        right_width,
        broad_left,
        broad_right,
        sigma,
        centroid,
        intensity,
        right - left,
    ])


def fast_fwhm(data, peaks, sort=False, wfactor=4, maxppmtol=1000):
    """Calculate width information for each requested peak."""
    data = np.asarray(data, dtype=float)
    peaks = np.asarray(peaks, dtype=float)
    if sort:
        data = data[np.argsort(data[:, 0])]
        peaks = peaks[np.argsort(peaks[:, 0])]
    indexes = np.searchsorted(data[:, 0], peaks[:, 0]).clip(0, len(data) - 1)
    previous = np.maximum(indexes - 1, 0)
    use_previous = np.abs(data[previous, 0] - peaks[:, 0]) <= np.abs(data[indexes, 0] - peaks[:, 0])
    indexes[use_previous] = previous[use_previous]
    return np.asarray([
        _find_fwhm(data, int(index), wfactor=wfactor, maxppmtol=maxppmtol)
        if peak[1] != 0 else np.full(9, -1.0)
        for index, peak in zip(indexes, peaks)
    ])


def remove_peaks(data, isotope_distribution, wfactor=6, maxppmtol=1000):
    """Remove profile-data regions occupied by an isotope distribution."""
    widths = fast_fwhm(data, isotope_distribution, wfactor=wfactor, maxppmtol=maxppmtol)
    keep = np.ones(len(data), dtype=bool)
    for width in widths:
        if width[0] > 0:
            keep[int(width[3]):int(width[4])] = False
    return np.asarray(data)[keep]


def ndis_std(x, mid, sig, a=1.0):
    """Evaluate a Gaussian parameterized by standard deviation."""
    return a * np.exp(-((np.asarray(x) - mid) ** 2) / (2 * sig ** 2))


def intensity_decon(data, distributions, n=10):
    """Iteratively apportion observed intensity among simulated peaks."""
    distributions = [np.asarray(distribution, dtype=float).copy() for distribution in distributions]
    observed = np.maximum(np.asarray(data)[:, 1], 0)
    simulated_sum = np.sum(distributions, axis=0)
    for _ in range(n):
        ratios = safedivide(observed, simulated_sum)
        for index, distribution in enumerate(distributions):
            total = np.sum(distribution)
            if total > 0:
                distributions[index] *= np.sum(ratios * distribution) / total
        simulated_sum = np.sum(distributions, axis=0)
    return simulated_sum


__all__ = ["fast_fwhm", "intensity_decon", "ndis_std", "remove_peaks"]
