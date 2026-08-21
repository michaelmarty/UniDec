"""Small numerical and file helpers formerly imported from UniDec."""

from copy import deepcopy
import os

import numpy as np
from numba import njit
from scipy.ndimage import gaussian_filter


def match_files_recursive(topdir, ending=".raw"):
    matches = []
    for root, _, files in os.walk(topdir):
        matches.extend(os.path.join(root, name) for name in files if name.endswith(ending))
    return np.asarray(matches)


def isempty(value):
    if value is None:
        return True
    try:
        return np.asarray(value, dtype=object).size == 0
    except (TypeError, ValueError, AttributeError):
        return False


def round_to_nearest(value, step):
    remainder = value % step
    return value + step - remainder if remainder * 2 >= step else value - remainder


def safedivide(numerator, denominator):
    numerator = np.asarray(numerator)
    denominator = np.asarray(denominator)
    output = np.zeros(np.broadcast_shapes(numerator.shape, denominator.shape), dtype=np.result_type(numerator, denominator, float))
    return np.divide(numerator, denominator, out=output, where=denominator != 0)


def datachop(data, lower, upper):
    data = np.asarray(data)
    return data[(data[:, 0] >= lower) & (data[:, 0] <= upper)]


def datacompsub(data, width):
    data = np.asarray(data)
    if len(data) == 0:
        return data
    width = abs(int(width))
    if width == 0:
        return data
    output = data.copy()
    minima = np.empty(len(output))
    for index in range(len(output)):
        start = max(0, index - width)
        stop = min(index + width, len(output))
        minima[index] = np.min(output[start:stop, 1])
    output[:, 1] -= gaussian_filter(minima, width * 2)
    return output


def dataremove(data, lower, upper):
    output = np.asarray(data).copy()
    output[(output[:, 0] >= lower) & (output[:, 0] <= upper), 1] = 0
    return output


@njit(fastmath=True)
def within_ppm(theoretical, experimental, tolerance):
    return np.abs(((theoretical - experimental) / theoretical) * 1e6) <= tolerance


def _nearest(array, target):
    index = np.searchsorted(array, target)
    if index == 0:
        return 0
    if index == len(array):
        return len(array) - 1
    return index if abs(array[index] - target) < abs(array[index - 1] - target) else index - 1


def lintegrate(data, new_x, fastmode=False):
    data = np.asarray(data)
    new_x = np.asarray(new_x)
    if len(new_x) < 2:
        raise ValueError("new_x must contain at least two points")
    if fastmode:
        bins = np.append(new_x, new_x[-1] + np.diff(new_x)[-1]) - np.diff(new_x)[0] / 2
        values, _ = np.histogram(data[:, 0], bins=bins, weights=data[:, 1])
        return np.column_stack((new_x, values))

    new_y = np.zeros_like(new_x, dtype=float)
    for x_value, intensity in data:
        if not new_x[0] < x_value < new_x[-1]:
            continue
        index = _nearest(new_x, x_value)
        if new_x[index] == x_value:
            new_y[index] += intensity
        elif new_x[index] < x_value and index < len(new_x) - 1:
            fraction = (x_value - new_x[index]) / (new_x[index + 1] - new_x[index])
            new_y[index] += (1 - fraction) * intensity
            new_y[index + 1] += fraction * intensity
        elif index > 0:
            fraction = (x_value - new_x[index - 1]) / (new_x[index] - new_x[index - 1])
            new_y[index - 1] += (1 - fraction) * intensity
            new_y[index] += fraction * intensity
    return np.column_stack((new_x, new_y))


__all__ = [
    "datachop", "datacompsub", "dataremove", "isempty", "lintegrate",
    "match_files_recursive", "round_to_nearest", "safedivide", "within_ppm",
]
