import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import scipy as sp
from scipy import signal
from numpy.fft import fft, ifft
import math


def periodic_corr(x, y):
    """Periodic correlation, implemented using the FFT.

    x and y must be real sequences with the same length.
    """
    return ifft(fft(x) * fft(y).conj()).real


def normVector(npArray, retFactor = False):
    '''
    Normalize any input vector.
    Drive maximum to 1.0
    '''
    if retFactor:
        return npArray/npArray.max(), npArray.max()
    else:
        return npArray/npArray.max()


def createSeqMTX(useSeq):
    '''
    useSeq is the target sequence used for modulate (not oversampled)
    '''
    seqMTX = np.zeros([len(useSeq), len(useSeq)])
    for i in range(len(useSeq)):
        seqMTX[i,:] = np.roll(useSeq,i)
    return seqMTX


def expandSequence(seq, osFactor=200, singlePulse=False, scaleFactor=0.5):
    '''
    Expand a PRBS by a target over sampling factor (osFactor)
    scaleFactor can be used to scale the matrix if needed.
    '''

    seqLen = len(seq)
    newSeq = np.zeros(osFactor * len(seq))
    j = -1
    maxVal = seq.max() * scaleFactor
    minVal = seq.min() * scaleFactor
    for i, v in enumerate(newSeq):  # go through each element
        if i % osFactor == 0:  # reset count when the oversampling factor is reached
            j += 1  # start counting upward
            curVal = seq[j]

            if singlePulse:
                if curVal > maxVal:
                    newSeq[i] = curVal
                    curVal = 0
                elif curVal < minVal:

                    newSeq[i] = curVal
                    curVal = 0
                else:
                    curVal = 0

            else:
                newSeq[i] = curVal
        else:
            newSeq[i] = curVal

    return newSeq

def next_power_of_2(x):
    return 1 if x == 0 else 2**math.ceil(math.log2(x))

def next_multiple_of_x(x, y):
    return x*math.ceil(y/x)

def spline_resample(x, y, osLen):
    interp_func = sp.interpolate.interp1d(x, y, kind='cubic')
    interpX = np.linspace(x[0], x[-1], osLen)
    interpY = interp_func(interpX)
    return interpY