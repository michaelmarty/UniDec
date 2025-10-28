# Notes
# I had to copy over lots of definitions
# Addition of zpseq to account for zero padding in the sequence
# Inserted osFactor calculation
# Changed OS factor in definitions from 200 to an earlier-defined osFactor
# Changed invExp (line 176) from expandSequence(tempInv, osFactor=osFactor, singlePulse=True) to two separate definitions:
#       invExp = np.zeros_like(yF)
#       invExp[:osFactor * len(seq)] = expandSequence(tempInv, osFactor=osFactor, singlePulse=True)  # expand sequence
#      The original was creating a sequence that didn't match the length of the data, causing matrix multiplication errors.
# Made similar adjustment for xSeq, replacing it with xseq = np.zeros_like(sumMtx[:,0]) and xSeq[:len(pSeq)] = np.arange(len(pSeq))*timeSpacing
# Original mzVal gets data from file name, but is not relevant to the names we've been using. But it's used in the plotDemux function.

# Imports
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
from scipy import signal
from numpy.fft import fft, ifft
import sys

# Target File
dir = 'Z:\\Group Share\\ONO\\CDMS\\5 bit HT data\\20240313 JDS Bgal 0o1g_l inj1o5 bit5 zp10 s1_unidecfiles'
os.chdir(dir)
csvFiles = '20240322_TIC_chrom.csv'

dir = 'Z:\\Group Share\\ONO\\CDMS\\5 bit HT data'
os.chdir(dir)
csvFiles = 'Raw_622.csv'

fileName = csvFiles

# Load Target csv
data = np.loadtxt(fileName, delimiter = ',')
x = data[:,0]
y = data[:,1]

# Sequence of interest
seqTxt = '0000100101100111110001101110101'
seq = []
for e in seqTxt:
    seq.append(int(e))

seq = np.array(seq, dtype = np.float64)#convert sequence
#print('Sequence array:', seq)
zpseq = len(seq)+10
#print('zpseq:', zpseq)

# Calculate OS factor
osFactor = round(len(x)/zpseq)
osFactor = 200
#print('OS Factor:', osFactor)
print('len(x):', len(x))


# Plots the chromatograph
def createSeqMTX(useSeq):
    '''
    useSeq is the target sequence used for modulate (not oversampled)
    '''
    seqMTX = np.zeros([len(useSeq), len(useSeq)])
    for i in range(len(useSeq)):
        seqMTX[i,:] = np.roll(useSeq,i)
    return seqMTX
mtx = createSeqMTX(seq)
invprs = np.linalg.inv(mtx)#create inverse of the matrix
invSeq = invprs[:,0] #grab the first row
mtxUS = createSeqMTX(seq)#unscaled MTX
invMtxUS = np.linalg.inv(mtxUS)
invSeqUS = invMtxUS[:,0]
fig = plt.figure(figsize = (8,3))
#print(fileName)
#smoothing of signal
ySmth = signal.savgol_filter(y, 5, 3)
plt.plot(x,y, alpha = 0.25, label = "%s"%fileName)
plt.plot(x,ySmth)
plt.legend()
plt.ylabel('Intensity')
plt.xlabel('Scan Number')
#plt.show()


def normVector(npArray, retFactor = False):
    '''
    Normalize any input vector.
    Drive maximum to 1.0
    '''
    if retFactor:
        return npArray/npArray.max(), npArray.max()
    else:
        return npArray/npArray.max()


def expandSequence(seq, osFactor=osFactor, singlePulse=False, scaleFactor=0.5):
    '''
    Expand a PRBS by a target over sampling factor (osFactor)
    scaleFactor can be used to scale the matrix if needed.
    '''

    #seqLen = len(seq)      # I think this might be vestigial code; it's not used anywhere
   # print('len(seq):', len(seq))       # For debugging. Repeats on a loop

    newSeq = np.zeros(osFactor * len(seq))
    #newSeq = np.zeros(osFactor * zpseq)
    #print('len(newSeq) line 98:', len(newSeq))     # For debugging. Repeats on a loop
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
    #print('len(newSeq):', len(newSeq))      # For debugging. Repeats on a loop
    return newSeq

# Fourier Transform
def periodic_corr(x, y):
    """Periodic correlation, implemented using the FFT.

    x and y must be real sequences with the same length.
    """
    #print('len(x):', len(x))   # For debugging, checking sequence length. Repeats on a loop
    #print('len(y):', len(y))   # For debugging. Repeats on a loop
    return ifft(fft(x) * fft(y).conj()).real

def createMaskedTransformMTX(intVec, invMtx, bitOrder=5, osFactor=osFactor, filterType='W', maskType='R', maskWindow=5,
                             randBool=False):
    '''
    Accepts:
        The intensity vector for transform (intVec)
        Inverted matrix of the masked matrix (invMtx)
        Bit order of the sequence under test (bitOrder)
        Oversampling Factor (osFactor)
        Mask Window (see Above, Algorithm #1)

    Returns the transformed matrix.

    '''
    sumDecon = np.zeros_like(intVec)
    sLen = int(2 ** bitOrder - 1)  # sequence length of the raw masks
    colMod = 1  # not used currently but left for posterity
    sumMtx = np.zeros((len(intVec), invSeq.shape[0] * colMod))


    if filterType == 'S':  # Savitsky Golay Smoothing to minimize spikes from TOF data
        yF = signal.savgol_filter(intVec, 27, 3)  # These are hard coded and it may be good to optimize
    if filterType == None:
        yF = intVec

    # see above for a longer descriptions
    win = maskWindow
    np.random.seed(42)  # done to make sure results are reproducible
    for i in range(invSeq.shape[0] * colMod):
        tempInv = invMtx[:, i % sLen]
        if i > 0:
            if randBool:  # use for the randomization of the mask
                tempInd = np.random.randint(sLen, size=int(sLen / 2))
                tempInv[tempInd] = 0
            else:
                indA = np.random.randint(sLen)
                indB = (indA + win) % sLen
                tempInv[indA:indB] = 0
        #print("T:", tempInv, np.sum(tempInv))

        invExp = np.zeros_like(yF)
        invExp[:osFactor * len(seq)] = expandSequence(tempInv, osFactor=osFactor, singlePulse=True)  # expand sequence
        tempCC = periodic_corr(yF, invExp)  # perform cross correlation and deconvolution
        sumDecon += np.roll(normVector(tempCC), osFactor * i)
        sumMtx[:, i % sLen] += np.roll(normVector(tempCC), osFactor * i)

    return sumMtx


#osFactor = 11
#print('invSeq:', invSeq)
sumMtx = createMaskedTransformMTX(normVector(y), invMtxUS, bitOrder = 5, osFactor = osFactor, filterType = 'S', maskWindow = 5)
#print('sumMtx:', sumMtx)
plt.figure()
plt.plot(sumMtx[:,0])
#plt.show()
pSeq = expandSequence(seq, osFactor = osFactor)
timeSpacing = 0.165576 #see above
xSeq = np.zeros_like(sumMtx[:,0])
xSeq[:len(pSeq)] = np.arange(len(pSeq))*timeSpacing
#print('len(xSeq):', len(xSeq))
#print('len(sumMtx):', len(sumMtx))
#print(type(sumMtx))

# Here's where deweighting happens
def demuxFromMatrix(dmxMTX, seq, cutOff=0):
    '''
    Accepts the matrix that has been demultiplexed and contains the amplified transform errors.
    It is possible to deweight by standard deviation or by noting the number of
    sign changes (which is an indirect capture of spectral variability induced by the artifacts)

    Returns:
        1) The raw transform zeroed at points where negative signs occur
        2) The locations where the sign occurs

    '''
    posVector = np.zeros_like(seq)      #Make a list of zeros the same length as the sequence
    negLocs = np.zeros_like(seq) + 1    #Make a list of ones the same length as the sequence
    #print('PosVector1: ', posVector)
    #print('NegLocs1: ', negLocs)
    mtxShape = dmxMTX.shape[0]          #The shape (rows?) of the matrix. So it's within the range of the target matrix.
    for i in range(mtxShape):
        negCount = np.sum(dmxMTX[i, :] < cutOff, axis=0)        #Identifies where negative values are in the matrix?
        if negCount == 0:
            negLocs[i] = 0                            #If there are no negative values, set the negative locations to 0
            posVector[i] += dmxMTX[:, 0][i]           #Add the first column of the matrix to posVector, which gets plotted?
    return posVector, negLocs

dmxVec, negLocs = demuxFromMatrix(sumMtx, xSeq, 0)
#np.set_printoptions(threshold=sys.maxsize)
print('sumMtx shape:', sumMtx.shape)
print('dmxVex shape:', dmxVec.shape)
print('dmxVec:', dmxVec)
print('negLocs:', negLocs)
#make new figure window
#plt.figure()
plt.plot(dmxVec)
#plt.show()

#mzVal = fileName.split('_')[-1]
#mzVal = int(mzVal.split('.csv')[0])
mzVal = 0

#Use dmxVec if you want to see the result where large degrees of variation in the demultiplexed spectra are zeroed
#   (i.e. the quick and an less than elegant solution) (Not too bad actually).
#Use the deweighted version for comparison. The manuscript uses the deweighting version as it preserves peak shape for
#   ion populations displaying larger tails (i.e. larger m/z ions).

# Here's where things are plotted
def plotDemux(dmxVector, xVec, rawVector, mzVal='', saveBool=False):
    fig, axList = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    fs = 18
    fig.suptitle("Extracted Ion Mobility Spectrum", fontsize=fs + 2) #\nm/z: %.2f" % mzVal, fontsize=fs + 2)

    newY = rawVector
    posVector = dmxVector       # posVector Used in demuxFromMatrix(dmxMTX, seq, cutOff=0):
    ax = axList[0]

    offset = posVector.max()
    ax.plot(xVec, newY, label="Raw Transform")
    ax.fill_between(xVec, 0, newY, alpha=0.25)
    ax.plot(xVec, posVector + offset, label="WSU DEMUX")
    ax.fill_between(xSeq, offset, posVector + offset, alpha=0.25)

    #ax.set_xlabel("Drift Time (ms)", fontsize=fs)
    ax.set_ylabel("Normalized Intensity (a.u.)", fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs - 4)
    ax.legend(fontsize=fs - 8)

    ax = axList[1]
    ax.plot(xVec, newY, label="Raw Transform")
    ax.fill_between(xVec, 0, newY, alpha=0.25)
    ax.plot(xVec, posVector + offset, label="WSU DEMUX")
    ax.fill_between(xVec, offset, posVector + offset, alpha=0.25)

    #ax.set_xlabel("Drift Time (ms)", fontsize=fs)
    ax.set_ylabel("Normalized Intensity (a.u.)", fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs - 4)
    maxInd = posVector.argmax()
    scale = int(0.05 * len(posVector))
    startInd = maxInd - scale
    if startInd < 0:
        startInd = 0
    endInd = maxInd + scale
    if endInd > len(posVector):
        endInd = len(posVector)

    ax.set_xlim(xSeq[startInd], xSeq[endInd])

    fig.subplots_adjust(wspace=.25)

    if saveBool:
        plt.savefig("DEMUX_ZOOM_%d" % mzVal)


rawVector = sumMtx[:,0]
#plotDemux(dmxVec, xSeq, rawVector, mzVal = mzVal)
#plt.show()
'''
#The following import is needed if you want to intelligently deweight rather than drive to 0.
from astropy.convolution import Gaussian1DKernel, convolve
gauss = Gaussian1DKernel(stddev=20)#Setup convolution kernel for deweighting spikes/transform artifacts
locG = np.convolve(1-negLocs,gauss.array, 'same') #local gaussian to deweight areas of high variability
fixedVec = normVector(rawVector*locG)
plotDemux(fixedVec, xSeq, rawVector, mzVal = mzVal, saveBool = False)
'''