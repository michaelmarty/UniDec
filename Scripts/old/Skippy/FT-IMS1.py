import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy as sp
from scipy import signal
import math

def next_power_of_2(x):
    return 1 if x == 0 else 2**math.ceil(math.log2(x))


# import text file
path = "Z:\\Group Share\\Skippy\\Projects\\FT IM CD MS\\01302024_GDH_ss3_r15_5t500Hz.txt"
#path = "Z:\\Group Share\\Skippy\\Projects\\FT IM CD MS\\TAA 5-5k.txt"
data = np.loadtxt(path, delimiter = '\t')
x = data[:,0]
y = data[:,1]



# Smooth and plot raw data
ySmth = signal.savgol_filter(y, 15, 3)
#plt.plot(x, y, alpha = 0.25, label = "Raw Data")
#plt.plot(x, ySmth)
#plt.legend()
#plt.show()

interpLen = int(next_power_of_2(len(ySmth))/2)
#print(interpLen)

# spline resample ysmth to interpLen
interp_func = sp.interpolate.interp1d(x, ySmth, kind='cubic')
interpX = np.linspace(x[0], x[-1], interpLen)
interpY = interp_func(interpX)
#print(len(interpY))


def FT_decon(x, y, startFreq, endFreq, scnComp=15, flattenTIC=True, apodize=False, plot=True, aFT=False, tOffset=0):

    if flattenTIC:
        # fit trendline to y and subtract to eliminate low frequency components
        ytrnd = signal.savgol_filter(y, 500, 3)
        y = y - ytrnd

    # compress data
    if scnComp > 1:
        ycomp = np.zeros(len(y)//scnComp)
        xcomp = np.zeros(len(x)//scnComp)
        for i in range(len(ycomp)):
            ycomp[i] = np.mean(y[i*scnComp:(i+1)*scnComp])
            xcomp[i] = x[i*scnComp]
        y = ycomp
        x = xcomp

   
    if apodize:
        # create hanning window
        N = len(y)*2
        hanning = np.hanning(N)
        y = y * hanning[len(y):]

    if plot:
        plt.subplot(2, 1, 1)
        plt.plot(x, y)
        # label x axis
        plt.xlabel('Time (min)')
        plt.subplots_adjust(hspace=0.3)

    # Fourier Transform
    Y = np.fft.fft(y)
    #print(np.real(Y[0:10]))
    #xmag = np.abs(Y[0:len(Y)//2])
    #xre = np.real(X[0:len(X)//2])

    sweepTime = x[-1] * 60
    sweepRate = (endFreq - startFreq)/sweepTime

    freq = np.fft.fftfreq(len(x), d=(x[1]-x[0])*60)
    #freq = (freq[0:len(freq)//2])

    # Create time axis
    driftTime = freq / sweepRate

    if aFT:
        deltaT = startFreq/sweepRate
        #print(len(freq))
        Y = Y*np.exp(-1j*(deltaT+tOffset)*np.pi*freq)
        Y = np.real(Y)
    else:
        Y = abs(Y)

    Y = Y[0:len(Y)//2]
    driftTime = driftTime[0:len(driftTime)//2]
    #print(len(driftTime))
    #print(len(Y))
    if plot:
        plt.subplot(2, 1, 2)
        plt.plot(driftTime, Y)
        # label x axis
        plt.xlabel('Drift Time (s)')
        plt.show()
    return driftTime, Y



#time, amp = FT_decon(interpX, interpY, 5, 500)
#time, amp = FT_decon(x, y, 5, 5000, scnComp=1, aFT=False, tOffset=40, plot=True, flattenTIC=True)
#plt.plot(time,  amp)
# set x range to 0 to 0.01
#plt.xlim(0, 0.5)
#plt.show()

def aFT_IMS(x, y, startFreq, endFreq, plot=True, nzp=0, apodize=True, flattenTIC=True, tdMax=0.2):

    if flattenTIC:
        # fit trendline to y and subtract to eliminate low frequency components
        ytrnd = signal.savgol_filter(y, 500, 3)
        y = y - ytrnd

    if apodize:
        # create half hanning window
        N = len(y)*2
        hanning = np.hanning(N)
        y = y * hanning[len(y):]

    if nzp > 0:
        # zero fill
        y = np.concatenate((y, np.zeros(len(y)*(nzp-1))))
        x = np.linspace(x[0], x[-1], len(y))

    if plot:
        plt.subplot(3, 1, 1)
        plt.plot(x, y)
        # label x axis
        plt.xlabel('Time (min)')
        plt.subplots_adjust(hspace=0.3)

    # Fourier Transform
    Y = np.fft.fft(y)

    # Create magnitude spectrum for reference
    mag = np.abs(Y)

    # Pick peaks in the magnitude spectrum
    peaks = sp.signal.find_peaks(mag[:int(len(mag)/2)], height=1, distance=1000, prominence=0.1, threshold=1000)

    # Plot magnitude spectrum with peaks


    print(peaks[0])

    # Calculate phase of Y at index of peaks
    phase = np.angle(Y[peaks[0]])


    sweepTime = x[-1] * 60
    sweepRate = (endFreq - startFreq)/sweepTime

    freq = np.fft.fftfreq(len(x), d=(x[1]-x[0])*60)

    # Create drift time axis
    driftTime = (freq / sweepRate)  # * (nzp-1)

    if plot:
        plt.subplot(3, 1, 2)
        plt.plot(driftTime[:int(len(x)/2)], mag[:int(len(mag)/2)])
        plt.plot(driftTime[peaks[0]], mag[peaks[0]], 'ro')
        plt.xlabel('Drift Time (s)')
        plt.xlim(0, tdMax)

    deltaT = startFreq/sweepRate
    #print(len(freq))
    #Y = Y*np.exp(-1j*(deltaT)*np.pi*freq)
    Y = Y * np.exp(-1j * phase[0])
    Y = np.real(Y)

    Y = Y[0:len(Y)//2]
    driftTime = driftTime[0:len(driftTime)//2]
    #print(len(driftTime))
    #print(len(Y))


    if plot:
        plt.subplot(3, 1, 3)
        plt.plot(driftTime, Y)
        # label x axis
        plt.xlabel('Drift Time (s)')
        plt.xlim(0, tdMax)
        plt.show()
        # set x axis range to 0 to 0.15

    return driftTime, Y

    return driftTime, Y

time, amp = aFT_IMS(interpX, interpY, 5, 500, plot=True)