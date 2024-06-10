import time

import numpy as np
import scipy.ndimage
from unidec.modules.fitting import *

from unidec.modules.CDEng import *
from unidec.modules.ChromEng import *
from unidec.modules.IM_functions import calc_linear_ccs, calc_linear_ccsconst
import unidec.tools as ud
from unidec.modules.unidecstructure import UniDecConfig

import matplotlib.pyplot as plt
import scipy.fft as fft
import math


# import warnings
# warnings.filterwarnings("error")


# import fast_histogram
# import pyfftw.interfaces.numpy_fft as fft

class HTEng:
    def __init__(self, *args, **kwargs):
        """
        Initialize the HTEng class. This class is used for handling Hadamard Transform (HT) related operations.
        :param args: Arguments (currently unused)
        :param kwargs: Keyword Arguments (currently unused)
        :return: None
        """

        super().__init__(*args, **kwargs)
        if not hasattr(self, "config"):
            self.config = UniDecConfig()
        self.config.htmode = True

        # Empty Arrays
        self.scans = []
        self.fullscans = []
        self.fulltime = []
        self.fulltic = []
        self.htkernel = []
        self.fftk = []
        self.htoutput = []
        self.indexrange = [0, 0]
        self.decontime = []
        self.deconscans = []
        self.dtype = float

        # Index Values
        self.padindex = 0
        self.shiftindex = 0
        self.kernelroll = 0
        self.cycleindex = 0
        self.rollindex = 0
        self.config.HTcycleindex = -1

        # Important Parameters that Should be Automatically Set by File Name
        self.config.htbit = 3  # HT bit, should automatically grab this
        self.config.htseq = ud.HTseqDict[str(self.config.htbit)]  # HT Sequence, should automatically set this

        self.config.HTcycletime = 2.0  # Total cycle time, should automatically grab this
        self.config.HTtimepad = 0  # Length of time pad, should automatically grab this

        # Important Parameters that the User may have to Set
        self.config.HTanalysistime = 38  # total analysis time in min, should automatically grab this but not for CDMS
        self.config.HTtimeshift = 7  # How far back in time to shift
        self.config.HTksmooth = 0  # Kernel Smooth

        # print("HT Engine")

    def parse_file_name(self, path):
        """
        Parse the file name to extract relevant parameters for HT processing.
        :param path: File path
        :return: None
        """
        path = os.path.split(path)[1]
        if "cyc" in path:
            # Find location of cyc and take next character
            cycindex = path.find("cyc")
            cyc = int(path[cycindex + 3])
            self.config.HTcycletime = float(cyc)

        # if "zp" in path:
        # Find location of zp and take next character
        #    zpindex = path.find("zp")
        #    zp = float(path[zpindex + 2])
        #    self.config.HTtimepad = zp * self.config.HTcycletime

        if "bit" in path:
            # Find location of bit and take next character
            bitindex = path.find("bit")
            try:
                firstbit = path[bitindex + 3]
            except Exception as e:
                print("Error setting HT bit:", e, bitindex, path[bitindex:bitindex + 4])
                firstbit = ""
            if firstbit.isdigit():
                self.config.htbit = int(firstbit)
            elif firstbit == "-":
                self.config.htbit = -1 * int(path[bitindex + 4])
            else:
                try:
                    self.config.htbit = int(path[bitindex + 3])
                except Exception as e:
                    print("Error setting HT bit:", e, bitindex, path[bitindex:bitindex + 4])
            try:
                self.config.htseq = ud.HTseqDict[str(self.config.htbit)]
            except KeyError as e:
                print("ERROR: Bit not found in dictionary. Using 3 bit sequence.")
                print(e)
                self.config.htseq = '1110100'
                self.config.htbit = 3

        if "frange" in path:
            # Split the file name by _
            parts = path.split("_")
            # Find the part that contains frange
            frangeindex = [i for i, s in enumerate(parts) if "frange" in s]
            if len(frangeindex) > 0:
                frangeindex = frangeindex[0]
                # get the next two fields
                ftstart = parts[frangeindex + 1]
                ftend = parts[frangeindex + 2]

                # Convert to float and set to field in config
                try:
                    self.config.FTstart = float(ftstart)
                except Exception as e:
                    print("Error setting FTstart:", e, ftstart)
                try:
                    self.config.FTend = float(ftend)
                except Exception as e:
                    print("Error setting FTend:", e, ftend)
                print("FT Start:", self.config.FTstart, "FT End:", self.config.FTend)

        print("Cycle Time:", self.config.HTcycletime, "Time Pad:", self.config.HTtimepad,
              "HT Sequence:", self.config.htseq)

    def setup_ht(self, cycleindex=None, seq=None):
        """
        Sets up the HT kernel for deconvolution. This is a binary sequence that is convolved with the data to
        deconvolve the data. The sequence is defined by the htseq variable. The timepad and timeshift variables
        shift the sequence in time. The config.HTksmooth variable smooths the kernel to reduce ringing.
        :param cycleindex: The length of a cycle in number of scans.
        If not specified, default is the full number of scans divided by the number of cycles.
        :return: None
        """
        if seq is None:
            self.config.htseq = ud.HTseqDict[str(self.config.htbit)]
        else:
            self.config.htseq = seq
        # Finds the number of scans that are padded
        if self.config.HTtimepad > 0:
            self.padindex = self.set_timepad_index(self.config.HTtimepad)

        if self.config.HTtimeshift >= 0:
            # Index for the first scan above the timepad
            padindex2 = self.set_timepad_index(self.config.HTtimeshift)
            # Describes the number of scans to shift the data back by
            self.shiftindex = self.padindex - padindex2

        # Check for negative shift index and fix if necessary
        if self.shiftindex < 0:
            self.shiftindex = 0
            print("Shift index is less than pad index. Check timepad and timeshift values. "
                  "Likely need to decrease timeshift. Setting to 0.")

        # Convert sequence to array
        seqarray = np.array([int(s) for s in self.config.htseq])

        seqarray = seqarray * 2 - 1
        # For any less than -1, set equal to 0
        seqarray[seqarray < -1] = 0

        shortkernel = seqarray
        # print(shortkernel)

        # Roll short kernel so that the first 1 is in index 0
        self.kernelroll = -np.argmax(shortkernel)
        shortkernel = np.roll(shortkernel, self.kernelroll)

        scaledkernel = np.arange(len(shortkernel)) / (len(shortkernel))
        # print("Scaled Kernel:", scaledkernel)
        kernelscans = scaledkernel * (np.amax(self.fullscans) - self.padindex)
        # print("Kernel Scans:", kernelscans)

        if cycleindex is None:
            if self.config.HTcycleindex <= 0:
                self.cycleindex = np.round(kernelscans[1])
            else:
                cycleindex = int(self.config.HTcycleindex)

        # Correct the cycle time to be not a division of the total sequence length but a fixed number of scans
        if cycleindex is not None:
            diff = np.amax(self.fullscans) - self.padindex - cycleindex * len(shortkernel)
            self.padindex += int(diff) + 1
            kernelscans = np.arange(len(shortkernel)) * cycleindex
            self.cycleindex = int(cycleindex)
            print("Correction Diff:", diff)

        # print("Cycle Index:", self.cycleindex)
        self.rollindex = int(-self.shiftindex + self.cycleindex * self.kernelroll)

        # Create the HT kernel
        self.htkernel = np.zeros_like(self.fullscans[int(self.padindex):]).astype(float)
        # print("HT Kernel Length:", len(self.htkernel), "Pad Index:", self.padindex, "Cycle Index:", self.cycleindex,
        #      "Shift Index:", self.shiftindex, "Roll Index:", self.rollindex)

        index = 0
        for i, s in enumerate(self.fullscans):
            # check if this is the first scan above the next scaledkernel
            if s >= kernelscans[index]:
                self.htkernel[i] = shortkernel[index]
                index += 1
                if index >= len(shortkernel):
                    break

        # Smooth the kernel if desired
        if self.config.HTksmooth > 0:
            gausskernel = np.zeros_like(self.htkernel)
            gausskernel += ndis(self.fullscans[self.padindex:], self.padindex, self.config.HTksmooth)
            gausskernel += ndis(self.fullscans[self.padindex:], np.amax(self.fullscans[self.padindex:]) + 1,
                                self.config.HTksmooth)
            gausskernel /= np.sum(gausskernel)
            fftg = fft.rfft(gausskernel)
            self.htkernel = fft.irfft(fft.rfft(self.htkernel) * fftg).real

        # Make fft of kernel for later use
        self.fftk = fft.rfft(self.htkernel).conj()
        self.decontime = self.fulltime
        self.deconscans = self.fullscans[:len(self.decontime)]

    def correct_cycle_time(self):
        """
        Correct the cycle time to be not a division of the total sequence length but a fixed number of scans
        :return: None
        """
        print("Correcting")
        self.get_cycle_time(data)
        self.setup_demultiplex(cycleindex=self.cycleindex)

    def run_demultiplex(self, data, mode=None, chop=True, **kwargs):
        """
        Overarching function to demultiplex data. Calls each demultiplexing type based on config parameter
        :param data: 1D array of data to be deconvolved. Should be same dimension as self.htkernel if HT mode.
        :param mode: Type of demultiplexing. Default is HT.
        :param kwargs: Additional keyword arguments.
        :return: Demultiplexed data. Same length as input.
        """
        if mode is None:
            mode = self.config.demultiplexmode

        if mode == "HT":
            output, data = self.htdecon(data, **kwargs)
        elif mode == "mHT":
            output, data = self.masked_demultiplex_ht(data, **kwargs)
        elif mode == "FT":
            output, data = self.ftdecon(data, aFT=False, **kwargs)
        elif mode == "aFT":
            output, data = self.ftdecon(data, aFT=True, **kwargs)
        else:
            output, data = data
            print("Demultiplexing mode not recognized. Returning original data.", mode)

        if chop:
            if self.config.HToutputub > 0:
                # find index of first peak greater than self.config.HToutputub
                i1 = np.argmax(self.decontime > self.config.HToutputub)
                output = output[:i1]
                self.decontime = self.decontime[:i1]
                self.deconscans = self.deconscans[:i1]
            if self.config.HToutputlb > 0:
                # find index of first peak greater than self.config.HToutputlb
                i2 = np.argmax(self.decontime > self.config.HToutputlb)
                output = output[i2:]
                self.decontime = self.decontime[i2:]
                self.deconscans = self.deconscans[i2:]

        return output, data

    def setup_demultiplex(self, mode=None, **kwargs):
        """
        Overarching function to set up demultiplexing. Calls each demultiplexing type based on config parameter
        :param mode: Type of demultiplexing. Default is HT.
        :param kwargs: Additional keyword arguments.
        :return: None
        """
        if mode is None:
            mode = self.config.demultiplexmode

        if mode == "HT":
            self.dtype = float
            self.setup_ht(**kwargs)
        elif mode == "FT" or mode == "aFT":
            self.dtype = complex
            self.setup_ft(**kwargs)

    def htdecon(self, data, **kwargs):
        """
        Deconvolve the data using the HT kernel. Need to call setup_ht first.
        :param data: 1D array of data to be deconvolved. Should be same dimension as self.htkernel.
        :param kwargs: Keyword arguments. Currently supports "normalize" which normalizes the output to the maximum
        value. Also supports "gsmooth" which smooths the data with a Gaussian filter before deconvolution.
        Also supports "sgsmooth" which smooths the data with a Savitzky-Golay filter before deconvolution.
        :return: Demultiplexed data. Same length as input.
        """
        # Whether to smooth the data before deconvolution
        if "gsmooth" in kwargs:
            data = scipy.ndimage.gaussian_filter1d(data, kwargs["gsmooth"])
            print("Gaussian Smoothing")
        if "sgsmooth" in kwargs:
            data = scipy.signal.savgol_filter(data, kwargs["sgsmooth"], 2)
            print("Savitzky-Golay Smoothing")

        # Set the range of indexes used in the deconvolution
        # Starts at the pad but shift will move it back
        self.indexrange = [self.padindex - self.shiftindex, len(data) - self.shiftindex]
        # print("Index Range:", self.indexrange, "Pad Index:", self.padindex, "Shift Index:", self.shiftindex, "Data Length:", len(data))

        # Do the convolution
        output = fft.irfft(fft.rfft(data[self.indexrange[0]:self.indexrange[1]]) * self.fftk).real
        # If output is odd, add a 0 to the end
        if len(output) < len(data[self.indexrange[0]:self.indexrange[1]]):
            output = np.append(output, 0)
        # print(len(output), len(data[self.indexrange[0]:self.indexrange[1]]))
        # Shift the output back to the original time
        if self.padindex > 0:
            # add zeros back on the front and roll to the correct index
            output = np.roll(np.concatenate((np.zeros(self.padindex), output)), self.rollindex)
        # print(np.shape(data), np.shape(output))
        if "normalize" in kwargs:
            if kwargs["normalize"]:
                output /= np.amax(output)
        # Return demultiplexed data
        return output, data

    def masked_demultiplex_ht(self, data, win=None, n=None, demult=None, mode="rand", **kwargs):
        if win is None:
            win = self.config.HTwin
        if n is None:
            n = self.config.HTmaskn
        win = int(win)
        n = int(n)

        self.config.htseq = ud.HTseqDict[str(self.config.htbit)]
        htseq = np.array([int(s) for s in self.config.htseq])
        if demult is None:
            self.setup_ht(seq=htseq)
            demult, fulltic = self.htdecon(data, **kwargs)
        outdata = []
        for i in range(n):
            seq = np.array(htseq).astype(int)
            # Alter HT seq
            if mode == "con":
                starti = np.random.randint(len(seq))
                endi = starti + win
                if endi > len(seq):
                    seq[starti:] = -1
                    seq[:endi % len(seq)] = -1
                else:
                    seq[starti:endi] = -1
            elif mode == "rand":
                seq[np.random.randint(len(seq), size=win)] = -1
            elif mode == "alt":
                starti = np.random.randint(len(seq))
                endi = starti + win * 2
                if endi > len(seq):
                    seq[starti::2] = -1
                    seq[:endi % len(seq):2] = -1
                else:
                    seq[starti:endi:2] = -1
            # Set window to -1
            self.setup_ht(seq=seq)
            out, _ = self.htdecon(data, **kwargs)
            outdata.append(out)
            # plt.plot(outdata[-1])
        outdata = np.array(outdata)
        # print('Shape of outdata: ', np.shape(outdata))
        # plt.show()
        # Calculate number of sign changes per data point
        negcount = np.sum(outdata < 0, axis=0) + 1
        # avg = np.average(outdata, axis=0)
        # negcount = ud.safedivide(np.std(outdata, axis=0),avg)
        # Divide the demultiplexed data by the negcount plus 1
        wavg = ud.safedivide(demult, negcount.astype(float))
        return wavg, fulltic

    def setup_ft(self, FTstart=None, FTend=None, nzp=None):
        """
        Set up the Fourier Transform Deconvolution. Sets self.decontime with drift time axis.
        :param FTstart: Starting frequence (hz)
        :param FTend: Ending frequence (hz)
        :return: None
        """
        if FTstart is None:
            FTstart = self.config.FTstart
        else:
            self.config.FTstart = FTstart
        if FTend is None:
            FTend = self.config.FTend
        else:
            self.config.FTend = FTend

        x = self.fulltime

        if nzp is None:
            nzp = self.config.HTtimepad

        x_len = int(len(x))
        sweep_time = x[-1]

        if nzp > 0:
            pad_len = int(2 ** math.ceil(math.log2(int(len(x)))) * nzp)
            time_step = (x[1] - x[0])
            x = np.linspace(x[0], x[-1] + time_step * (pad_len - x_len), pad_len)
            # self.fulltime = x
            # self.fullscans = np.arange(len(self.fulltime))

        freq = np.fft.rfftfreq(len(x), d=(x[1] - x[0]) * 60)

        # Create time axis
        sweepRate = (FTend - FTstart) / (sweep_time * 60)
        self.decontime = freq / sweepRate * 1000.
        self.deconscans = self.fullscans[:len(self.decontime)]

    def ftdecon(self, data, flatten=None, apodize=None, aFT=False, normalize=False, nzp=None, keepcomplex=False):
        """
        Perform Fourier Transform Deconvolution
        :param data: 1D data array of y-data only
        :param flatten: Whether to flatten the TIC before demultiplexing to remove low frequency components
        :param apodize: Whether to apodize the data with a Hanning window
        :param aFT: Whether to use Absorption FT mode
        :param normalize: Whether to normalize the output to the maximum value
        :param nzp: Zero padding factor
        :return: 1D array of demultiplexed data. Same length as input.
        """
        if np.amax(data) == 0:
            if keepcomplex:
                return np.zeros(len(self.decontime), dtype=self.dtype), data
            else:
                return data[:len(self.decontime)], data

        if flatten is not None:
            self.config.FTflatten = flatten
        if apodize is not None:
            self.config.FTapodize = apodize

        y = data

        if self.config.HTksmooth > 0:
            if self.config.HTksmooth < 4:
                self.config.HTksmooth = 4
                print("Warning: FT smoothing kernel too small. Setting to 4.")
            if self.config.HTksmooth > len(y) - 1:
                self.config.HTksmooth = 10
                print("Warning: FT smoothing kernel too long. Setting to 10.")
            y = scipy.signal.savgol_filter(y, int(np.round(float(self.config.HTksmooth))), 3)

        if nzp is None:
            nzp = self.config.HTtimepad
        if nzp > 0:
            if self.config.FTapodize == 0:
                print("Warning: Zero padding requires apodization. Setting apodization to True.")
                self.config.FTapodize = 1

        if self.config.FTapodize == 1:
            # create hanning window
            hanning = np.hanning(len(y) * 2)
            y = y * hanning[len(y):]

        if self.config.FTflatten:
            # fit trendline to y and subtract to eliminate low frequency components
            ytrnd = scipy.signal.savgol_filter(y, 15, 3)
            y = y - ytrnd

        original_len = len(y)
        if nzp > 0 and self.config.FTapodize:
            pad_len = int(2 ** math.ceil(math.log2(int(len(y)))) * nzp)
            z = np.zeros(pad_len)
            z[:len(y)] = y
            y = z

        # Fourier Transform
        Y = np.fft.rfft(y)

        if aFT:
            maxindex = np.argmax(np.abs(Y[5:])) + 5
            phase = np.angle(Y[maxindex])
            Y = Y * np.exp(-1j * phase)
            if not keepcomplex:
                Y = np.real(Y)
        else:
            if not keepcomplex:
                Y = np.abs(Y)

        if normalize:
            Y /= np.amax(Y)

        return Y, y[:original_len]

    def set_timepad_index(self, timepad):
        """
        Find the index of the first scan above a timepad
        :param timepad: Time value
        :return: Index of first scan above timepad
        """
        # find first index above timepad in self.fulltimes
        padindex = np.argmax(self.fulltime >= timepad)
        return padindex

    def get_first_peak(self, data, threshold=0.5):
        """
        Get the index of the first peak in the data above a threshold.
        :param data: Data array 1D
        :param threshold: Threshold as ratio of highest peak
        :return: Index of first data point above threshold
        """
        # Find the first peak above threshold
        maxindex = np.argmax(data)
        maxval = np.amax(data)
        if maxval > 0:
            # Find the first peak above threshold
            peakindex = np.argmax(data[maxindex:] < threshold * maxval) + maxindex
        else:
            peakindex = 0
        return peakindex

    def get_cycle_time(self, data=None, cycleindexguess=None, widthguess=110):
        """
        Get the cycle time from the data. This is the time between the first peak and the next peak.
        Uses an autocorrelation and then peak picking on the autocorrelation.
        May need to adjust the peak picking parameters to get it to work right.
        :param data: Input data
        :param cycleindexguess: Guess for the cycle index
        :param widthguess: Guess for peak width in number of scans
        :return: autocorrelation of the data
        """
        if data is None:
            data = self.fulltic
        # autocorrelation of the data
        ac = np.correlate(data, data, mode="same")
        # Get peak after largest peak
        maxindex = np.argmax(ac)
        ac = ac[maxindex:]
        if cycleindexguess is not None:
            maxindex = cycleindexguess
        else:
            maxindex = len(ac) - 2 * widthguess - 1
        # Get first peak
        self.cycleindex = np.argmax(ac[widthguess:widthguess + maxindex]) + widthguess
        timespacing = np.amax(self.fulltime) / len(self.fulltime)
        self.config.HTcycletime = self.cycleindex * timespacing
        print("Cycle Index:", self.cycleindex, "Cycle Time:", self.config.HTcycletime)
        return ac

    '''

        def htdecon_speedy(self, data):
            """
            Deconvolve the data using the HT kernel. Need to call setup_ht first. Currently unused.
            :param data: 1D data array. Should be same dimension as self.htkernel.
            :return: Demultiplexed data. Same length as input.
            """
            # Do the convolution, and only the convolution... :)
            return fft.irfft(fft.rfft(data) * self.fftk).real, data

        def decon_3d_fft(self, array, **kwargs):
            """
            Developed this to see if it would speed things up. It turns out not to. About half as slow. Leaving in for
            legacy reasons and because it's super cool code.
            :param array: 3D array of data to be deconvolved.
                Should be same length as self.htkernel with the other dimensions set by the harray size.
            :param kwargs: Keyword arguments. Currently supports "normalize" which normalizes the output to the maximum
                value.
            :return: Demultiplexed data array. Same length as input array.
            """
            starttime = time.perf_counter()
            # Slice data to appropriate range
            dims = np.shape(array)
            self.indexrange = [self.padindex - self.shiftindex, dims[0] - self.shiftindex]
            data = array[self.indexrange[0]:self.indexrange[1]]
            dims2 = np.shape(data)
            print(dims2)

            # HT Kernel 3D FFT
            # Create 3D array of copies of 1D kernel
            kernel3d = np.broadcast_to(self.htkernel[:, np.newaxis, np.newaxis], dims2)

            print("3D Kernel", np.shape(kernel3d), time.perf_counter() - starttime)

            # FFT of kernel3d
            fftk3d = fft.rfftn(kernel3d).conj()
            print("3D FFT of Kernel", np.shape(fftk3d), time.perf_counter() - starttime)

            data_fft = fft.rfftn(data)
            print("3D FFT of Data", np.shape(data_fft), time.perf_counter() - starttime)

            # Deconvolve
            output = fft.irfftn(data_fft * fftk3d)
            print("Decon", np.shape(output), time.perf_counter() - starttime)
            output = np.real(output)
            if "normalize" in kwargs:
                if kwargs["normalize"]:
                    output /= np.amax(output)
            # Return demultiplexed data
            return output'''


'''
class UniChromHT(HTEng, ChromEngine):
    def __init__(self, *args, **kwargs):
        """
        Initialize the UniChromHT class. This class is used for handling Hadamard Transform (HT) related operations
        on chromatograms.
        :param args: Arguments
        :param kwargs: Keyword Arguments
        """
        super().__init__(*args, **kwargs)
        print("HT Chromatogram Engine")

    def open_file(self, path):
        """
        Open file and set up the time domain.
        :param path: File path
        :return: None
        """
        self.open_chrom(path)
        times = self.get_minmax_times()
        self.config.HTanalysistime = np.amax(times[1])
        self.fullscans -= np.amin(self.fullscans)
        self.scans = np.array(self.fullscans)
        self.parse_file_name(path)
        print("Loaded File:", path)

    def get_eic(self, massrange):
        """
        Get the EIC from the chromatogram.
        :param massrange: Mass range for EIC selection [low, high]
        :return:
        """
        return self.chromdat.get_eic(mass_range=np.array(massrange))

    def eic_ht(self, massrange):
        """
        Get the EIC and run HT on it.
        :param massrange: Mass range for EIC selection [low, high]
        :return: Demultiplexed data output
        """
        eic = self.get_eic(massrange)
        print(eic.shape)
        self.fulltic = eic[:, 1]
        self.fulltime = eic[:, 0]
        self.setup_demultiplex()
        self.htoutput = self.run_demultiplex(self.fulltic)[0]
        return self.htoutput

    def tic_ht(self, correct=False, **kwargs):
        """
        Get the TIC and run HT on it.
        :param correct: Whether to correct the data for the first peak
        :param kwargs: Deconvolution keyword arguments
        :return: Demultiplexed data output
        """
        self.fulltic = self.ticdat[:, 1]
        self.fulltime = self.ticdat[:, 0]
        self.setup_demultiplex()
        self.htoutput = self.run_demultiplex(self.fulltic, correct=correct, **kwargs)[0]
        return self.htoutput
'''


class UniChromCDEng(HTEng, UniDecCD):
    def __init__(self, *args, **kwargs):
        """
        Initialize the UniChromCDEng class. This class is used for handling Hadamard Transform (HT) related operations
        on CDMS data.
        :param args: Arguments
        :param kwargs: Keyword Arguments
        :return: None
        """
        super(UniChromCDEng, self).__init__(*args, **kwargs)
        print("HT-CD-MS Engine")
        self.config.poolflag = 0

        self.hstack = None
        self.fullhstack = None
        self.topfarray = None
        self.topzarray = None
        self.topharray = None
        self.mzaxis = None
        self.zaxis = None
        self.ccsaxis = None
        self.X = None
        self.Y = None
        self.mass = None
        self.fullhstack_ht = None
        self.mstack = None
        self.mstack_ht = None
        self.fullmstack = None
        self.fullmstack_ht = None
        self.mass_tic = None
        self.mass_tic_ht = None
        self.ccsstack_ht = None
        self.mz = None
        self.ztab = None
        self.sarray = None  # Params for swoop selection

    def open_file(self, path, refresh=False):
        """
        Open CDMS file and set up the time domain.
        :param path: Path to file
        :param refresh: Whether to refresh the data. Default False.
        :return: None
        """
        self.clear_arrays()
        self.open_cdms_file(path, refresh=refresh)
        self.parse_file_name(path)

    def clear_arrays(self, massonly=False):
        """
        Clear arrays to reset.
        :param massonly: Whether to only reset mass arrays. Default False.
        :return: None
        """
        if not massonly:
            self.fullhstack = None
            self.fullhstack_ht = None
        self.mstack = None
        self.mstack_ht = None
        self.fullmstack = None
        self.fullmstack_ht = None
        self.mass_tic = None
        self.mass_tic_ht = None
        self.ccsstack_ht = None

    def prep_time_domain(self):
        """
        Prepare the time domain for CDMS data. Creates scans, fullscans, fulltime arrays.
        Need to set self.config.HTanalysistime before calling this function.
        :return: None
        """
        self.scans = np.unique(self.farray[:, 2])
        if self.config.HTmaxscans < 0:
            maxscans = np.amax(self.scans)
            print("Max Scans: ", maxscans)
        else:
            maxscans = self.config.HTmaxscans
            if self.config.CDScanCompress > 1:
                maxscans = int(maxscans / self.config.CDScanCompress)

            if maxscans < np.amax(self.scans):
                print("Error with setting Max Scans, less than total scan number:", maxscans, np.amax(self.scans))
                maxscans = np.amax(self.scans)

        self.fullscans = np.arange(1, maxscans + 1)
        self.fulltime = self.fullscans * self.config.HTanalysistime / np.amax(self.fullscans)
        self.decontime = self.fulltime
        self.deconscans = self.fullscans

    def process_data_scans(self, transform=True):
        """
        Main function for processing CDMS data, includes filtering, converting intensity to charge, histogramming,
        processing the histograms, and transforming (if specified). Transforming is sometimes unnecessary, which is why
        it can be turned off for speed.

        :param transform: Sets whether to transform the data from m/z to mass. Default True.
        :return: None
        """
        starttime = time.perf_counter()
        # self.process_data()
        self.clear_arrays()
        self.prep_time_domain()

        self.topfarray = deepcopy(self.farray)
        self.topzarray = deepcopy(self.zarray)

        # Prepare histogram
        self.prep_hist(mzbins=self.config.mzbins, zbins=self.config.CDzbins)

        # Create a stack with one histogram for each scan
        self.hstack = np.zeros((len(self.scans), self.topharray.shape[0], self.topharray.shape[1]))

        print("Creating Histograms for Each Scan", time.perf_counter() - starttime)
        # Loop through scans and create histograms
        for i, s in enumerate(self.scans):
            # Pull out subset of data for this scan
            b1 = self.topfarray[:, 2] == s

            # Create histogram
            harray = self.histogramLC(x=self.topfarray[b1, 0], y=self.topzarray[b1], w=self.topfarray[b1, 3])
            harray = self.hist_data_prep(harray)

            # Add histogram to stack
            self.hstack[i] = harray

        self.topharray = np.sum(self.hstack, axis=0)

        # Normalize the histogram
        if self.config.datanorm == 1:
            maxval = np.amax(self.topharray)
            if maxval > 0:
                self.topharray /= maxval
        # Reshape the data into a 3 column array
        # self.data.data3 = np.transpose([np.ravel(self.X, order="F"), np.ravel(self.Y, order="F"),
        #                                np.ravel(self.topharray, order="F")])

        # create full hstack to fill any holes in the data for scans with no ions
        self.fullhstack = np.zeros((len(self.fullscans), self.topharray.shape[0], self.topharray.shape[1]))
        for i, s in enumerate(self.scans):
            self.fullhstack[int(s) - 1] = self.hstack[i]

        print("Process Time HT:", time.perf_counter() - starttime)

    def prep_hist(self, mzbins=1, zbins=1, mzrange=None, zrange=None):
        """
        Prepare the histogram for process_data_scans CDMS data.
        :param mzbins: Bin size for m/z
        :param zbins: Bin size for charge
        :param mzrange: m/z range
        :param zrange: charge range
        :return: None
        """
        # Set up parameters
        if mzbins < 0.001:
            print("Error, mzbins too small. Changing to 1", mzbins)
            mzbins = 1
            self.config.mzbins = 1
        self.config.mzbins = mzbins
        self.config.CDzbins = zbins

        x = self.farray[:, 0]
        y = self.zarray
        # Set Up Ranges
        if mzrange is None:
            mzrange = [np.floor(np.amin(x)), np.amax(x)]
        if zrange is None:
            zrange = [np.floor(np.amin(y)), np.amax(y)]

        # Create Axes
        mzaxis = np.arange(mzrange[0] - mzbins / 2., mzrange[1] + mzbins / 2, mzbins)
        # Weird fix to make this axis even is necessary for CuPy fft for some reason...
        if len(mzaxis) % 2 == 1:
            mzaxis = np.arange(mzrange[0] - mzbins / 2., mzrange[1] + 3 * mzbins / 2, mzbins)
        zaxis = np.arange(zrange[0] - zbins / 2., zrange[1] + zbins / 2, zbins)
        self.mzaxis = mzaxis
        self.zaxis = zaxis

        # Create an empty histogram with the right dimensions
        self.topharray, self.mz, self.ztab = np.histogram2d([], [], [self.mzaxis, self.zaxis])
        self.topharray = np.transpose(self.topharray)

        # Calculate the charge and m/z arrays for histogram
        self.mz = self.mz[1:] - self.config.mzbins / 2.
        self.ztab = self.ztab[1:] - self.config.CDzbins / 2.
        self.data.ztab = self.ztab
        # Create matrices of m/z and z
        self.X, self.Y = np.meshgrid(self.mz, self.ztab, indexing='xy')
        # Calculate the mass of every point on histogram
        self.mass = (self.X - self.config.adductmass) * self.Y

    def histogramLC(self, x=None, y=None, w=None):
        """
        Histogram function used for LC-CD-MS data.
        :param x: x-axis (m/z)
        :param y: y-axis (charge)
        :return: Histogram array
        """
        # X is m/z
        if x is None:
            x = self.farray[:, 0]
        # Y is charge
        if y is None:
            y = self.zarray
        # Look for empty data and return zero array if so
        if len(x) <= 1 or len(y) <= 1:
            harray = np.zeros_like(self.topharray)
            return harray

        if self.config.CDiitflag:
            if w is None:
                weights = self.farray[:, 3]
                if len(weights) != len(x):
                    weights = None
            else:
                weights = w
        else:
            weights = None

        # Create histogram
        harray, mz, ztab = np.histogram2d(x, y, [self.mzaxis, self.zaxis], weights=weights)
        # harray = fast_histogram.histogram2d(x, y, [len(self.mzaxis)-1, len(self.zaxis)-1],
        #                                    [(np.amin(self.mzaxis), np.amax(self.mzaxis)),
        #                                     (np.amin(self.zaxis), np.amax(self.zaxis))])
        # Transpose and return
        harray = np.transpose(harray)
        return harray

    def create_chrom(self, farray, **kwargs):
        """
        Create a chromatogram from the farray.
        :param farray: Data array of features
        :param kwargs: Keyword arguments. Currently supports "normalize" which normalizes the output to the maximum.
        :return: TIC/EIC in 2D array (time, intensity)
        """
        # Count of number of time each scans appears in farray
        scans, counts = np.unique(farray[:, 2], return_counts=True)

        fulleic = np.zeros_like(self.fullscans)
        for i, s in enumerate(scans):
            fulleic[int(s) - 1] = counts[i]
        # Normalize
        if "normalize" in kwargs:
            if kwargs["normalize"]:
                fulleic = fulleic / np.amax(fulleic)
        else:
            fulleic /= np.amax(fulleic)
        return np.transpose([self.fulltime, fulleic])

    def get_tic(self, farray=None, **kwargs):
        """
        Get the TIC from the farray.
        :param farray: Optional input feature array
        :param kwargs: Keywords to be passed down to create_chrom
        :return: 2D array of TIC (time, intensity)
        """
        self.prep_time_domain()
        if farray is None:
            farray = self.farray
        fulltic = self.create_chrom(farray, **kwargs)
        self.fulltic = fulltic[:, 1]
        return fulltic

    def tic_ht(self, **kwargs):
        """
        Get the TIC and run HT on it.
        :param kwargs: Keyword arguments. Passed down to create_chrom and htdecon.
        :return: Demultiplexed data output. 2D array (time, intensity)
        """
        self.get_tic(**kwargs)
        self.setup_demultiplex()
        # print(np.shape(self.decontime), np.shape(self.fulltic))
        self.htoutput, self.fulltic = self.run_demultiplex(self.fulltic, **kwargs)
        if self.config.datanorm:
            self.htoutput /= np.amax(self.htoutput)
            self.fulltic /= np.amax(self.fulltic)
        # print(np.shape(self.decontime), np.shape(self.htoutput))
        if len(self.htoutput) != len(self.decontime):
            print("ERROR: Length of HT output does not match time domain.", len(self.htoutput), len(self.decontime))
            raise ValueError("Length of HT output does not match time domain.")
            # self.decontime = np.arange(len(self.htoutput))
            # self.deconscans = np.arange(len(self.htoutput))
        else:
            output = np.transpose([self.decontime, self.htoutput])
        return output

    def extract_swoop_subdata(self, sarray):
        """
        Extract a subdata object based on the Swoop selection.
        :param sarray: Swoop array, m/z mid, z mid, z spread (vertical), z width (horizontal)
        :return: Data Object
        """
        # Calculate the Swoop m/z range, zrange, upper charge, and lower charge bounds
        mz, z, zup, zdown = ud.calc_swoop(sarray, adduct_mass=self.config.adductmass)
        # Create Boolean array
        bsum = np.zeros(self.harray.shape)
        # Loop over all charge states
        for i, zval in enumerate(z):
            # For each charge state, filter z values within the bounds
            b1 = self.ztab >= zdown[i]
            b2 = self.ztab <= zup[i]
            bz = b1 & b2

            # Filter m/z values within the bounds of that charge state
            mzmin, mzmax = ud.get_swoop_mz_minmax(mz, i)
            b1 = self.mz >= mzmin
            b2 = self.mz <= mzmax
            bmz = b1 & b2

            # Take everything that is within the charge and m/z range for that charge state
            # Add rather than multiple because it's OR for each charge state
            bsum += np.outer(bz, bmz)

        # Create new data object
        newd = deepcopy(self.data)
        # Filter Harray with the boolean array
        newh2 = self.harray * bsum
        # Transform and populate the new data object
        if np.sum(newh2) != 0:
            newd = self.transform(newh2, newd)
        # Return data object
        return newd

    def get_swoop_eic(self, sarray, **kwargs):
        """
        Extract an EIC based on the Swoop selection.
        :param sarray: Swoop array, m/z mid, z mid, z spread (vertical), z width (horizontal)
        :param kwargs: Keywords to be passed down to create_chrom
        :return: Data Object
        """
        # Calculate the Swoop m/z range, zrange, upper charge, and lower charge bounds
        mz, z, zup, zdown = ud.calc_swoop(sarray, adduct_mass=self.config.adductmass)
        # Create Boolean array
        bsum = np.zeros(len(self.farray))
        # Loop over all charge states
        for i, zval in enumerate(z):
            # For each charge state, filter z values within the bounds
            b1 = self.zarray >= zdown[i]
            b2 = self.zarray <= zup[i]
            bz = b1 & b2

            # Filter m/z values within the bounds of that charge state
            mzmin, mzmax = ud.get_swoop_mz_minmax(mz, i)
            b1 = self.farray[:, 0] >= mzmin
            b2 = self.farray[:, 0] <= mzmax
            bmz = b1 & b2

            # Take everything that is within the charge and m/z range for that charge state
            # Add rather than multiple because it's OR for each charge state
            bsum += bmz * bz
        # Filter farray
        farray2 = self.farray[bsum.astype(bool)]

        # Create EIC
        eic = self.create_chrom(farray2, **kwargs)
        return eic

    def get_eic(self, mzrange, zrange, **kwargs):
        """
        Get the EIC from the farray.
        :param mzrange: m/z range
        :param zrange: charge range
        :param kwargs: Keywords to be passed down to create_chrom
        :return: 2D array of EIC (time, intensity)
        """
        # Filter farray
        b1 = self.farray[:, 0] >= mzrange[0]
        b2 = self.farray[:, 0] <= mzrange[1]
        b3 = self.zarray >= zrange[0]
        b4 = self.zarray <= zrange[1]
        b = np.logical_and(b1, b2)
        b = np.logical_and(b, b3)
        b = np.logical_and(b, b4)
        farray2 = self.farray[b]

        # Create EIC
        eic = self.create_chrom(farray2, **kwargs)
        return eic

    def extract_subdata(self, mzrange, zrange):
        """
        Extract a subdata object based on m/z and charge range.
        :param mzrange: m/z range, [low, high]
        :param zrange: z range, [low, high]
        :return: Data Object
        """
        b1 = self.mz >= mzrange[0]
        b2 = self.mz <= mzrange[1]
        b3 = self.ztab >= zrange[0]
        b4 = self.ztab <= zrange[1]
        b = np.logical_and(b1, b2)
        b2 = np.logical_and(b3, b4)

        b3 = np.outer(b2, b)
        newd = deepcopy(self.data)

        newh2 = self.harray * b3

        if np.sum(newh2) != 0:
            newd = self.transform(newh2, newd)
        # newd = self.transform(newh, newd, ztab=ztab, mz=mz, mass=mass)

        return newd

    def eic_ht(self, mzrange, zrange, sarray=None, **kwargs):
        """
        Get the EIC and run HT on it.
        :param mzrange: m/z range
        :param zrange: charge range
        :param sarray: Swoop array, m/z mid, z mid, z spread (vertical), z width (horizontal), default None
        :param kwargs: Keyword arguments. Passed down to create_chrom and htdecon.
        :return: Demultiplexed data output. 2D array (time, intensity)
        """
        if sarray is not None and sarray[0] != -1:
            eic = self.get_swoop_eic(sarray, **kwargs)
        else:
            eic = self.get_eic(mzrange, zrange, **kwargs)
        self.setup_demultiplex()
        self.htoutput, eic[:, 1] = self.run_demultiplex(eic[:, 1], **kwargs)

        if len(self.htoutput) != len(self.decontime):
            print("ERROR: Length of HT output does not match time domain.", len(self.htoutput), len(self.decontime))
            raise ValueError("Length of HT output does not match time domain.")
            # self.decontime = np.arange(len(self.htoutput))
            # self.deconscans = np.arange(len(self.htoutput))
        else:
            output = np.transpose([self.decontime, self.htoutput])

        return output, eic

    def run_all_ht(self):
        """
        Run HT on all data in full 3D array. Will call process_data_scans if necessary.
        :return: TIC based on demultiplexed data. 2D array (time, intensity)
        """
        starttime = time.perf_counter()
        if self.fullhstack is None:
            self.process_data_scans()
        self.clear_arrays(massonly=True)

        # Setup HT
        self.setup_demultiplex()

        # Run the HT on each track in the stack
        self.fullhstack_ht = np.empty((len(self.decontime), self.topharray.shape[0], self.topharray.shape[1])
                                      , dtype=self.dtype)

        processed_tic = np.zeros_like(self.fulltime)
        for i in range(len(self.mz)):
            for j in range(len(self.ztab)):
                trace = self.fullhstack[:, j, i]
                tracesum = self.topharray[j, i]
                if tracesum <= self.config.intthresh:
                    htoutput = np.zeros_like(self.decontime)
                else:
                    htoutput, trace = self.run_demultiplex(trace, chop=False, keepcomplex=True)
                    processed_tic += trace
                self.fullhstack_ht[:, j, i] = htoutput

        # Clip all values below 1e-6 to zero
        # self.fullhstack_ht[np.abs(self.fullhstack_ht) < 1e-6] = 0

        if self.config.HToutputub > 0:
            # find index of first peak greater than self.config.HToutputub
            i1 = np.argmax(self.decontime > self.config.HToutputub)
            self.fullhstack_ht = self.fullhstack_ht[:i1]
            self.decontime = self.decontime[:i1]
            self.deconscans = self.deconscans[:i1]
        if self.config.HToutputlb > 0:
            # find index of first peak greater than self.config.HToutputlb
            i2 = np.argmax(self.decontime > self.config.HToutputlb)
            self.fullhstack_ht = self.fullhstack_ht[i2:]
            self.decontime = self.decontime[i2:]
            self.deconscans = self.deconscans[i2:]

        tic = np.sum(self.fullhstack_ht, axis=(1, 2))
        if self.config.demultiplexmode != "FT":
            tic = np.real(tic)
        else:
            tic = np.abs(tic)
        ticdat = np.transpose(np.vstack((self.decontime, tic)))
        if self.config.datanorm == 1:
            ticdat[:, 1] /= np.amax(ticdat[:, 1])
            processed_tic /= np.amax(processed_tic)

        print("Full HT Demultiplexing Done:", time.perf_counter() - starttime)
        return ticdat, processed_tic

    def select_ht_range(self, range=None):
        """
        Select a range of time from HT stack and processes it as a histogram array
        :param range: Time range
        :return: 2D histogram array
        """
        if range is None:
            range = [np.amin(self.decontime), np.amax(self.decontime)]
        b1 = self.decontime >= range[0]
        b2 = self.decontime <= range[1]
        b = np.logical_and(b1, b2)
        substack_ht = np.real(self.fullhstack_ht[b])
        self.harray = np.sum(substack_ht, axis=0)
        self.harray = np.clip(self.harray, 0, np.amax(self.harray))
        self.harray_process()
        return self.harray

    def select_raw_range(self, range=None):
        """
        Select a range of time from raw data and processes it as a histogram array
        :param range: Time range
        :return: 2D histogram array
        """
        if self.fullhstack is None:
            self.process_data_scans()
        if range is None:
            range = [np.amin(self.fulltime), np.amax(self.fulltime)]
        b1 = self.fulltime >= range[0]
        b2 = self.fulltime <= range[1]
        b = np.logical_and(b1, b2)
        substack = self.fullhstack[b]
        self.harray = np.sum(substack, axis=0)
        self.harray = np.clip(self.harray, 0, np.amax(self.harray))
        self.harray_process()
        return self.harray

    def transform_array(self, array, dtype=float):
        """
        Transforms a histogram stack from m/z to mass
        :param array: Histogram stack. Shape is time vs. charge vs. m/z.
        :return: Transformed array
        """
        mlen = len(self.massaxis)
        zlen = len(self.ztab)
        outarray = np.zeros((len(array), mlen, zlen), dtype=dtype)

        for j in range(len(self.ztab)):
            indexes = np.array([ud.nearest(self.massaxis, m) for m in self.mass[j]])
            uindexes = np.unique(indexes)

            subarray = array[:, j]
            # Sum together everything with the same index
            subarray = np.transpose([np.sum(subarray[:, indexes == u], axis=1) for u in uindexes])
            # Add to the output array
            outarray[:, uindexes, j] += subarray

        return np.sum(outarray, axis=2), outarray

    def transform_stacks(self):
        """
        Transform the histogram stacks from m/z to mass. Calls transform_array function on each stack.
        :return: None
        """
        if self.massaxis is None:
            self.process_data(transform=True)
        if self.fullhstack is None:
            self.process_data_scans(transform=True)

        self.ccsstack_ht = None

        if self.config.poolflag == 1:
            print("Transforming Stacks by Interpolation")
        else:
            print("Transforming Stacks by Integration")
        starttime = time.perf_counter()
        mlen = len(self.massaxis)

        self.mstack, self.fullmstack = self.transform_array(self.fullhstack)

        self.mass_tic = np.transpose([self.fulltime, np.sum(self.mstack, axis=1)])

        if self.config.datanorm == 1:
            norm = np.amax(self.mass_tic[:, 1])
            self.mass_tic[:, 1] /= norm
            self.mstack /= norm

        print("Full Mass 1 Transform Done:", time.perf_counter() - starttime)

        if self.fullhstack_ht is None:
            return

        self.mstack_ht, self.fullmstack_ht = self.transform_array(self.fullhstack_ht, dtype=self.dtype)

        self.mass_tic_ht = np.transpose([self.decontime, np.real(np.sum(self.mstack_ht, axis=1))])

        if self.config.datanorm == 1:
            norm = np.amax(self.mass_tic_ht[:, 1])
            self.mass_tic_ht[:, 1] /= norm
            # self.mstack_ht /= norm
        print("Full Mass 2 Transform Done:", time.perf_counter() - starttime)

    def get_mass_eic(self, massrange, zrange=None, ht=False):
        """
        Get the EIC for a mass range after transforming the data to mass. Can be either HT or not.
        :param massrange: Mass range
        :param zrange: Charge range. Default None sets to full range.
        :param ht: Boolean whether to use HT or not
        :return: 2D array of EIC (time, intensity)
        """
        if zrange is None:
            zrange = [np.amin(self.ztab), np.amax(self.ztab)]

        # Filter mstack
        b1 = self.massaxis >= massrange[0]
        b2 = self.massaxis <= massrange[1]
        b = np.logical_and(b1, b2)

        # Filter ztab
        b3 = self.ztab >= zrange[0]
        b4 = self.ztab <= zrange[1]
        b5 = np.logical_and(b4, b3)

        if ht:
            array = self.fullmstack_ht
            xvals = self.decontime
        else:
            array = self.fullmstack
            xvals = self.fulltime

        substack = array[:, b, :][:, :, b5]
        mass_eic = np.sum(substack, axis=(1, 2))
        return np.transpose([xvals, mass_eic])

    def get_ccs_eic(self, massrange=None, zrange=None, mzrange=None, normalize=False):
        """
        Get the EIC for a mass range after transforming the data to CCS.
        :param massrange: Mass range
        :param zrange: Charge range. Default None sets to full range.
        :param mzrange: m/z range. Default None sets to full range.
        :param normalize: Whether to normalize the output to the maximum.
        :return: 2D array of EIC (time, intensity)
        """
        array = self.ccsstack_ht
        print("Ranges:", massrange, zrange, mzrange)
        if zrange is not None:
            # Filter ztab
            b3 = self.ztab >= zrange[0]
            b4 = self.ztab <= zrange[1]
            b5 = np.logical_and(b4, b3)
            # Project b5 into array of same shape
            b5 = np.array([[b5 for i in range(len(self.massaxis))] for j in range(len(self.ccsaxis))])
        else:
            b5 = np.ones_like(array, dtype=bool)

        if mzrange is not None:
            dt3d, mass3d, ztab3d = np.meshgrid(self.ccsaxis, self.massaxis, self.ztab, indexing='ij')
            mz3d = (mass3d + ztab3d * self.config.adductmass) / ztab3d
            b6 = mz3d >= mzrange[0]
            b7 = mz3d <= mzrange[1]
            b8 = np.logical_and(b7, b6)
        else:
            b8 = np.ones_like(array, dtype=bool)

        if massrange is not None:
            # Filter mstack
            b1 = self.massaxis >= massrange[0]
            b2 = self.massaxis <= massrange[1]
            b = np.logical_and(b1, b2)
            b = np.transpose([[b for i in range(len(self.ccsaxis))] for j in range(len(self.ztab))], axes=[1, 2, 0])
        else:
            b = np.ones_like(array, dtype=bool)

        # Merge B8, B5 and B into one large 3D array
        b8 = b8 * b5 * b

        substack = array * b8

        ccs_eic = np.abs(np.sum(substack, axis=(1, 2)))

        if self.config.FTsmooth > 0:
            ccs_eic = scipy.signal.savgol_filter(ccs_eic, int(self.config.FTsmooth), 3)

        if normalize:
            ccs_eic /= np.amax(ccs_eic)

        return np.transpose([self.ccsaxis, ccs_eic])

    def transform_dt_ccs_array(self, array, keepcomplex=True):
        """
        Transforms a histogram stack from drift time to CCS
        :param array: Histogram Mass stack. Shape is time vs. mass vs. charge.
        :return: Transformed array
        """
        # Create dt, mass, and z arrays in 3D
        dt3d, mass3d, ztab3d = np.meshgrid(self.decontime, self.massaxis, self.ztab, indexing='ij')
        # Parallel calc all CCS values into new 3d array
        calc_linear_ccsconst(self.config)
        ccs3d = calc_linear_ccs(mass3d, ztab3d, dt3d, self.config)
        # Create new CCS axis
        minccs = np.amin(ccs3d)
        maxccs = np.amax(ccs3d)
        if self.config.ccsbins == -1:
            binsize = (maxccs - minccs) / len(self.decontime)
        else:
            binsize = self.config.ccsbins

        ccsaxis = np.arange(minccs, maxccs, binsize)
        # Create bins shifted by half a bin
        ccsbins = np.arange(minccs - binsize / 2., maxccs + binsize / 2., binsize)
        self.ccsaxis = ccsaxis

        mlen = len(self.massaxis)
        zlen = len(self.ztab)
        outarray = np.zeros((len(ccsaxis), mlen, zlen), dtype=self.dtype)
        if not keepcomplex:
            array = np.abs(array)
        # Loop through array and paste back onto the new CCS axis
        for i in range(mlen):
            for j in range(zlen):
                # indexes = np.array([ud.nearest(self.ccsaxis, c) for c in ccs3d[:, i, j]])
                # outlist = np.zeros_like(self.ccsaxis, dtype=self.dtype)
                # outlist[indexes] += array[:, i, j]

                outlist = np.histogram(ccs3d[:, i, j], bins=ccsbins, weights=array[:, i, j])[0]
                outarray[:, i, j] = outlist
                # interpolate the existing arrray onto the new ccs axis
                # outarray[:, i, j] = np.interp(self.ccsaxis, ccs3d[:, i, j], array[:, i, j])

        return outarray

    def ccs_transform_stacks(self):
        """
        Transform the histogram stacks from drift time to CCS.
        :return: None
        """
        if self.fullmstack_ht is None:
            self.run_all_ht()
            self.transform_stacks()
        starttime = time.perf_counter()
        self.ccsstack_ht = self.transform_dt_ccs_array(self.fullmstack_ht)
        print("Full CCS Transform Done:", time.perf_counter() - starttime)

        ccs_tic = np.sum(self.ccsstack_ht, axis=(1, 2))
        ccs_tic = np.abs(ccs_tic)
        b1 = ccs_tic > 0
        ccs_tic = ccs_tic[b1]
        self.ccsaxis = self.ccsaxis[b1]
        self.ccsstack_ht = self.ccsstack_ht[b1]
        if self.config.FTsmooth > 0:
            if self.config.FTsmooth > len(ccs_tic):
                self.config.FTsmooth = len(ccs_tic)
            if len(ccs_tic) > 4:
                print("Smoothing", self.config.FTsmooth)
                ccs_tic = scipy.signal.savgol_filter(ccs_tic, int(self.config.FTsmooth), 3)

        ccs_tic = np.transpose(np.vstack((self.ccsaxis, ccs_tic)))
        return ccs_tic

    def convert_trace_to_ccs(self, trace, mzrange, zrange, sarray=None, normalize=False):
        """
        Convert a trace to CCS.
        :param trace: 1D array of intensity
        :param mzrange: m/z range
        :param zrange: charge range
        :param sarray: Swoop array, m/z mid, z mid, z spread (vertical), z width (horizontal), default None
        :param normalize: Whether to normalize the output to the maximum.
        :return: 2D array of CCS (time, intensity)
        """
        # If needed, process the scans
        if self.topharray is None:
            self.process_data_scans()

        # Swoop Array Filtering
        if sarray is not None and sarray[0] != -1:
            b = np.zeros_like(self.topharray)
            # Calculate the Swoop m/z range, zrange, upper charge, and lower charge bounds
            mz, z, zup, zdown = ud.calc_swoop(sarray, adduct_mass=self.config.adductmass)
            # Loop over all charge states
            for i, zval in enumerate(z):
                # For each charge state, filter z values within the bounds
                b1 = self.Y >= zdown[i]
                b2 = self.Y <= zup[i]
                bz = b1 & b2

                # Filter m/z values within the bounds of that charge state
                mzmin, mzmax = ud.get_swoop_mz_minmax(mz, i)
                b1 = self.X >= mzmin
                b2 = self.X <= mzmax
                bmz = b1 & b2

                # Take everything that is within the charge and m/z range for that charge state
                # Add rather than multiple because it's OR for each charge state
                b += bz * bmz
            b = b.astype(bool)
        else:
            # Basic Rectangle Filtering
            b1 = self.X >= mzrange[0]
            b2 = self.X <= mzrange[1]
            b3 = self.Y >= zrange[0]
            b4 = self.Y <= zrange[1]
            b = np.logical_and(b1, b2)
            b = np.logical_and(b, b3)
            b = np.logical_and(b, b4)
        # Filter int, mz, z
        intarray = self.topharray[b]
        mzarray = self.X[b]
        zarray = self.Y[b]

        # Calc avg mz, z, and mass
        avgmz = np.sum(intarray * mzarray) / np.sum(intarray)
        avgz = np.sum(intarray * zarray) / np.sum(intarray)
        avgz = np.round(avgz)
        avgmass = (avgmz - self.config.adductmass) * avgz

        # Convert DT to CCS
        calc_linear_ccsconst(self.config)
        trace = deepcopy(trace)
        ccs_axis = calc_linear_ccs(avgmass, avgz, trace[:, 0], self.config)

        # Process Intensity Data
        if normalize:
            trace[:, 1] /= np.amax(trace[:, 1])
        if self.config.FTsmooth > 0:
            if self.config.FTsmooth > len(ccs_axis):
                self.config.FTsmooth = 10.
            trace[:, 1] = scipy.signal.savgol_filter(trace[:, 1], int(self.config.FTsmooth), 3)

        return np.transpose([ccs_axis, trace[:, 1]]), avgz, avgmz


if __name__ == '__main__':
    eng = UniChromCDEng()
    # eng = UniChromHT()

    dir = "C:\\Data\\HT-CD-MS"
    dir = "Z:\\Group Share\\Skippy\\Projects\\HT\\Example data for MTM\\2023-10-26"
    # dir = "Z:\\Group Share\\Skippy\\Projects\\HT\\Example data for MTM"
    os.chdir(dir)
    path = "C:\\Data\\HT-CD-MS\\20230906 JDS BSA SEC f22 10x dilute STORI high flow 1_20230906171314_2023-09-07-01-43-26.dmt"
    path = "C:\\Data\\HT-CD-MS\\20230906 JDS BSA SEC f22 10x dilute STORI high flow 1_20230906171314.raw"
    path = "2023103 JDS BSA inj5s cyc2m bit3 zp3 rep2.raw"
    path = "2023103 JDS BSA inj5s cyc2m bit3 zp5 rep1.raw"
    path = "20231026 JDS BSA cyc2s inj5s bit3 zp1 rep1.raw"
    # path = '20231102_ADH_BSA SEC 15k cyc2m inj3s bit3 zp3 no1.raw'
    # path = '20231202 JDS 0o1uMBgal 0o4uMgroEL shortCol 300ul_m 6_1spl bit5 zp7 inj4s cyc1m AICoff IIT100.RAW'
    # path = "Z:\\Group Share\\Skippy\Projects\HT\\2023-10-13 BSA ADH\\20231013 BSA STORI inj2s cyc2m 5bit_2023-10-13-05-06-03.dmt"
    # path = "20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56.dmt"
    # path = "20231202 JDS 0o1uMBgal 0o4uMgroEL shortCol 300ul_m 6_1spl bit3 zp4 inj5s cyc1m AICoff IIT200.RAW"
    path = "Z:\\Group Share\\Skippy\\Projects\\FT IM CD MS\\GDH\\01302024_GDH_stepsize3_repeat15_5to500_2024-02-06-04-44-16.dmt"
    pathft = "Z:\\Group Share\\Skippy\\Projects\\FT IM CD MS\\01302024_GDH_stepsize3_repeat15_5to500.dmt"
    eng.open_file(pathft)
    eng.process_data()
    eng.process_data_scans()
    eng.run_all_ht()
    eng.transform_stacks()
    eng.config.ccsbins = 1
    eng.ccs_transform_stacks()

    ccs_tic = np.sum(eng.ccsstack_ht, axis=(1, 2))
    plt.plot(eng.ccsaxis, ccs_tic)
    plt.show()
    # eng.get_mass_eic([8500, 10500], [1, 100])
    # eng.process_data_scans()
    # xicdata = eng.get_eic([8500, 10500])
    # eng.eic_ht([8500, 10500], [1, 100])
    # eng.tic_ht()
    # np.savetxt("tic.txt", np.transpose([eng.fulltime, eng.fulltic]))
    # ac = eng.get_cycle_time(eng.fulltic)

    # import matplotlib.pyplot as plt
    # plt.plot(ac)
    # plt.show()
