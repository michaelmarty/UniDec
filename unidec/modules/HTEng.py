import scipy.ndimage
from modules.fitting import *

from modules.CDEng import *
from modules.ChromEng import *
import unidec.tools as ud
from modules.unidecstructure import DataContainer

import matplotlib.pyplot as plt


class HTEng:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.config.htmode = True

        # Empty Arrays
        self.scans = []
        self.fullscans = []
        self.fulltime = []
        self.fulltic = []
        self.htkernel = []
        self.fftk = []
        self.htoutput = []

        # Index Values
        self.padindex = 0
        self.shiftindex = 0
        self.kernelroll = 0
        self.cycleindex = 0
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

        print("HT Engine")

    def parse_file_name(self, path):
        if "cyc" in path:
            # Find location of cyc and take next character
            cycindex = path.find("cyc")
            cyc = int(path[cycindex + 3])
            self.config.HTcycletime = float(cyc)

        if "zp" in path:
            # Find location of zp and take next character
            zpindex = path.find("zp")
            zp = float(path[zpindex + 2])
            self.config.HTtimepad = zp * self.config.HTcycletime

        if "bit" in path:
            # Find location of bit and take next character
            bitindex = path.find("bit")
            self.config.htbit = int(path[bitindex + 3])
            try:
                self.config.htseq = ud.HTseqDict[str(self.config.htbit)]
            except KeyError as e:
                print("ERROR: Bit not found in dictionary. Using 3 bit sequence.")
                print(e)
                self.config.htseq = '1110100'
                self.config.htbit = 3

        print("Cycle Time:", self.config.HTcycletime, "Time Pad:", self.config.HTtimepad,
              "HT Sequence:", self.config.htseq)

    def setup_ht(self, cycleindex=None):
        """
        Sets up the HT kernel for deconvolution. This is a binary sequence that is convolved with the data to
        deconvolve the data. The sequence is defined by the htseq variable. The timepad and timeshift variables
        shift the sequence in time. The config.HTksmooth variable smooths the kernel to reduce ringing.
        :param cycleindex: The length of a cycle in number of scans.
            If not specified, default is the full number of scans divided by the number of cycles.
        :return: None
        """
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
        shortkernel = seqarray

        # Roll short kernel so that the first 1 is in index 0
        self.kernelroll = -np.argmax(shortkernel)
        shortkernel = np.roll(shortkernel, self.kernelroll)

        scaledkernel = np.arange(len(shortkernel)) / (len(shortkernel))
        kernelscans = scaledkernel * (np.amax(self.scans) - self.padindex)

        if cycleindex is None:
            if self.config.HTcycleindex <= 0:
                self.cycleindex = np.round(kernelscans[1])
            else:
                cycleindex = int(self.config.HTcycleindex)

        # Correct the cycle time to be not a division of the total sequence length but a fixed number of scans
        if cycleindex is not None:
            diff = np.amax(self.scans) - self.padindex - cycleindex * len(shortkernel)
            self.padindex += int(diff) + 1
            kernelscans = np.arange(len(shortkernel)) * cycleindex
            self.cycleindex = int(cycleindex)
            print("Correction Diff:", diff)

        # Create the HT kernel
        self.htkernel = np.zeros_like(self.fullscans[int(self.padindex):]).astype(float)
        print("HT Kernel Length:", len(self.htkernel), "Pad Index:", self.padindex, "Cycle Index:", self.cycleindex)

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
            self.gausskernel = np.zeros_like(self.htkernel)
            self.gausskernel += ndis(self.fullscans[self.padindex:], self.padindex, self.config.HTksmooth)
            self.gausskernel += ndis(self.fullscans[self.padindex:], np.amax(self.fullscans[self.padindex:]) + 1,
                                     self.config.HTksmooth)
            self.gausskernel /= np.sum(self.gausskernel)
            self.fftg = np.fft.fft(self.gausskernel)
            self.htkernel = np.fft.ifft(np.fft.fft(self.htkernel) * self.fftg)

        # Make fft of kernel for later use
        self.fftk = np.fft.fft(self.htkernel).conj()

    def htdecon(self, data, *args, **kwargs):
        # Whether to correct the kernel to optimize the cycle time
        if "correct" in kwargs:
            if kwargs["correct"]:
                print("Correcting")
                self.get_cycle_time(data)
                self.setup_ht(cycleindex=self.cycleindex)
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
        # print("Index Range:", self.indexrange, "Pad Index:", self.padindex, "Shift Index:", self.shiftindex, "Len:",
        #      len(data))
        # Do the convolution
        output = np.fft.ifft(
            np.fft.fft(data[self.indexrange[0]:self.indexrange[1]]) * self.fftk)

        # Shift the output back to the original time
        if self.padindex > 0:
            # add zeros back on the front
            rollindex = int(-self.shiftindex + self.cycleindex * self.kernelroll)
            output = np.roll(np.concatenate((np.zeros(self.padindex), output)), rollindex)

        output = np.real(output)
        if "normalize" in kwargs:
            if kwargs["normalize"]:
                output /= np.amax(output)
        # Return demultiplexed data
        return output

    def set_timepad_index(self, timepad):
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


class UniChromHT(HTEng, ChromEngine):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        print("HT Chromatogram Engine")

    def open_file(self, path):
        self.open_chrom(path)
        times = self.get_minmax_times()
        self.config.HTanalysistime = np.amax(times[1])
        self.fullscans -= np.amin(self.fullscans)
        self.scans = np.array(self.fullscans)
        self.parse_file_name(path)
        print("Loaded File:", path)

    def eic_ht(self, massrange):
        eic = self.chromdat.get_eic(mass_range=np.array(massrange))
        print(eic.shape)
        self.fulltic = eic[:, 1]
        self.fulltime = eic[:, 0]
        self.setup_ht()
        self.htoutput = self.htdecon(self.fulltic)
        return self.htoutput

    def tic_ht(self, correct=False, *args, **kwargs):
        self.fulltic = self.ticdat[:, 1]
        self.fulltime = self.ticdat[:, 0]
        self.setup_ht()
        self.htoutput = self.htdecon(self.fulltic, correct=correct, *args, **kwargs)
        return self.htoutput


class UniDecCDHT(HTEng, UniDecCD):
    def __init__(self, *args, **kwargs):
        super(UniDecCDHT, self).__init__(*args, **kwargs)
        print("HT-CD-MS Engine")
        self.hstack = None
        self.fullhstack = None
        self.topfarray = None
        self.topzarray = None
        self.topharray = None
        self.mzaxis = None
        self.zaxis = None
        self.X = None
        self.Y = None
        self.mass = None
        self.fullstack_ht = None
        self.mz = None
        self.ztab = None

    def open_file(self, path, refresh=False):
        self.open_cdms_file(path, refresh=refresh)
        self.parse_file_name(path)

    def prep_time_domain(self):
        self.scans = np.unique(self.farray[:, 2])
        self.fullscans = np.arange(1, np.amax(self.scans) + 1)
        self.fulltime = self.fullscans * self.config.HTanalysistime / np.amax(self.fullscans)

    def process_data_scans(self, transform=True):
        """
        Main function for processing CDMS data, includes filtering, converting intensity to charge, histogramming,
        processing the histograms, and transforming (if specified). Transforming is sometimes unnecessary, which is why
        it can be turned off for speed.

        :param transform: Sets whether to transform the data from m/z to mass. Default True.
        :return: None
        """
        starttime = time.perf_counter()
        self.process_data()

        self.prep_time_domain()

        self.topfarray = deepcopy(self.farray)
        self.topzarray = deepcopy(self.zarray)

        # Prepare histogram
        self.prep_hist(mzbins=self.config.mzbins, zbins=self.config.CDzbins)

        # Create a stack with one histogram for each scan
        self.hstack = np.zeros((len(self.scans), self.topharray.shape[0], self.topharray.shape[1]))

        print("Creating Histograms for Each Scan")
        # Loop through scans and create histograms
        for i, s in enumerate(self.scans):
            # Pull out subset of data for this scan
            b1 = self.topfarray[:, 2] == s

            # Create histogram
            harray = self.histogramLC(x=self.topfarray[b1, 0], y=self.topzarray[b1])

            # Transform histogram m/z to mass
            # NEED TO WORK ON THIS SECTION
            if len(harray) > 0 and np.amax(harray) > 0 and False:
                print("TRANSFORMING")
                harray = self.hist_data_prep(harray)
                dataobj = DataContainer()
                if transform:
                    dataobj = self.transform(harray=harray, dataobj=dataobj)

            # Add histogram to stack
            self.topharray += harray
            self.hstack[i] = harray

        # Normalize the histogram
        if self.config.datanorm == 1:
            maxval = np.amax(self.topharray)
            if maxval > 0:
                self.topharray /= maxval

        # Reshape the data into a 3 column array
        self.data.data3 = np.transpose([np.ravel(self.X, order="F"), np.ravel(self.Y, order="F"),
                                        np.ravel(self.topharray, order="F")])

        print("Process Time:", time.perf_counter() - starttime)

    def prep_hist(self, mzbins=1, zbins=1, mzrange=None, zrange=None):
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

    def histogramLC(self, x=None, y=None):
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
        # Create histogram
        harray, mz, ztab = np.histogram2d(x, y, [self.mzaxis, self.zaxis])
        # Transpose and return
        harray = np.transpose(harray)
        return harray

    def create_chrom(self, farray, *args, **kwargs):
        # Count of number of time each scans appears in farray
        scans, counts = np.unique(farray[:, 2], return_counts=True)

        fulleic = np.zeros_like(self.fullscans)
        for i, s in enumerate(scans):
            fulleic[int(s) - 1] = counts[i]
        # Normalize
        if "normalize" in kwargs:
            if kwargs["normalize"]:
                fulleic /= np.amax(fulleic)
        else:
            fulleic /= np.amax(fulleic)
        return np.transpose([self.fulltime, fulleic])

    def get_tic(self, farray=None, *args, **kwargs):
        self.prep_time_domain()
        if farray is None:
            farray = self.farray
        fulltic = self.create_chrom(farray, *args, **kwargs)
        self.fulltic = fulltic[:, 1]
        return fulltic

    def tic_ht(self, *args, **kwargs):
        self.get_tic(*args, **kwargs)
        self.setup_ht()
        self.htoutput = self.htdecon(self.fulltic, *args, **kwargs)
        return np.transpose([self.fulltime, self.htoutput])

    def get_eic(self, mzrange, zrange, *args, **kwargs):
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
        eic = self.create_chrom(farray2, *args, **kwargs)
        return eic

    def eic_ht(self, mzrange, zrange, *args, **kwargs):
        eic = self.get_eic(mzrange, zrange, *args, **kwargs)
        self.setup_ht()
        self.htoutput = self.htdecon(eic[:, 1], *args, **kwargs)
        return np.transpose([self.fulltime, self.htoutput]), eic

    def run_ht(self):
        # create full hstack to fill any holes in the data for scans with no ions
        self.fullhstack = np.zeros((len(self.fullscans), self.topharray.shape[0], self.topharray.shape[1]))
        self.fullhstack_ht = np.zeros((len(self.fullscans), self.topharray.shape[0], self.topharray.shape[1]))
        for i, s in enumerate(self.scans):
            self.fullhstack[int(s) - 1] = self.hstack[i]

        # Run the HT on each track in the stack
        # TEST THIS!
        self.setup_ht()
        for i, x in enumerate(self.mz):
            for j, y in enumerate(self.ztab):
                trace = self.fullhstack[:, j, i]
                htoutput = self.htdecon(trace)
                self.fullhstack_ht[:, j, i] = htoutput


if __name__ == '__main__':

    eng = UniDecCDHT()

    dir = "C:\Data\HT-CD-MS"
    # dir = "Z:\\Group Share\Skippy\Projects\HT\Example data for MTM\\2023-10-26"
    dir = "Z:\\Group Share\\Skippy\\Projects\\HT\\Example data for MTM"
    os.chdir(dir)
    path = "C:\\Data\\HT-CD-MS\\20230906 JDS BSA SEC f22 10x dilute STORI high flow 1_20230906171314_2023-09-07-01-43-26.dmt"
    path = "C:\\Data\\HT-CD-MS\\20230906 JDS BSA SEC f22 10x dilute STORI high flow 1_20230906171314.raw"
    path = "2023103 JDS BSA inj5s cyc2m bit3 zp3 rep2.raw"
    path = "2023103 JDS BSA inj5s cyc2m bit3 zp5 rep1.raw"
    path = "20231026 JDS BSA cyc2s inj5s bit3 zp1 rep1.raw"
    # path = '20231102_ADH_BSA SEC 15k cyc2m inj3s bit3 zp3 no1.raw'
    path = '20231202 JDS 0o1uMBgal 0o4uMgroEL shortCol 300ul_m 6_1spl bit5 zp7 inj4s cyc1m AICoff IIT100.RAW'
    path = "Z:\\Group Share\\Skippy\Projects\HT\\2023-10-13 BSA ADH\\20231013 BSA STORI inj2s cyc2m 5bit_2023-10-13-05-06-03.dmt"
    path = "20231202 JDS Bgal groEL bit5 zp7 inj4s cyc1m_2023-12-07-03-46-56.dmt"
    # path = "20231202 JDS 0o1uMBgal 0o4uMgroEL shortCol 300ul_m 6_1spl bit3 zp4 inj5s cyc1m AICoff IIT200.RAW"

    eng.open_file(path)
    eng.process_data_scans()
    # eng.eic_ht([8500, 10500])
    eng.tic_ht()
    # np.savetxt("tic.txt", np.transpose([eng.fulltime, eng.fulltic]))
    ac = eng.get_cycle_time(eng.fulltic)

    # import matplotlib.pyplot as plt
    # plt.plot(ac)
    # plt.show()
    # exit()
    plt.plot(eng.fulltime, eng.fulltic / np.amax(eng.fulltic))
    plt.plot(eng.fulltime[eng.padindex:], np.roll(eng.htkernel, 0))
    plt.plot(eng.fulltime, eng.htoutput / np.amax(eng.htoutput) - 1)
    plt.show()
    exit()

    eng.run_ht()
    print(np.shape(eng.hstack))
    plt.figure()

    plt.subplot(121)
    for i, x in enumerate(eng.mz):
        for j, y in enumerate(eng.ztab):
            plt.plot(eng.fullscans, eng.fullhstack_ht[:, j, i])

    plt.subplot(122)
    plt.imshow(np.sum(eng.fullhstack[:150], axis=0), aspect="auto", origin="lower",
               extent=[eng.mz[0], eng.mz[-1], eng.ztab[0], eng.ztab[-1]])

    plt.show()

    from unidec.modules import PlotAnimations as PA
    import wx

    app = wx.App(False)
    PA.AnimationWindow(None, eng.fullhstack[50:], mode="2D")
    app.MainLoop()
