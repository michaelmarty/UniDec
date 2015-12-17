import os
import numpy as np
import unidectools as ud
import platform
import matplotlib.cm as cm
from matplotlib.pyplot import colormaps

__author__ = 'Michael.Marty'


# noinspection PyAttributeOutsideInit
class UniDecConfig:
    """
    Class containing all options and configurations for UniDec GUI and Program. Contains methods to export and import
    config to text file for running UniDec core binaries and for storing parameters for GUI.
    """

    def __init__(self):
        """
        Initialize Everything. Set default paths and run self.initialize
        :return: UniDecConfig object
        """
        self.infname = "input.dat"
        self.outfname = "output"
        self.mfile = "mass.dat"
        self.manualfile = "man.dat"
        self.confname = "conf.dat"
        self.ofile = "ofile.dat"
        self.matchfile = "match.csv"
        self.peaksfile = "peaks.dat"
        self.dirname = ''
        self.filename = ''
        self.extension = ""
        self.imflag = 0
        self.initialize()

    def initialize(self):
        """
        Initialize configuration parameters but not paths. Runs self.default_colormaps
        :return: None
        """
        # plotting
        self.publicationmode = 0
        self.discreteplot = 0
        self.cmap = "spectral"
        self.peakcmap = "rainbow"
        self.rawflag = 0

        # data prep
        self.detectoreffva = 0
        self.mzbins = 1
        self.smooth = 0
        self.subbuff = 0
        self.subtype = 0
        self.intthresh = 0
        self.minmz = ''
        self.maxmz = ''

        # UniDec
        self.numit = 100
        self.zzsig = 1
        self.startz = 1
        self.endz = 100
        self.numz = 100
        self.mzsig = 20
        self.psfun = 0
        self.massub = 5000000
        self.masslb = 100
        self.msig = 0
        self.molig = 0
        self.massbins = 100
        self.adductmass = 1.007276467
        self.damp = 0
        self.aggressiveflag = 0
        self.suppression = 0
        self.isotopemode = 0

        # Peak Selection and plotting
        self.peakwindow = 500
        self.peakthresh = 0.1
        self.peakplotthresh = 0.1
        self.separation = 0.025
        self.peaknorm = 1

        # Results
        self.error = 0

        # Other
        self.mtabsig = 0
        self.poolflag = 1
        self.nativezub = 100
        self.nativezlb = -100
        self.inflate = 1
        self.linflag = 0
        self.integratelb = ""
        self.integrateub = ""

        # IM specific
        self.mindt = ''
        self.maxdt = ''
        self.smoothdt = 0
        self.subbufdt = 0
        self.ccslb = 100
        self.ccsub = 25000
        self.nativeccsub = 20000
        self.nativeccslb = -20000
        self.dtsig = 0.2
        self.ccsbins = 100
        self.csig = 0
        self.pusher = 0
        self.zout = 0

        self.temp = 25
        self.pressure = 2
        self.volt = 50
        self.to = 0
        self.driftlength = 0.18202
        self.tcal1 = 0.3293
        self.tcal2 = 6.3597
        self.edc = 1.57
        self.gasmass = 4.002602
        self.twaveflag = 0

        # Misc
        self.batchflag = 0
        self.procflag = 0
        self.runtime = 0
        self.massdatnormtop = 0

        self.mfileflag = False
        self.manualfileflag = False

        self.kendrickmass = None

        self.masslist = []
        self.matchlist = []
        self.oligomerlist = []
        self.manuallist = []
        self.zoffs = []
        self.massoffset = 0
        self.extractshape = 0
        self.gridparams = None
        self.griddecon = None
        self.matchtolerance = 1000

        self.default_colormaps()

    def default_colormaps(self):
        """
        Get default matplotlib colormaps and set names to self.cmaps.
        A more list is set as self.cmaps2. For some reason, matplotlib v.1.5.0 allows viridis and other new colormaps
        with the colormaps() but not cm.datad. Even weirder, they work for some things and not others. Solution, two
        lists, one for the 2D plots and the other for the linear colors.
        :return: None
        """
        m = [i for i in cm.datad]
        self.cmaps = sorted(m)

        m2 = colormaps()
        self.cmaps2 = sorted(m2)

    def config_export(self, name):
        """
        Writes config to file give in name. Typically in format: name value.

        Also exports manuallist, masslist, and oligomerlist to text files.
        :param name: File name to write to.
        :return: None
        """
        self.numz = self.endz - self.startz + 1
        f = open(name, 'w')
        f.write("input " + str(self.infname) + "\n")
        f.write("output " + str(self.outfname) + "\n")
        f.write("numit " + str(self.numit) + "\n")
        f.write("numz " + str(self.numz) + "\n")
        f.write("endz " + str(self.endz) + "\n")
        f.write("startz " + str(self.startz) + "\n")
        f.write("zzsig " + str(self.zzsig) + "\n")
        f.write("mzsig " + str(self.mzsig) + "\n")
        f.write("psfun " + str(self.psfun) + "\n")
        f.write("discreteplot " + str(self.discreteplot) + "\n")
        f.write("massub " + str(self.massub) + "\n")
        f.write("masslb " + str(self.masslb) + "\n")
        f.write("msig " + str(self.msig) + "\n")
        f.write("molig " + str(self.molig) + "\n")
        f.write("massbins " + str(self.massbins) + "\n")
        f.write("mtabsig " + str(self.mtabsig) + "\n")
        if self.mfileflag:
            if not ud.isempty(self.masslist):
                f.write("mfile " + str(self.mfile) + "\n")
            else:
                print "Need to specify mass list. Running without mass list."
        if self.manualfileflag:
            if not ud.isempty(self.manuallist):
                f.write("manualfile " + str(self.manualfile) + "\n")
            else:
                print "Need to specify manual assignments. Running without assignments."
        f.write("minmz " + str(self.minmz) + "\n")
        f.write("maxmz " + str(self.maxmz) + "\n")
        f.write("subbuff " + str(self.subbuff) + "\n")
        f.write("subtype " + str(self.subtype) + "\n")
        f.write("smooth " + str(self.smooth) + "\n")
        f.write("mzbins " + str(self.mzbins) + "\n")
        f.write("peakwindow " + str(self.peakwindow) + "\n")
        f.write("peakthresh " + str(self.peakthresh) + "\n")
        f.write("peakplotthresh " + str(self.peakplotthresh) + "\n")
        f.write("plotsep " + str(self.separation) + "\n")
        f.write("intthresh " + str(self.intthresh) + "\n")
        if self.linflag != 2:
            f.write("speedy " + str(1) + "\n")
        else:
            f.write("speedy " + str(0) + "\n")
        f.write("aggressive " + str(self.aggressiveflag) + "\n")
        f.write("rawflag " + str(self.rawflag) + "\n")
        f.write("adductmass " + str(self.adductmass) + "\n")
        f.write("nativezub " + str(self.nativezub) + "\n")
        f.write("nativezlb " + str(self.nativezlb) + "\n")
        f.write("poolflag " + str(self.poolflag) + "\n")
        f.write("accvol " + str(self.detectoreffva) + "\n")
        f.write("peakshapeinflate " + str(self.inflate) + "\n")
        f.write("damp " + str(self.damp) + "\n")
        f.write("linflag " + str(self.linflag) + "\n")
        f.write("cmap " + str(self.cmap) + "\n")
        f.write("peakcmap " + str(self.peakcmap) + "\n")
        f.write("publicationmode " + str(self.publicationmode) + "\n")
        f.write("isotopemode " + str(self.isotopemode) + "\n")
        f.write("peaknorm " + str(self.peaknorm) + "\n")
        f.write("suppression " + str(self.suppression) + "\n")
        if self.integratelb != "" and self.integrateub != "":
            try:
                f.write("integratelb " + str(self.integratelb) + "\n")
                f.write("integrateub " + str(self.integrateub) + "\n")
            except ValueError:
                print "Failed to write integation areas:", self.integratelb, self.integrateub
                pass

        if self.mindt != '' or self.maxdt != '':
            f.write("zout " + str(self.zout) + "\n")
            f.write("pusher " + str(self.pusher) + "\n")
            f.write("mindt " + str(self.mindt) + "\n")
            f.write("maxdt " + str(self.maxdt) + "\n")
            f.write("ccsub " + str(self.ccsub) + "\n")
            f.write("ccslb " + str(self.ccslb) + "\n")
            f.write("dtsig " + str(self.dtsig) + "\n")
            f.write("csig " + str(self.csig) + "\n")
            f.write("ccsbins " + str(self.ccsbins) + "\n")
            f.write("subbufdt " + str(self.subbufdt) + "\n")
            f.write("smoothdt " + str(self.smoothdt) + "\n")
            f.write("ubnativeccs " + str(self.nativeccsub) + "\n")
            f.write("lbnativeccs " + str(self.nativeccslb) + "\n")
            f.write("twaveflag " + str(self.twaveflag) + "\n")
            if self.twaveflag == 0:  # and (self.pressure!=0 or self.temp!=0 or self.volt!=0 or self.to!=0):
                f.write("temp " + str(self.temp) + "\n")
                f.write("pressure " + str(self.pressure) + "\n")
                f.write("volt " + str(self.volt) + "\n")
                f.write("gasmass " + str(self.gasmass) + "\n")
                f.write("tnaught " + str(self.to) + "\n")
                f.write("driftlength " + str(self.driftlength) + "\n")

            if self.twaveflag == 1:  # and (self.tcal1!=0 or self.tcal2!=0 or self.edc!=0):
                f.write("tcal1 " + str(self.tcal1) + "\n")
                f.write("tcal2 " + str(self.tcal2) + "\n")
                f.write("edc " + str(self.edc) + "\n")
                f.write("gasmass " + str(self.gasmass) + "\n")

        f.close()
        if not ud.isempty(self.masslist):
            ud.dataexport(self.masslist, self.mfile)
        if not ud.isempty(self.manuallist):
            # print self.manuallist
            if self.imflag == 0:
                if self.manuallist.shape[1] == 3:
                    ud.dataexport(self.manuallist, self.manualfile)
                else:
                    print "Manual List Shape is wrong. Try using manual list tool again."
                    print self.manuallist.shape
            else:
                if self.manuallist.shape[1] == 5:
                    ud.dataexport(self.manuallist, self.manualfile)
                else:
                    print "Manual List Shape is wrong. Try using manual list tool again."
                    print self.manuallist.shape
        if not ud.isempty(self.oligomerlist):
            np.savetxt(self.ofile, self.oligomerlist, fmt='%s')

    def config_import(self, name):
        """
        Imports configuration from txt file. Also imports masslist, manuallist, and oligomerlist.
        :param name: File to import from.
        :return: None
        """
        if self.batchflag != 1:
            f = open(name, 'r')
            self.manualfileflag = 0
            self.mfileflag = 0
            for line in f:
                if len(line.split()) > 1:
                    if line.startswith("minmz"):
                        self.minmz = ud.string_to_value(line.split()[1])
                    if line.startswith("maxmz"):
                        self.maxmz = ud.string_to_value(line.split()[1])
                    if line.startswith("mindt"):
                        self.mindt = ud.string_to_value(line.split()[1])
                    if line.startswith("maxdt"):
                        self.maxdt = ud.string_to_value(line.split()[1])
                    if self.batchflag == 0:
                        '''
                        if (line.startswith("input")):
                            self.infname = line.split()[1]
                        if (line.startswith("output")):
                            self.outfname = line.split()[1]
                        '''
                        if line.startswith("numit"):
                            self.numit = ud.string_to_int(line.split()[1])
                        if line.startswith("numz"):
                            self.numz = ud.string_to_int(line.split()[1])
                        if line.startswith("endz"):
                            self.endz = ud.string_to_int(line.split()[1])
                        if line.startswith("startz"):
                            self.startz = ud.string_to_int(line.split()[1])
                        if line.startswith("zzsig"):
                            self.zzsig = ud.string_to_value(line.split()[1])
                        if line.startswith("mzsig"):
                            self.mzsig = ud.string_to_value(line.split()[1])
                        if line.startswith("psfun"):
                            self.psfun = ud.string_to_int(line.split()[1])
                        if line.startswith("discreteplot"):
                            self.discreteplot = ud.string_to_int(line.split()[1])
                        if line.startswith("massub"):
                            self.massub = ud.string_to_value(line.split()[1])
                        if line.startswith("masslb"):
                            self.masslb = ud.string_to_value(line.split()[1])
                        if line.startswith("msig"):
                            self.msig = ud.string_to_value(line.split()[1])
                        if line.startswith("molig"):
                            self.molig = ud.string_to_value(line.split()[1])
                        if line.startswith("massbins"):
                            self.massbins = ud.string_to_value(line.split()[1])
                        if line.startswith("mtabsig"):
                            self.mtabsig = ud.string_to_value(line.split()[1])
                        if line.startswith("mfile"):
                            # self.mfile = line.split()[1]
                            self.mfileflag = True
                        if line.startswith("manualfile"):
                            # self.manualfile = line.split()[1]
                            self.manualfileflag = True
                        if line.startswith("subbuff"):
                            self.subbuff = ud.string_to_value(line.split()[1])
                            if self.subbuff < 0:
                                self.subtype = 2
                                self.subbuff = abs(self.subbuff)
                        if line.startswith("subtype"):
                            self.subtype = ud.string_to_int(line.split()[1])
                        if line.startswith("smooth") and not line.startswith("smoothdt"):
                            self.smooth = ud.string_to_value(line.split()[1])
                        if line.startswith("mzbins"):
                            self.mzbins = ud.string_to_value(line.split()[1])
                        if line.startswith("peakwindow"):
                            self.peakwindow = ud.string_to_value(line.split()[1])
                        if line.startswith("peakthresh"):
                            self.peakthresh = ud.string_to_value(line.split()[1])
                        if line.startswith("peakplotthresh"):
                            self.peakplotthresh = ud.string_to_value(line.split()[1])
                        if line.startswith("plotsep"):
                            self.separation = ud.string_to_value(line.split()[1])
                        if line.startswith("intthresh"):
                            self.intthresh = ud.string_to_value(line.split()[1])
                        if line.startswith("aggressive"):
                            self.aggressiveflag = ud.string_to_int(line.split()[1])
                        if line.startswith("rawflag"):
                            self.rawflag = ud.string_to_int(line.split()[1])
                        if line.startswith("adductmass"):
                            self.adductmass = ud.string_to_value(line.split()[1])
                        if line.startswith("nativezub"):
                            self.nativezub = ud.string_to_value(line.split()[1])
                        if line.startswith("nativezlb"):
                            self.nativezlb = ud.string_to_value(line.split()[1])
                        if line.startswith("poolflag"):
                            self.poolflag = ud.string_to_int(line.split()[1])
                        if line.startswith("accvol"):
                            self.detectoreffva = ud.string_to_value(line.split()[1])
                        if line.startswith("peakshapeinflate"):
                            self.inflate = ud.string_to_value(line.split()[1])
                        if line.startswith("damp"):
                            self.damp = ud.string_to_value(line.split()[1])
                        if line.startswith("linflag"):
                            self.linflag = ud.string_to_int(line.split()[1])
                        if line.startswith("cmap"):
                            self.cmap = str(line.split()[1])
                        if line.startswith("peakcmap"):
                            self.peakcmap = str(line.split()[1])
                        if line.startswith("publicationmode"):
                            self.publicationmode = ud.string_to_int(line.split()[1])
                        if line.startswith("isotopemode"):
                            self.isotopemode = ud.string_to_int(line.split()[1])
                        if line.startswith("integratelb"):
                            self.integratelb = ud.string_to_value(line.split()[1])
                        if line.startswith("integrateub"):
                            self.integrateub = ud.string_to_value(line.split()[1])
                        if line.startswith("peaknorm"):
                            self.peaknorm = ud.string_to_value(line.split()[1])
                        if line.startswith("suppression"):
                            self.suppression = ud.string_to_value(line.split()[1])

                        # IM Imports
                        if line.startswith("ccsub"):
                            self.ccsub = ud.string_to_value(line.split()[1])
                        if line.startswith("ccslb"):
                            self.ccslb = ud.string_to_value(line.split()[1])
                        if line.startswith("dtsig"):
                            self.dtsig = ud.string_to_value(line.split()[1])
                        if line.startswith("csig"):
                            self.csig = ud.string_to_value(line.split()[1])
                        if line.startswith("ccsbins"):
                            self.ccsbins = ud.string_to_value(line.split()[1])
                        if line.startswith("subbufdt"):
                            self.subbufdt = ud.string_to_value(line.split()[1])
                        if line.startswith("smoothdt"):
                            self.smoothdt = ud.string_to_value(line.split()[1])
                        if line.startswith("temp"):
                            self.temp = ud.string_to_value(line.split()[1])
                        if line.startswith("pressure"):
                            self.pressure = ud.string_to_value(line.split()[1])
                        if line.startswith("volt"):
                            self.volt = ud.string_to_value(line.split()[1])
                        if line.startswith("gasmass"):
                            self.gasmass = ud.string_to_value(line.split()[1])
                        if (line.startswith("to")) or (line.startswith("tnaught")):
                            self.to = ud.string_to_value(line.split()[1])
                        if line.startswith("driftlength"):
                            self.driftlength = ud.string_to_value(line.split()[1])
                        if line.startswith("tcal1"):
                            self.tcal1 = ud.string_to_value(line.split()[1])
                        if line.startswith("tcal2"):
                            self.tcal2 = ud.string_to_value(line.split()[1])
                        if line.startswith("edc"):
                            self.edc = ud.string_to_value(line.split()[1])
                        if line.startswith("twaveflag"):
                            self.twaveflag = ud.string_to_int(line.split()[1])
                        if line.startswith("pusher"):
                            self.pusher = ud.string_to_value(line.split()[1])
                        if line.startswith("zout"):
                            self.zout = ud.string_to_value(line.split()[1])
                        if line.startswith("nativeccsub") or line.startswith("ubnativeccs"):
                            self.nativeccsub = ud.string_to_value(line.split()[1])
                        if line.startswith("nativeccslb") or line.startswith("lbnativeccs"):
                            self.nativeccslb = ud.string_to_value(line.split()[1])
            f.close()
            self.endz = self.startz + self.numz - 1

            if os.path.isfile(self.mfile):
                self.masslist = np.loadtxt(self.mfile)
                if self.masslist.size == 1:
                    self.masslist = np.array([self.masslist])
            else:
                self.masslist = np.array([])

            if os.path.isfile(self.ofile):
                self.oligomerlist = np.genfromtxt(self.ofile, dtype='str')
                if self.oligomerlist.shape == (4,) or self.oligomerlist.shape == (5,):
                    self.oligomerlist = np.array([self.oligomerlist])
            else:
                self.oligomerlist = np.array([])

            if os.path.isfile(self.manualfile):
                self.manuallist = np.loadtxt(self.manualfile)
                if self.imflag == 0:
                    if self.manuallist.shape == (3,):
                        self.manuallist = np.array([self.manuallist])
                else:
                    if self.manuallist.shape == (5,):
                        self.manuallist = np.array([self.manuallist])
            else:
                self.manuallist = np.array([])

    def print_config(self):
        """
        Simple debugging command to read in the config file and print out its contents.
        :return: None
        """
        f = open(self.confname)
        print f.read()
        f.close()

    def default_file_names(self):
        """
        Sets the default file names. For things comming into and out of the program. In theory these can be modified,
         but it might be risky.
        :return: None
        """
        self.infname = self.outfname + "_input.dat"
        self.confname = self.outfname + "_conf.dat"
        self.mfile = self.outfname + "_mfile.dat"
        self.manualfile = self.outfname + "_manualfile.dat"
        self.ofile = self.outfname + "_ofile.dat"
        self.matchfile = self.outfname + "_match.dat"
        self.peaksfile = self.outfname + "_peaks.dat"

    def check_badness(self):
        """
        Test for a few things that will crash the program:
            Min is greater than Max for m/z, charge, mass, native charge, ccs, native ccs, dt
            Bad IM-MS calibration values.
            Peak width is zero
            m/z resolution is really small.
        :return: None
        """
        self.badtest = 0
        self.warning = ""

        # Check that mz min and max are not reversed
        if self.maxmz <= self.minmz:
            self.warning = "Max m/z is less or equal to Min m/z\nFix m/z range."
            self.badtest = 1

        # Check that bin size is not too small
        if self.mzbins < 0.01 and self.linflag != 2:
            self.warning = "Bin size is really small!\nGo ahead, but know that it will be really slow."

        # Check Charge range
        if self.endz < self.startz:
            self.warning = "Max charge is less than Min charge\nFix charge range."
            self.badtest = 1

        # Check Mass Range
        if abs(self.massub) <= abs(self.masslb):
            self.warning = "Max mass is less or equal to Min mass\nFix mass range."
            self.badtest = 1

        # Check Native Z range
        if self.nativezub <= self.nativezlb:
            self.warning = "Max native z offset is less or equal to Min native z offset\nFix native z offset range."
            self.badtest = 1

        # Check peak width
        if self.mzsig == 0:
            self.warning = "Peak width is zero\nFix peak width to be positive number"
            self.badtest = 1

        if self.imflag == 1:
            if self.nativeccsub <= self.nativeccslb:
                self.warning = "Max native ccs offset is less or equal to Min native ccs offset" \
                               "\nFix native ccs offset range."
                self.badtest = 1

            if self.maxdt < self.mindt:
                self.warning = "Max Arrival Time is less or equal to Min Arrival Time\nFix Arrival Time range."
                self.badtest = 1

            if self.ccsub < self.ccslb:
                self.warning = "CCS Upper Bound lower then CCS Lower Bound\nFix CCS range."
                self.badtest = 1

            if self.twaveflag == 0 and (self.pressure == 0 or self.driftlength == 0):
                self.warning = "Pressure and Cell Length must be nonzero for Linear Cell"
                self.badtest = 1

            if self.twaveflag == 1 and (self.tcal1 == 0 or self.tcal2 == 0 or self.edc == 0):
                self.warning = "Note: One or more T-wave calibration parameters has been set to 0" \
                               "\nCheck to make sure this is correct"

        return self.badtest, self.warning

    def default_high_res(self):
        """
        Sets some defaults for high resolution spectra. Leaves other values unchanged.
        :return: None
        """
        # Interpolate Spectrum at higher resolution
        self.mzbins = 0.05
        self.linflag = 3
        self.adductmass = 1.007276467
        # Narrow Charge and Mass Range
        self.startz = 1
        self.endz = 30
        self.massub = 30000
        self.masslb = 100
        # Increase mass resolution and lower peak shape
        self.massbins = 0.1
        self.mzsig = 0.1
        self.psfun = 0

    def default_zero_charge(self):
        """
        Sets some defaults for when the zero-charge mass spectrum itself is to be deisotoped. Leaves other values
        unchanged.
        :return: None
        """
        self.mzbins = 0
        self.linflag = 2
        self.startz = 1
        self.endz = 1
        self.massbins = 0.1
        self.mzsig = 0.4
        self.psfun = 0
        self.isotopemode = 1
        self.adductmass = 0
        self.minmz = ''
        self.maxmz = ''

    def initialize_system_paths(self):
        """
        Initialize initial paths for UniDec directories
        :return: None
        """

        print "\nInitial File Locations..."
        self.system = platform.system()
        pathtofile = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.defaultUnidecDir = os.path.join(pathtofile, 'unidec_bin')

        if self.system == 'Windows':
            self.defaultUnidecName = "UniDec.exe"
            self.defaultIMName = "UniDecIM.exe"
            print "Windows: ", self.defaultUnidecName
        elif self.system == 'Darwin':
            self.defaultUnidecName = "unidecmac"
            self.defaultIMName = "unidecmacIM"
            self.defaultUnidecDir = '/Applications/GUniDecMac.app/Contents/MacOS'
            print "Mac:", self.defaultUnidecName
        else:
            self.defaultUnidecName = "unideclinux"
            self.defaultIMName = "unideclinuxIM"
            print "Linux or other: unidec"

        self.UniDecPath = os.path.join(self.defaultUnidecDir, self.defaultUnidecName)
        if not os.path.isfile(self.UniDecPath):
            self.defaultUnidecDir = pathtofile
            self.UniDecPath = os.path.join(self.defaultUnidecDir, self.defaultUnidecName)
            if not os.path.isfile(self.UniDecPath):
                self.defaultUnidecDir = ""
                self.UniDecPath = self.defaultUnidecName
                print "Assuming " + self.defaultUnidecName + " is in system path."

        self.UniDecName = self.defaultUnidecName
        self.UniDecIMName = self.defaultIMName
        self.UniDecDir = self.defaultUnidecDir
        # self.UniDecPath = os.path.join(self.UniDecDir, self.UniDecName)
        self.UniDecIMPath = os.path.join(self.UniDecDir, self.UniDecIMName)
        self.rawreaderpath = os.path.join(self.UniDecDir, "rawreader.exe")
        self.cdcreaderpath = os.path.join(self.UniDecDir, "CDCreader.exe")
        self.defaultconfig = os.path.join(self.UniDecDir, "default_conf.dat")

        print "UniDec Path:", self.UniDecPath
