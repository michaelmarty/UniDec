import os
import numpy as np
from unidec_modules import unidectools as ud
import platform
import matplotlib.cm as cm
from matplotlib.pyplot import colormaps
import h5py
from unidec_modules.hdf5_tools import *
from unidec_modules.unidec_enginebase import version as version

__author__ = 'Michael.Marty'


def ofile_reader(path):
    oligos = []
    for line in open(path):
        a = np.array(line.split())
        string = " ".join(a[4:])
        oarray = np.array([a[0], a[1], a[2], a[3], string])
        oligos.append(oarray)
    return np.array(oligos)


# noinspection PyAttributeOutsideInit
class UniDecConfig(object):
    """
    Class containing all options and configurations for UniDec GUI and Program. Contains methods to export and import
    config to text file for running UniDec core binaries and for storing parameters for GUI.
    """

    def __init__(self):
        """
        Initialize Everything. Set default paths and run self.initialize
        :return: UniDecConfig object
        """
        self.version = version
        self.inputversion = None
        self.dtype = np.single
        self.infname = "input.dat"
        self.outfname = ""
        self.mfile = "mass.dat"
        self.manualfile = "man.dat"
        self.confname = "conf.dat"
        self.hdf_file = "default.hdf5"
        self.ofile = "ofile.dat"
        self.matchfile = "match.csv"
        self.peaksfile = "peaks.dat"
        self.dirname = ''
        self.udir = ''
        self.filename = ''
        self.extension = ''
        self.imflag = 0
        self.metamode = -2
        self.filetype = 0
        self.autotune = False

        self.time_window = 2
        self.chrom_peak_width = 2
        self.sw_time_window = 1
        self.sw_scan_offset = 10

        self.twavedict = {1: "Logarithmic", 2: "Linear", 3: "Power Law"}
        self.backgroundchoices = ["Subtract Minimum", "Subtract Line",
                                  "Subtract Curved", "Subtract Constant"]  # , "Subtract SavGol","Subtract Polynomial"]
        self.isotopechoices = ["Isotopes: Off", "Isotopes: Mono", "Isotopes: Average"]
        self.figsize = (6, 5)
        self.mass_proton = 1.007276467
        self.mass_diff_carbon = 1.0033
        self.initialize()

    def initialize(self):
        """
        Initialize configuration parameters but not paths. Runs self.default_colormaps
        :return: None
        """
        # plotting
        self.publicationmode = 1
        self.discreteplot = 0
        self.cmap = u"nipy_spectral"
        self.peakcmap = u"rainbow"
        self.spectracmap = u"rainbow"
        self.rawflag = 0

        # data prep
        self.minmz = ''
        self.maxmz = ''
        self.intscale = u"Linear"

        # Results
        self.error = 0
        self.runtime = 0

        # IM specific
        self.mindt = ''
        self.maxdt = ''
        self.zout = 0
        self.pusher = 0
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
        self.massdatnormtop = 0
        self.mfileflag = 0
        self.manualfileflag = 0
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
        self.defectparams = None
        self.defectcomparefiles = [None, None]
        self.matchtolerance = 1000

        self.avgscore = 0

        self.default_decon_params()

        self.default_colormaps()

    def default_decon_params(self):
        # Data Prep
        self.detectoreffva = 0
        self.mzbins = 0
        self.smooth = 0
        self.subbuff = 0
        self.subtype = 2
        self.intthresh = 0
        self.reductionpercent = 0

        # UniDec
        self.numit = 100
        self.zzsig = 1
        self.psig = 1
        self.beta = 0
        self.startz = 1
        self.endz = 50
        self.numz = 50
        self.mzsig = 0.85
        self.automzsig = 0
        self.psfun = 0
        self.autopsfun = 0
        self.massub = 500000
        self.masslb = 5000
        self.msig = 0
        self.molig = 0
        self.massbins = 10
        self.adductmass = 1.007276467
        self.baselineflag = 1
        self.aggressiveflag = 0
        self.noiseflag = 0
        self.isotopemode = 0
        self.orbimode = 0

        # Other
        self.mtabsig = 0
        self.poolflag = 2
        self.nativezub = 1000
        self.nativezlb = -1000
        self.inflate = 1
        self.linflag = 2
        self.integratelb = ""
        self.integrateub = ""
        self.filterwidth = 20
        self.zerolog = -12

        # Peak Selection and plotting
        self.peakwindow = 500
        self.peakthresh = 0.1
        self.peakplotthresh = 0.1
        self.separation = 0.025
        self.peaknorm = 1
        self.exwindow = 0
        self.exchoice = 0
        self.exchoicez = 1
        self.exthresh = 10
        self.exnorm = 1
        self.exnormz = 0
        self.datanorm = 1
        self.numtot = 20
        self.crossover = 100

        # IM Specific
        self.smoothdt = 0
        self.subbufdt = 0
        self.ccslb = 100
        self.ccsub = 25000
        self.nativeccsub = 20000
        self.nativeccslb = -20000
        self.dtsig = 0.2
        self.ccsbins = 100
        self.csig = 0

        self.doubledec = False
        self.kernel = ""

    def default_colormaps(self):
        """
        Get default matplotlib colormaps and set names to self.cmaps.
        A more list is set as self.cmaps2. For some reason, matplotlib v.1.5.0 allows viridis and other new colormaps
        with the colormaps() but not cm.datad. Even weirder, they work for some things and not others. Solution, two
        lists, one for the 2D plots and the other for the linear colors.
        :return: None
        """
        m = [i for i in cm.datad]
        m2 = colormaps()

        self.cmaps = sorted(m2)
        self.cmaps2 = sorted(m2)

    def config_export(self, name):
        """
        Writes config to file give in name. Typically in format: name value.

        Also exports manuallist, masslist, and oligomerlist to text files.
        :param name: File name to write to.
        :return: None
        """

        self.numz = self.endz - self.startz + 1
        f = open(name, 'w+')
        '''
        except:
            path = os.path.join(os.getcwd(), name)
            path = "\\\\?\\%s" % path
            print(path)
            f = open(path, 'w+')
            print("NOTE: Your path length might exceed the limit for Windows. Please shorten your file name.")'''
        f.write("version " + str(self.version) + "\n")
        f.write("imflag " + str(self.imflag) + "\n")
        f.write("input " + str(self.infname) + "\n")
        f.write("output " + str(self.outfname) + "\n")
        f.write("numit " + str(self.numit) + "\n")
        f.write("numz " + str(self.numz) + "\n")
        f.write("endz " + str(self.endz) + "\n")
        f.write("startz " + str(self.startz) + "\n")
        f.write("zzsig " + str(self.zzsig) + "\n")
        f.write("psig " + str(self.psig) + "\n")
        f.write("beta " + str(self.beta) + "\n")
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
                print("Need to specify mass list. Running without mass list.")
        if self.manualfileflag:
            if not ud.isempty(self.manuallist):
                f.write("manualfile " + str(self.manualfile) + "\n")
            else:
                print("Need to specify manual assignments. Running without assignments.")
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
        f.write("reductionpercent " + str(self.reductionpercent) + "\n")
        f.write("aggressive " + str(self.aggressiveflag) + "\n")
        f.write("rawflag " + str(self.rawflag) + "\n")
        f.write("adductmass " + str(self.adductmass) + "\n")
        f.write("nativezub " + str(self.nativezub) + "\n")
        f.write("nativezlb " + str(self.nativezlb) + "\n")
        f.write("poolflag " + str(self.poolflag) + "\n")
        f.write("accvol " + str(self.detectoreffva) + "\n")
        f.write("peakshapeinflate " + str(self.inflate) + "\n")
        f.write("noiseflag " + str(self.noiseflag) + "\n")
        f.write("linflag " + str(self.linflag) + "\n")
        f.write("cmap " + str(self.cmap) + "\n")
        f.write("peakcmap " + str(self.peakcmap) + "\n")
        f.write("spectracmap " + str(self.spectracmap) + "\n")
        f.write("publicationmode " + str(self.publicationmode) + "\n")
        f.write("isotopemode " + str(self.isotopemode) + "\n")
        f.write("peaknorm " + str(self.peaknorm) + "\n")
        f.write("datanorm " + str(self.datanorm) + "\n")
        f.write("baselineflag " + str(self.baselineflag) + "\n")
        f.write("orbimode " + str(self.orbimode) + "\n")
        if self.integratelb != "" and self.integrateub != "":
            try:
                f.write("integratelb " + str(self.integratelb) + "\n")
                f.write("integrateub " + str(self.integrateub) + "\n")
            except ValueError:
                print("Failed to write integation areas:", self.integratelb, self.integrateub)
                pass
        f.write("filterwidth " + str(self.filterwidth) + "\n")
        f.write("zerolog " + str(self.zerolog) + "\n")

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

            if self.twaveflag > 0:  # and (self.tcal1!=0 or self.tcal2!=0 or self.edc!=0):
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
                    print("Manual List Shape is wrong. Try using manual list tool again.")
                    print(self.manuallist.shape)
            else:
                if self.manuallist.shape[1] == 5:
                    ud.dataexport(self.manuallist, self.manualfile)
                else:
                    print("Manual List Shape is wrong. Try using manual list tool again.")
                    print(self.manuallist.shape)
        if not ud.isempty(self.oligomerlist):
            np.savetxt(self.ofile, self.oligomerlist, fmt='%s')

    def config_import(self, name):
        """
        Imports configuration from txt file. Also imports masslist, manuallist, and oligomerlist.
        :param name: File to import from.
        :return: None
        """
        if self.batchflag != 1:
            f = open(name, 'r', encoding="utf-8")
            self.manualfileflag = 0
            self.mfileflag = 0
            for line in f:
                if len(line.split()) > 1:
                    if line.startswith("version"):
                        self.inputversion = line.split()[1]
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
                        if line.startswith("psig"):
                            self.psig = ud.string_to_value(line.split()[1])
                        if line.startswith("beta"):
                            self.beta = ud.string_to_value(line.split()[1])
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
                            self.mfileflag = 1
                        if line.startswith("manualfile"):
                            # self.manualfile = line.split()[1]
                            self.manualfileflag = 1
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
                        if line.startswith("reductionpercent"):
                            self.reductionpercent = ud.string_to_value(line.split()[1])
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
                        if line.startswith("noiseflag"):
                            self.noiseflag = ud.string_to_value(line.split()[1])
                        if line.startswith("linflag"):
                            self.linflag = ud.string_to_int(line.split()[1])
                        if line.startswith("cmap"):
                            try:
                                self.cmap = str(line.split()[1], encoding="utf-8")
                            except:
                                self.cmap = str(line.split()[1])
                        if line.startswith("peakcmap"):
                            try:
                                self.peakcmap = str(line.split()[1], encoding="utf-8")
                            except:
                                self.peakcmap = str(line.split()[1])
                        if line.startswith("spectracmap"):
                            try:
                                self.spectracmap = str(line.split()[1], encoding="utf-8")
                            except:
                                self.spectracmap = str(line.split()[1])
                        if line.startswith("publicationmode"):
                            self.publicationmode = ud.string_to_int(line.split()[1])
                        if line.startswith("isotopemode"):
                            self.isotopemode = ud.string_to_int(line.split()[1])
                        if line.startswith("integratelb"):
                            self.integratelb = ud.string_to_value(line.split()[1])
                        if line.startswith("integrateub"):
                            self.integrateub = ud.string_to_value(line.split()[1])
                        if line.startswith("filterwidth"):
                            self.filterwidth = ud.string_to_value(line.split()[1])
                        if line.startswith("zerolog"):
                            self.zerolog = ud.string_to_value(line.split()[1])
                        if line.startswith("peaknorm"):
                            self.peaknorm = ud.string_to_value(line.split()[1])
                        if line.startswith("datanorm"):
                            self.datanorm = ud.string_to_value(line.split()[1])
                        if line.startswith("baselineflag"):
                            self.baselineflag = ud.string_to_value(line.split()[1])
                        if line.startswith("orbimode"):
                            self.orbimode = ud.string_to_value(line.split()[1])

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
                self.oligomerlist = ofile_reader(self.ofile)
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

    def write_hdf5(self, file_name=None):
        if file_name is None:
            file_name = self.hdf_file
        hdf = h5py.File(file_name, 'a')
        config_group = hdf.require_group("config")

        if self.metamode != -2:
            self.linflag = 2
        try:
            self.maxmz = float(self.maxmz)
        except:
            self.maxmz = 1000000
        try:
            self.minmz = float(self.minmz)
        except:
            self.minmz = 0

        cdict = {  # "input": self.infname, "output": self.outfname,
            "numit": self.numit, "version": self.version,
            "endz": self.endz, "startz": self.startz, "zzsig": self.zzsig, "psig": self.psig, "mzsig": self.mzsig,
            "beta": self.beta,
            "psfun": self.psfun, "discreteplot": self.discreteplot, "massub": self.massub, "masslb": self.masslb,
            "msig": self.msig, "molig": self.molig, "massbins": self.massbins, "mtabsig": self.mtabsig,
            "minmz": self.minmz, "maxmz": self.maxmz, "subbuff": self.subbuff, "smooth": self.smooth,
            "mzbins": self.mzbins, "peakwindow": self.peakwindow, "peakthresh": self.peakthresh,
            "peakplotthresh": self.peakplotthresh, "plotsep": self.separation, "intthresh": self.intthresh,
            "reductionpercent": self.reductionpercent,
            "aggressive": self.aggressiveflag, "rawflag": self.rawflag, "adductmass": self.adductmass,
            "nativezub": self.nativezub, "nativezlb": self.nativezlb, "poolflag": self.poolflag,
            "accvol": self.detectoreffva, "peakshapeinflate": self.inflate, "noiseflag": self.noiseflag,
            "linflag": self.linflag, "cmap": self.cmap, "peakcmap": self.peakcmap, "spectracmap": self.spectracmap,
            "publicationmode": self.publicationmode, "isotopemode": self.isotopemode, "peaknorm": self.peaknorm,
            "baselineflag": self.baselineflag, "orbimode": self.orbimode, "zout": self.zout, "pusher": self.pusher,
            "mindt": self.mindt,
            "maxdt": self.maxdt, "ccsub": self.ccsub, "ccslb": self.ccslb, "dtsig": self.dtsig, "csig": self.csig,
            "ccsbins": self.ccsbins, "subbufdt": self.subbufdt, "smoothdt": self.smoothdt,
            "ubnativeccs": self.nativeccsub, "lbnativeccs": self.nativeccslb, "twaveflag": self.twaveflag,
            "temp": self.temp, "pressure": self.pressure, "volt": self.volt,
            "tnaught": self.to, "driftlength": self.driftlength, "tcal1": self.tcal1, "tcal2": self.tcal2,
            "edc": self.edc, "gasmass": self.gasmass, "integratelb": self.integratelb,
            "integrateub": self.integrateub, "filterwidth": self.filterwidth, "zerolog": self.zerolog,
            "manualfileflag": self.manualfileflag, "mfileflag": self.mfileflag, "imflag": self.imflag,
            "exwindow": self.exwindow, "exchoice": self.exchoice, "exchoicez": self.exchoicez,
            "exthresh": self.exthresh,
            "exnorm": self.exnorm, "exnormz": self.exnormz, "metamode": self.metamode,
            "datanorm": self.datanorm, "chrom_time_window": self.time_window, "chrom_peak_width": self.chrom_peak_width,
            "sw_time_window": self.sw_time_window, "sw_scan_offset": self.sw_scan_offset
        }

        for key, value in cdict.items():
            try:
                config_group.attrs[key] = value
            except:
                print("Error with key, value:", key, value)

        if not ud.isempty(self.masslist):
            replace_dataset(config_group, "masslist", data=self.masslist)
        if not ud.isempty(self.manuallist):
            replace_dataset(config_group, "manuallist", data=self.manuallist)
        if not ud.isempty(self.oligomerlist):
            replace_dataset(config_group, "oligomerlist", data=self.oligomerlist.astype(np.string_))

        hdf.close()
        pass

    def read_attr(self, thing, string, config):
        try:
            if string in list(config.attrs.keys()):
                val = config.attrs.get(string)
                if isinstance(val, np.ndarray):
                    return val[0]
                else:
                    return val
            else:
                return thing
        except:
            return thing

    def read_hdf5(self, file_name=None):
        if file_name is None:
            file_name = self.hdf_file
        hdf = h5py.File(file_name, 'r')
        config_group = hdf.get("config")
        # self.infname = self.read_attr(self.infname, "input", config_group)
        self.inputversion = self.read_attr(self.inputversion, "version", config_group)
        self.maxmz = self.read_attr(self.maxmz, "maxmz", config_group)
        self.minmz = self.read_attr(self.minmz, "minmz", config_group)
        self.metamode = self.read_attr(self.metamode, "metamode", config_group)
        # self.outfname = self.read_attr(self.outfname, "output", config_group)
        self.numit = self.read_attr(self.numit, "numit", config_group)
        self.endz = self.read_attr(self.endz, "endz", config_group)
        self.startz = self.read_attr(self.startz, "startz", config_group)
        self.zzsig = self.read_attr(self.zzsig, "zzsig", config_group)
        self.psig = self.read_attr(self.psig, "psig", config_group)
        self.beta = self.read_attr(self.beta, "beta", config_group)
        self.mzsig = self.read_attr(self.mzsig, "mzsig", config_group)
        self.psfun = self.read_attr(self.psfun, "psfun", config_group)
        self.discreteplot = self.read_attr(self.discreteplot, "discreteplot", config_group)
        self.massub = self.read_attr(self.massub, "massub", config_group)
        self.masslb = self.read_attr(self.masslb, "masslb", config_group)
        self.molig = self.read_attr(self.molig, "molig", config_group)
        self.msig = self.read_attr(self.msig, "msig", config_group)
        self.massbins = self.read_attr(self.massbins, "massbins", config_group)
        self.mtabsig = self.read_attr(self.mtabsig, "mtabsig", config_group)
        self.subbuff = self.read_attr(self.subbuff, "subbuff", config_group)
        self.smooth = self.read_attr(self.smooth, "smooth", config_group)
        self.mzbins = self.read_attr(self.mzbins, "mzbins", config_group)
        self.peakwindow = self.read_attr(self.peakwindow, "peakwindow", config_group)
        self.peakthresh = self.read_attr(self.peakthresh, "peakthresh", config_group)
        self.peakplotthresh = self.read_attr(self.peakplotthresh, "peakplotthresh", config_group)
        self.separation = self.read_attr(self.separation, "separation", config_group)
        self.intthresh = self.read_attr(self.intthresh, "intthresh", config_group)
        self.reductionpercent = self.read_attr(self.reductionpercent, "reductionpercent", config_group)
        self.aggressiveflag = self.read_attr(self.aggressiveflag, "aggressiveflag", config_group)
        self.rawflag = self.read_attr(self.rawflag, "rawflag", config_group)
        self.adductmass = self.read_attr(self.adductmass, "adductmass", config_group)
        self.nativezub = self.read_attr(self.nativezub, "nativezub", config_group)
        self.nativezlb = self.read_attr(self.nativezlb, "nativezlb", config_group)
        self.poolflag = self.read_attr(self.poolflag, "poolflag", config_group)
        self.detectoreffva = self.read_attr(self.detectoreffva, "accvol", config_group)
        self.inflate = self.read_attr(self.inflate, "peakshapeinflate", config_group)
        self.noiseflag = self.read_attr(self.noiseflag, "noiseflag", config_group)
        self.linflag = self.read_attr(self.linflag, "linflag", config_group)
        self.cmap = self.read_attr(self.cmap, "cmap", config_group)
        self.peakcmap = self.read_attr(self.peakcmap, "peakcmap", config_group)
        self.spectracmap = self.read_attr(self.spectracmap, "spectracmap", config_group)
        self.publicationmode = self.read_attr(self.publicationmode, "publicationmode", config_group)
        self.isotopemode = self.read_attr(self.isotopemode, "isotopemode", config_group)
        self.peaknorm = self.read_attr(self.peaknorm, "peaknorm", config_group)
        self.baselineflag = self.read_attr(self.baselineflag, "baselineflag", config_group)
        self.orbimode = self.read_attr(self.orbimode, "orbimode", config_group)
        self.zout = self.read_attr(self.zout, "zout", config_group)

        self.pusher = self.read_attr(self.pusher, "pusher", config_group)
        self.mindt = self.read_attr(self.mindt, "mindt", config_group)
        self.maxdt = self.read_attr(self.maxdt, "maxdt", config_group)
        self.ccsub = self.read_attr(self.ccsub, "ccsub", config_group)
        self.ccslb = self.read_attr(self.ccslb, "ccslb", config_group)
        self.dtsig = self.read_attr(self.dtsig, "dtsig", config_group)
        self.csig = self.read_attr(self.csig, "csig", config_group)
        self.ccsbins = self.read_attr(self.ccsbins, "ccsbins", config_group)
        self.subbufdt = self.read_attr(self.subbufdt, "subbufdt", config_group)
        self.smoothdt = self.read_attr(self.smoothdt, "smoothdt", config_group)
        self.nativeccsub = self.read_attr(self.nativeccsub, "nativeccsub", config_group)
        self.nativeccslb = self.read_attr(self.nativeccslb, "nativeccslb", config_group)
        self.twaveflag = self.read_attr(self.twaveflag, "twaveflag", config_group)
        self.temp = self.read_attr(self.temp, "temp", config_group)
        self.pressure = self.read_attr(self.pressure, "pressure", config_group)
        self.volt = self.read_attr(self.volt, "volt", config_group)
        self.to = self.read_attr(self.to, "tnaught", config_group)
        self.driftlength = self.read_attr(self.driftlength, "driftlength", config_group)
        self.tcal1 = self.read_attr(self.tcal1, "tcal1", config_group)
        self.tcal2 = self.read_attr(self.tcal2, "tcal2", config_group)
        self.edc = self.read_attr(self.edc, "edc", config_group)
        self.gasmass = self.read_attr(self.gasmass, "gasmass", config_group)

        self.integratelb = self.read_attr(self.integratelb, "integratelb", config_group)
        self.integrateub = self.read_attr(self.integrateub, "integrateub", config_group)
        self.filterwidth = self.read_attr(self.filterwidth, "filterwidth", config_group)
        self.zerolog = self.read_attr(self.zerolog, "zerolog, config_group", config_group)
        self.mfileflag = self.read_attr(self.mfileflag, "mfileflag", config_group)
        self.manualfileflag = self.read_attr(self.manualfileflag, "manualfileflag", config_group)
        self.imflag = self.read_attr(self.imflag, "imflag", config_group)

        self.exchoice = self.read_attr(self.exchoice, "exchoice", config_group)
        self.exchoicez = self.read_attr(self.exchoicez, "exchoicez", config_group)
        self.exthresh = self.read_attr(self.exthresh, "exthresh", config_group)
        self.exnorm = self.read_attr(self.exnorm, "exnorm", config_group)
        self.exnormz = self.read_attr(self.exnormz, "exnormz", config_group)
        self.datanorm = self.read_attr(self.datanorm, "datanorm", config_group)
        self.exwindow = self.read_attr(self.exwindow, "exwindow", config_group)

        self.masslist = get_dataset(config_group, "masslist")
        self.manuallist = get_dataset(config_group, "manuallist")
        self.oligomerlist = get_dataset(config_group, "oligomerlist").astype(np.unicode_)

        self.time_window = self.read_attr(self.time_window, "chrom_time_window", config_group)
        self.chrom_peak_width = self.read_attr(self.chrom_peak_width, "chrom_peak_width", config_group)
        self.sw_time_window = self.read_attr(self.sw_time_window, "sw_time_window", config_group)
        self.sw_scan_offset = self.read_attr(self.sw_scan_offset, "sw_scan_offset", config_group)

        hdf.close()

    def print_config(self):
        """
        Simple debugging command to read in the config file and print out its contents.
        :return: None
        """
        f = open(self.confname)
        print(f.read())
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
        self.massdatfile = self.outfname + "_mass.txt"
        self.massgridfile = self.outfname + "_massgrid.bin"
        self.fitdatfile = self.outfname + "_fitdat.bin"
        self.errorfile = self.outfname + "_error.txt"
        self.mzgridfile = self.outfname + "_grid.bin"
        if self.filetype == 0:
            self.hdf_file = self.outfname + ".hdf5"

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

        try:
            x = float(self.maxmz)
        except:
            self.maxmz = 100000
        try:
            x = float(self.minmz)
        except:
            self.minmz = 0

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
        # if self.mzsig == 0:
        #    self.warning = "Peak width is zero\nFix peak width to be positive number"
        #    self.badtest = 1

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

            if self.twaveflag > 0 and (self.tcal1 == 0 or self.tcal2 == 0 or self.edc == 0):
                self.warning = "Note: One or more T-wave calibration parameters has been set to 0" \
                               "\nCheck to make sure this is correct"

        path = os.path.join(os.getcwd(), self.confname)
        if len(path) > 250:
            self.badtest = 1
            self.warning = "ERROR: PATH LENGTH TOO LONG: Windows sets a 260 character path length." \
                           "\nPlease shorten your file and folder path lengths."

        return self.badtest, self.warning

    def default_high_res(self):
        """
        Sets some defaults for high resolution spectra. Leaves other values unchanged.
        :return: None
        """
        # Interpolate Spectrum at higher resolution
        self.mzbins = 0
        self.linflag = 2
        self.adductmass = 1.007276467
        # Narrow Charge and Mass Range
        self.startz = 1
        self.endz = 50
        self.massub = 500000
        self.masslb = 1000
        # Increase mass resolution and lower peak shape
        self.massbins = 10
        self.mzsig = 1
        self.psfun = 0
        self.isotopemode = 0
        self.molig = 0
        self.msig = 0
        self.numit = 50
        self.zzsig = 1
        self.psig = 1
        self.poolflag = 2

    def default_low_res(self):
        """
        Sets some defaults for high resolution spectra. Leaves other values unchanged.
        :return: None
        """
        # Interpolate Spectrum at higher resolution
        self.mzbins = 1
        self.linflag = 0
        self.adductmass = 1.007276467
        # Narrow Charge and Mass Range
        self.startz = 1
        self.endz = 100
        self.massub = 1000000
        self.masslb = 10000
        # Increase mass resolution and lower peak shape
        self.massbins = 100
        self.mzsig = 10
        self.psfun = 2
        self.isotopemode = 0
        self.molig = 0
        self.msig = 0
        self.numit = 50
        self.zzsig = 1
        self.psig = 1
        self.poolflag = 2

    def default_nanodisc(self):
        """
        Sets some defaults for high resolution spectra. Leaves other values unchanged.
        :return: None
        """
        # Interpolate Spectrum at higher resolution
        self.mzbins = 0
        self.linflag = 2
        self.adductmass = 1.007276467
        # Narrow Charge and Mass Range
        self.startz = 1
        self.endz = 30
        self.massub = 200000
        self.masslb = 20000
        # Increase mass resolution and lower peak shape
        self.massbins = 10
        self.mzsig = 10
        self.psfun = 0
        self.molig = 760
        self.msig = 1
        self.zzsig = 1
        self.psig = 1
        self.isotopemode = 0
        self.poolflag = 2

    def default_isotopic_res(self):
        """
        Sets some defaults for isotopic resolution spectra. Leaves other values unchanged.
        :return: None
        """
        # Interpolate Spectrum at higher resolution
        self.mzbins = 0
        self.linflag = 2
        self.adductmass = 1.007276467
        # Narrow Charge and Mass Range
        self.startz = 1
        self.endz = 30
        self.massub = 30000
        self.masslb = 100
        # Increase mass resolution and lower peak shape
        self.massbins = 0.1
        self.mzsig = 0
        self.psfun = 0
        self.isotopemode = 1
        self.psig = 0
        self.poolflag = 2

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
        self.poolflag = 2

    def initialize_system_paths(self):
        """
        Initialize initial paths for UniDec directories
        :return: None
        """

        # print "\nInitial File Locations..."
        self.system = platform.system()
        pathtofile = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.defaultUnidecDir = os.path.join(pathtofile, 'unidec_bin')

        if self.system == 'Windows':
            self.defaultUnidecName = "UniDec.exe"
            self.h5repackfile = "h5repack.exe"
            self.opencommand = "start \"\" "
            # print "Windows: ", self.defaultUnidecName
        elif self.system == 'Darwin':
            self.defaultUnidecName = "unidecmac"
            self.h5repackfile = "h5repack"
            # self.defaultUnidecDir = '/Applications/GUniDecMac.app/Contents/MacOS'
            self.opencommand = "open "
            # print "Mac:", self.defaultUnidecName
        else:
            self.defaultUnidecName = "unideclinux"
            self.h5repackfile = "h5repack"
            self.opencommand = "gnome-open "  # TODO: Test whether this is right
            # print "Linux or other: unidec"

        def giveup():
            self.defaultUnidecDir = ""
            self.UniDecPath = self.defaultUnidecName
            print("Assuming " + self.defaultUnidecName + " is in system path.")

        self.UniDecPath = os.path.join(self.defaultUnidecDir, self.defaultUnidecName)
        if not os.path.isfile(self.UniDecPath):
            self.defaultUnidecDir = pathtofile
            self.UniDecPath = os.path.join(self.defaultUnidecDir, self.defaultUnidecName)

            if not os.path.isfile(self.UniDecPath):
                if self.system == 'Darwin':
                    self.defaultUnidecDir = '/Applications/GUniDecMac.app/Contents/MacOS'
                    self.UniDecPath = os.path.join(self.defaultUnidecDir, self.defaultUnidecName)
                    if not os.path.isfile(self.UniDecPath):
                        giveup()
                else:
                    giveup()

        self.UniDecName = self.defaultUnidecName
        self.UniDecDir = self.defaultUnidecDir
        self.UniDecPath = os.path.join(self.UniDecDir, self.UniDecName)
        # self.rawreaderpath = os.path.join(self.UniDecDir, "rawreader.exe")
        self.cdcreaderpath = os.path.join(self.UniDecDir, "CDCreader.exe")
        self.defaultconfig = os.path.join(self.UniDecDir, "default_conf.dat")
        self.masstablefile = os.path.join(self.UniDecDir, "mass_table.csv")
        self.recentfile = os.path.join(self.UniDecDir, "recent.txt")
        self.h5repackfile = os.path.join(self.UniDecDir, self.h5repackfile)
        self.presetdir = os.path.join(self.UniDecDir, "Presets")
        self.exampledatadir = os.path.join(self.UniDecDir, "Example Data")

        print("\nUniDec Path:", self.UniDecPath)

    def check_new(self, other):
        flag = False
        try:
            for item in self.__dict__:
                value = self.__dict__[item]
                value2 = other.__dict__[item]
                try:
                    try:
                        if not np.array_equal(value, value2):
                            flag = True
                            break
                    except RuntimeWarning as e:
                        print(e)
                        print(value, value2)
                except:
                    pass
        except KeyError:
            flag = True
        return flag  # not self.__dict__ == other.__dict__

    def get_colors(self, n, cmap=None):
        if cmap is None:
            cmap = self.peakcmap
        if cmap[:2] == "b'":
            cmap = cmap[2:-1]
        try:
            cmap = str(cmap, encoding="utf-8")
        except:
            pass

        colormap = cm.get_cmap(cmap, n)
        if colormap is None:
            colormap = cm.get_cmap(u"rainbow", n)
        colors = colormap(np.arange(n))
        return colors

    def get_preset_list(self):
        for dirpath, dirnames, files in os.walk(self.presetdir):
            print(files)


class DataContainer:
    def __init__(self):
        """
        Initialize DataContainer with empty arrays.
        :return: None
        """
        self.fitdat = np.array([])
        self.baseline = np.array([])
        self.fitdat2d = np.array([])
        self.rawdata = np.array([])
        self.rawdata3 = np.array([])
        self.data2 = np.array([])
        self.data3 = np.array([])
        self.massdat = np.array([])
        self.mzgrid = np.array([])
        self.massgrid = np.array([])
        self.ztab = np.array([])
        self.massccs = np.array([])
        self.ccsz = np.array([])
        self.ccsdata = np.array([])
        self.tscore = 0

    def write_hdf5(self, file_name):
        hdf = h5py.File(file_name)
        config_group = hdf.require_group("ms_data")
        immsdata = hdf.require_group("imms_data")
        replace_dataset(config_group, "raw_data_ms", data=self.rawdata)
        replace_dataset(immsdata, "raw_data_imms", data=self.rawdata3)
        replace_dataset(config_group, "fit_data", data=self.fitdat)
        replace_dataset(immsdata, "fit_data_2d", data=self.fitdat2d)
        replace_dataset(config_group, "processed_data", data=self.data2)
        replace_dataset(immsdata, "processed_data_2d", data=self.data3)
        replace_dataset(config_group, "mass_data", data=self.massdat)
        replace_dataset(config_group, "mz_grid", data=self.mzgrid)
        replace_dataset(config_group, "mass_grid", data=self.massgrid)
        replace_dataset(config_group, "charges", data=self.ztab)
        replace_dataset(immsdata, "mass_ccs_grid", data=self.massccs)
        replace_dataset(immsdata, "ccs_charge_grid", data=self.ccsz)
        replace_dataset(immsdata, "ccs_data", data=self.ccsdata)
        replace_dataset(config_group, "baseline", data=self.baseline)
        hdf.close()

    def read_hdf5(self, file_name):
        hdf = h5py.File(file_name)
        config_group = hdf.get("ms_data")
        immsdata = hdf.get("imms_data")
        self.rawdata = config_group.get("raw_data_ms")
        self.rawdata3 = immsdata.get("raw_data_imms")
        self.fitdat = config_group.get("fit_data")
        self.fitdat2d = immsdata.get("fit_data_2d")
        self.data2 = config_group.get("processed_data")
        self.data3 = immsdata.get("processed_data_2d")
        self.massdat = config_group.get("mass_data")
        self.mzgrid = config_group.get("mz_grid")
        self.massgrid = config_group.get("mass_grid")
        self.ztab = config_group.get("charges")
        self.massccs = immsdata.get("mass_ccs_grid")
        self.ccsz = immsdata.get("ccs_charge_grid")
        self.ccsdata = immsdata.get("ccs_data")
        self.baseline = config_group.get("baseline")
        hdf.close()


if __name__ == '__main__':
    fname = "test.hdf5"
    data = DataContainer()
    config = UniDecConfig()
    config.initialize_system_paths()
    config.get_preset_list()
    # config.write_hdf5(fname)
    # data.write_hdf5(fname)
    # config.read_hdf5(fname)
    # data.read_hdf5(fname)
