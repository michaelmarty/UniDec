import os
import time

import numpy as np
import unidec.tools as ud
import platform
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.pyplot import colormaps
import h5py
from unidec.modules.unidec_enginebase import version as version
from unidec.modules.hdf5_tools import replace_dataset, get_dataset

__author__ = 'Michael.Marty'


def ofile_reader(path):
    oligos = []
    for line in open(path):
        if len(line) > 1:
            a = np.array(line.split())
            string = " ".join(a[4:])
            oarray = np.array([a[0], a[1], a[2], a[3], string])
            oligos.append(oarray)
    return np.array(oligos)


def read_attr(thing, string, config_group):
    try:
        if string in list(config_group.attrs.keys()):
            val = config_group.attrs.get(string)
            if isinstance(val, np.ndarray):
                return val[0]
            else:
                return val
        else:
            return thing
    except Exception:
        return thing


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

        # File names and paths
        self.system = platform.system()
        self.infname = "input.dat"
        self.outfname = ""
        self.mfile = "mass.dat"
        self.manualfile = "man.dat"
        self.smashfile = "smash.dat"
        self.confname = "conf.dat"
        self.hdf_file = "default.hdf5"
        self.ofile = "ofile.dat"
        self.matchfile = "match.csv"
        self.peaksfile = "peaks.dat"
        self.dirname = ''
        self.udir = ''
        self.filename = ''
        self.extension = ''
        self.deconfile = ''
        self.errorfile = ''
        self.fitdatfile = ''
        self.massgridfile = ''
        self.massdatfile = ''
        self.cdrawextracts = ''
        self.cdchrom = ''
        self.mzgridfile = ''
        self.cdcreaderpath = ''
        self.UniDecPath = ''
        self.UniDecDir = ''
        self.UniDecName = ''
        self.defaultUnidecDir = ''
        self.mplstyledir = ''
        self.mplstylefile = 'UniDec.mplstyle'
        self.opencommand = ''
        self.defaultUnidecName = ''
        self.iconfile = ''
        self.toplogofile = ''
        self.exampledatadirUC = ''
        self.exampledatadirCD = ''
        self.exampledatadir = ''
        self.presetdirCD = ''
        self.presetdir = ''
        self.h5repackfile = ''
        self.recentfileCD = ''
        self.recentfile = ''
        self.masstablefile = ''
        self.defaultconfig = ''

        self.imflag = 0
        self.cdmsflag = 0
        self.metamode = -2
        self.filetype = 0
        self.autotune = False

        self.time_window = 2
        self.chrom_peak_width = 2
        self.sw_time_window = 1
        self.sw_scan_offset = 10
        self.time_start = ""
        self.time_end = ""

        self.twavedict = {1: "Logarithmic", 2: "Linear", 3: "Power Law", 4: "SLIM poly2", 5: "SLIM poly3"}
        self.backgroundchoices = ["Subtract Minimum", "Subtract Line",
                                  "Subtract Curved", "Subtract Constant"]  # , "Subtract SavGol","Subtract Polynomial"]
        self.isotopechoices = ["Isotopes: Off", "Isotopes: Mono", "Isotopes: Average"]
        self.figsize = (6, 5)
        self.mass_proton = 1.007276467
        self.mass_diff_carbon = 1.0033
        self.polarity = "Positive"

        # Initialize Parameters
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
        self.tcal3 = 0
        self.tcal4 = 0
        self.edc = 1.57
        self.gasmass = 4.002602
        self.twaveflag = 0
        self.ccsconstant = 0

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
        self.smashlist = []
        self.smashflag = 0
        self.zoffs = []
        self.smashrange = []
        self.massoffset = 0
        self.extractshape = 0
        self.gridparams = None
        self.griddecon = None
        self.defectparams = None
        self.defectcomparefiles = [None, None]
        self.matchtolerance = 1000

        self.avgscore = 0

        # Data Prep
        self.detectoreffva = 0
        self.mzbins = 0
        self.smooth = 0
        self.subbuff = 0
        self.subtype = 2
        self.intthresh = 0
        self.reductionpercent = 0

        # unidec
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
        self.psfunz = 0
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
        self.normthresh = 1
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
        self.compressflag = 1

        # Reused by IM and CD
        self.csig = 0
        # Charge Detection
        self.CDslope = 0.2074
        self.CDzbins = 1
        self.CDres = 0
        self.CDScanStart = -1
        self.CDScanEnd = -1
        self.CDiitflag = False

        # Hadamard Transform Parameters
        self.htmode = False
        self.htseq = ""
        self.htbit = "5"
        self.HTksmooth = 0
        self.HTanalysistime = 38.0
        self.HTcycletime = 2.0
        self.HTtimepad = 0.0
        self.HTtimeshift = 7.0
        self.HTcycleindex = -1
        self.HTxaxis = "Time"
        self.CDScanCompress = 3
        self.HTmaxscans = -1
        self.HToutputub = -1
        self.HToutputlb = -1

        self.demultiplexmode = "HT"
        self.demultiplexchoices = ["HT", "mHT", "FT", "aFT"]
        self.FTstart = 5
        self.FTend = 1000
        self.FTflatten = True
        self.FTapodize = 1
        self.FTsmooth = 0
        self.HTmaskn = 1000
        self.HTwin = 5
        self.showlegends = True

        self.doubledec = False
        self.kernel = ""

        self.cmaps = None
        self.cmaps2 = None

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
        self.tcal3 = 0
        self.tcal4 = 0
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
        self.smashlist = []
        self.smashflag = 0
        self.zoffs = []
        self.smashrange = []
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

        # unidec
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
        self.psfunz = 0
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
        self.normthresh = 1
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
        self.compressflag = 1

        # Reused by IM and CD
        self.csig = 0
        # Charge Detection
        self.CDslope = 0.2074
        self.CDzbins = 1
        self.CDres = 0
        self.CDScanStart = -1
        self.CDScanEnd = -1
        self.HTksmooth = 0
        self.CDScanCompress = 0
        self.HTcycleindex = -1
        self.HTxaxis = "Time"

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
        f.write("cdmsflag " + str(self.cdmsflag) + "\n")
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
        f.write("zpsfn " + str(self.psfunz) + "\n")
        f.write("discreteplot " + str(self.discreteplot) + "\n")
        f.write("smashflag " + str(self.smashflag) + "\n")
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
        f.write("normthresh" + str(int(self.normthresh)) + "\n")
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

        f.write("doubledec " + str(int(self.doubledec)) + "\n")
        f.write("kernel " + str(self.kernel) + "\n")

        f.write("CDslope " + str(self.CDslope) + "\n")
        f.write("CDzbins " + str(self.CDzbins) + "\n")
        f.write("CDres " + str(self.CDres) + "\n")
        f.write("CDScanStart " + str(self.CDScanStart) + "\n")
        f.write("CDScanEnd " + str(self.CDScanEnd) + "\n")
        f.write("HTksmooth " + str(self.HTksmooth) + "\n")
        f.write("CDScanCompress " + str(self.CDScanCompress) + "\n")
        f.write("HTmaxscans " + str(self.HTmaxscans) + "\n")
        f.write("HToutlb " + str(self.HToutputlb) + "\n")
        f.write("HToutub " + str(self.HToutputub) + "\n")
        f.write("HTtimepad " + str(self.HTtimepad) + "\n")
        f.write("HTanalysistime " + str(self.HTanalysistime) + "\n")
        f.write("HTxaxis " + str(self.HTxaxis) + "\n")
        f.write("HTcycleindex " + str(self.HTcycleindex) + "\n")
        f.write("HTcylctime " + str(self.HTcycletime) + "\n")
        f.write("HTtimeshift " + str(self.HTtimeshift) + "\n")
        f.write("htbit " + str(self.htbit) + "\n")
        f.write("FTstart " + str(self.FTstart) + "\n")
        f.write("FTend " + str(self.FTend) + "\n")
        f.write("FTflatten " + str(int(self.FTflatten)) + "\n")
        f.write("FTapodize " + str(int(self.FTapodize)) + "\n")
        f.write("demultiplexmode " + str(self.demultiplexmode) + "\n")
        f.write("FTsmooth " + str(self.FTsmooth) + "\n")
        f.write("HTmaskn " + str(self.HTmaskn) + "\n")
        f.write("HTwin " + str(self.HTwin) + "\n")

        f.write("csig " + str(self.csig) + "\n")
        f.write("smoothdt " + str(self.smoothdt) + "\n")
        f.write("subbufdt " + str(self.subbufdt) + "\n")
        if self.mindt != '' or self.maxdt != '':
            f.write("zout " + str(self.zout) + "\n")
            f.write("pusher " + str(self.pusher) + "\n")
            f.write("mindt " + str(self.mindt) + "\n")
            f.write("maxdt " + str(self.maxdt) + "\n")
            f.write("ccsub " + str(self.ccsub) + "\n")
            f.write("ccslb " + str(self.ccslb) + "\n")
            f.write("dtsig " + str(self.dtsig) + "\n")
            f.write("ccsbins " + str(self.ccsbins) + "\n")
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
                f.write("tcal3 " + str(self.tcal3) + "\n")
                f.write("tcal4 " + str(self.tcal4) + "\n")
                f.write("edc " + str(self.edc) + "\n")
                f.write("gasmass " + str(self.gasmass) + "\n")

        f.close()
        # Export Mass List
        if not ud.isempty(self.masslist):
            ud.dataexport(self.masslist, self.mfile)
        # Export Manual List
        if not ud.isempty(self.manuallist):
            self.manuallist = np.array(self.manuallist)
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
        # Export Smash list
        if not ud.isempty(self.smashlist):
            print("Exporting Smashlist:", self.smashlist, self.smashfile)
            self.smashlist = np.array(self.smashlist)
            if self.imflag == 0 and self.cdmsflag == 0:
                if self.smashlist.shape[1] == 2:
                    ud.dataexport(self.smashlist, self.smashfile)
                else:
                    print("Manual List Shape is wrong. Try using manual list tool again.")
                    print(self.smashlist.shape)
            else:
                if self.smashlist.shape[1] == 4:
                    ud.dataexport(self.smashlist, self.smashfile)
                else:
                    print("Manual List Shape is wrong. Try using manual list tool again.")
                    print(self.smashlist.shape)

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

                        if line.startswith("CDslope"):
                            self.CDslope = ud.string_to_value(line.split()[1])
                        if line.startswith("CDzbins"):
                            self.CDzbins = ud.string_to_value(line.split()[1])
                        if line.startswith("CDres"):
                            self.CDres = ud.string_to_value(line.split()[1])
                        if line.startswith("CDScanStart"):
                            self.CDScanStart = ud.string_to_int(line.split()[1])
                        if line.startswith("CDScanEnd"):
                            self.CDScanEnd = ud.string_to_int(line.split()[1])
                        if line.startswith("HTksmooth"):
                            self.HTksmooth = ud.string_to_value(line.split()[1])
                        if line.startswith("HTtimeshift"):
                            self.HTtimeshift = ud.string_to_value(line.split()[1])
                        if line.startswith("HTtimepad"):
                            self.HTtimepad = ud.string_to_value(line.split()[1])
                        if line.startswith("HTanalysistime"):
                            self.HTanalysistime = ud.string_to_value(line.split()[1])
                        if line.startswith("HTxaxis"):
                            self.HTxaxis = line.split()[1]
                        if line.startswith("HTcycleindex"):
                            self.HTcycleindex = ud.string_to_value(line.split()[1])
                        if line.startswith("HTcylctime"):
                            self.HTcycletime = ud.string_to_value(line.split()[1])
                        if line.startswith("htbit"):
                            self.htbit = str(line.split()[1])
                        if line.startswith("CDScanCompress"):
                            self.CDScanCompress = ud.string_to_value(line.split()[1])
                        if line.startswith("HTmaxscans"):
                            self.HTmaxscans = ud.string_to_value(line.split()[1])
                        if line.startswith("HToutlb"):
                            self.HToutputlb = ud.string_to_value(line.split()[1])
                        if line.startswith("HToutub"):
                            self.HToutputub = ud.string_to_value(line.split()[1])
                        if line.startswith("demultiplexmode"):
                            self.demultiplexmode = line.split()[1]
                        if line.startswith("FTstart"):
                            self.FTstart = ud.string_to_value(line.split()[1])
                        if line.startswith("FTend"):
                            self.FTend = ud.string_to_value(line.split()[1])
                        if line.startswith("FTflatten"):
                            self.FTflatten = ud.string_to_int(line.split()[1])
                        if line.startswith("FTapodize"):
                            self.FTapodize = ud.string_to_int(line.split()[1])
                        if line.startswith("FTsmooth"):
                            self.FTsmooth = ud.string_to_value(line.split()[1])
                        if line.startswith("HTmaskn"):
                            self.HTmaskn = ud.string_to_value(line.split()[1])
                        if line.startswith("HTwin"):
                            self.HTwin = ud.string_to_value(line.split()[1])
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
                        if line.startswith("zpsfn"):
                            self.psfunz = ud.string_to_int(line.split()[1])
                        if line.startswith("discreteplot"):
                            self.discreteplot = ud.string_to_int(line.split()[1])
                        if line.startswith("smashflag"):
                            self.smashflag = ud.string_to_int(line.split()[1])
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
                        if line.startswith("normthresh"):
                            self.normthresh = ud.string_to_int(line.split()[1])
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
                            self.cmap = str(line.split()[1])
                        if line.startswith("peakcmap"):
                            self.peakcmap = str(line.split()[1])
                        if line.startswith("spectracmap"):
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
                            self.peaknorm = ud.string_to_int(line.split()[1])
                        if line.startswith("datanorm"):
                            self.datanorm = ud.string_to_int(line.split()[1])
                        if line.startswith("baselineflag"):
                            self.baselineflag = ud.string_to_value(line.split()[1])
                        if line.startswith("orbimode"):
                            self.orbimode = ud.string_to_int(line.split()[1])
                        if line.startswith("doubledec"):
                            ddval = ud.string_to_int(line.split()[1])
                            if ddval == 0:
                                self.doubledec = False
                            else:
                                self.doubledec = True
                        if line.startswith("kernel"):
                            self.kernel = line.strip()[7:]

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
                        if line.startswith("tcal3"):
                            self.tcal3 = ud.string_to_value(line.split()[1])
                        if line.startswith("tcal4"):
                            self.tcal4 = ud.string_to_value(line.split()[1])
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

            # Import massfile
            if os.path.isfile(self.mfile):
                self.masslist = np.loadtxt(self.mfile)
                if self.masslist.size == 1:
                    self.masslist = np.array([self.masslist])
            else:
                self.masslist = np.array([])

            # Import Ofile
            if os.path.isfile(self.ofile):
                self.oligomerlist = ofile_reader(self.ofile)
                if self.oligomerlist.shape == (4,) or self.oligomerlist.shape == (5,):
                    self.oligomerlist = np.array([self.oligomerlist])
            else:
                self.oligomerlist = np.array([])

            # Import manual file
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

            # Import smash file
            self.import_smashfile()

    def import_smashfile(self, filename=None):
        if filename is None:
            filename = self.smashfile
        else:
            self.smashfile = filename
        # Import smash file
        if os.path.isfile(self.smashfile):
            self.smashlist = np.loadtxt(self.smashfile)
            if self.imflag == 0 and self.cdmsflag == 0:
                if self.smashlist.shape == (2,):
                    self.smashlist = np.array([self.smashlist])
            else:
                if self.smashlist.shape == (4,):
                    self.smashlist = np.array([self.smashlist])
        else:
            self.smashlist = np.array([])

    def get_config_dict(self):
        if self.metamode != -2:
            self.linflag = 2
        try:
            self.maxmz = float(self.maxmz)
        except Exception:
            self.maxmz = 1000000
        try:
            self.minmz = float(self.minmz)
        except Exception:
            self.minmz = 0

        if self.tcal3 == "":
            self.tcal3 = 0
        if self.tcal4 == "":
            self.tcal4 = 0

        cdict = {  # "input": self.infname, "output": self.outfname,
            "numit": self.numit, "version": self.version,
            "endz": self.endz, "startz": self.startz, "zzsig": self.zzsig, "psig": self.psig, "mzsig": self.mzsig,
            "beta": self.beta,
            "psfun": self.psfun, "psfunz": self.psfunz, "discreteplot": self.discreteplot, "massub": self.massub,
            "masslb": self.masslb,
            "msig": self.msig, "molig": self.molig, "massbins": self.massbins, "mtabsig": self.mtabsig,
            "minmz": self.minmz, "maxmz": self.maxmz, "subbuff": self.subbuff, "smooth": self.smooth,
            "mzbins": self.mzbins, "peakwindow": self.peakwindow, "peakthresh": self.peakthresh,
            "normthresh": self.normthresh,
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
            "tcal3": self.tcal3, "tcal4": self.tcal4,
            "edc": self.edc, "gasmass": self.gasmass, "integratelb": self.integratelb,
            "integrateub": self.integrateub, "filterwidth": self.filterwidth, "zerolog": self.zerolog,
            "manualfileflag": self.manualfileflag, "mfileflag": self.mfileflag, "imflag": self.imflag,
            "cdmsflag": self.cdmsflag,
            "exwindow": self.exwindow, "exchoice": self.exchoice, "exchoicez": self.exchoicez,
            "exthresh": self.exthresh,
            "exnorm": self.exnorm, "exnormz": self.exnormz, "metamode": self.metamode,
            "datanorm": self.datanorm, "chrom_time_window": self.time_window, "chrom_peak_width": self.chrom_peak_width,
            "sw_time_window": self.sw_time_window, "sw_scan_offset": self.sw_scan_offset, "time_start": self.time_start,
            "time_end": self.time_end
        }
        return cdict

    def write_hdf5(self, file_name=None):
        if file_name is None:
            file_name = self.hdf_file
        hdf = h5py.File(file_name, 'a')
        config_group = hdf.require_group("config")

        cdict = self.get_config_dict()

        for key, value in cdict.items():
            if value is not None:
                try:
                    # if value is None:
                    #    value = ""
                    config_group.attrs[key] = value
                except Exception as e:
                    print("Error with key, value:", key, value, e)

        if not ud.isempty(self.masslist):
            replace_dataset(config_group, "masslist", data=self.masslist)
        if not ud.isempty(self.manuallist):
            replace_dataset(config_group, "manuallist", data=self.manuallist)
        if not ud.isempty(self.oligomerlist):
            self.oligomerlist = np.array(self.oligomerlist)
            replace_dataset(config_group, "oligomerlist", data=self.oligomerlist.astype(np.string_))
            np.savetxt(self.ofile, self.oligomerlist, fmt='%s')

        hdf.close()
        pass

    def read_hdf5(self, file_name=None):
        if file_name is None:
            file_name = self.hdf_file
        hdf = h5py.File(file_name, 'r')
        config_group = hdf.get("config")
        # self.infname = self.read_attr(self.infname, "input", config_group)
        self.inputversion = read_attr(self.inputversion, "version", config_group)
        self.maxmz = read_attr(self.maxmz, "maxmz", config_group)
        self.minmz = read_attr(self.minmz, "minmz", config_group)
        self.metamode = read_attr(self.metamode, "metamode", config_group)
        # self.outfname = self.read_attr(self.outfname, "output", config_group)
        self.numit = read_attr(self.numit, "numit", config_group)
        self.endz = read_attr(self.endz, "endz", config_group)
        self.startz = read_attr(self.startz, "startz", config_group)
        self.zzsig = read_attr(self.zzsig, "zzsig", config_group)
        self.psig = read_attr(self.psig, "psig", config_group)
        self.beta = read_attr(self.beta, "beta", config_group)
        self.mzsig = read_attr(self.mzsig, "mzsig", config_group)
        self.psfun = read_attr(self.psfun, "psfun", config_group)
        self.psfunz = read_attr(self.psfunz, "psfunz", config_group)
        self.discreteplot = read_attr(self.discreteplot, "discreteplot", config_group)
        self.massub = read_attr(self.massub, "massub", config_group)
        self.masslb = read_attr(self.masslb, "masslb", config_group)
        self.molig = read_attr(self.molig, "molig", config_group)
        self.msig = read_attr(self.msig, "msig", config_group)
        self.massbins = read_attr(self.massbins, "massbins", config_group)
        self.mtabsig = read_attr(self.mtabsig, "mtabsig", config_group)
        self.subbuff = read_attr(self.subbuff, "subbuff", config_group)
        self.smooth = read_attr(self.smooth, "smooth", config_group)
        self.mzbins = read_attr(self.mzbins, "mzbins", config_group)
        self.peakwindow = read_attr(self.peakwindow, "peakwindow", config_group)
        self.peakthresh = read_attr(self.peakthresh, "peakthresh", config_group)
        self.normthresh = read_attr(self.normthresh, "normthresh", config_group)
        self.peakplotthresh = read_attr(self.peakplotthresh, "peakplotthresh", config_group)
        self.separation = read_attr(self.separation, "separation", config_group)
        self.intthresh = read_attr(self.intthresh, "intthresh", config_group)
        self.reductionpercent = read_attr(self.reductionpercent, "reductionpercent", config_group)
        self.aggressiveflag = read_attr(self.aggressiveflag, "aggressiveflag", config_group)
        self.rawflag = read_attr(self.rawflag, "rawflag", config_group)
        self.adductmass = read_attr(self.adductmass, "adductmass", config_group)
        self.nativezub = read_attr(self.nativezub, "nativezub", config_group)
        self.nativezlb = read_attr(self.nativezlb, "nativezlb", config_group)
        self.poolflag = read_attr(self.poolflag, "poolflag", config_group)
        self.detectoreffva = read_attr(self.detectoreffva, "accvol", config_group)
        self.inflate = read_attr(self.inflate, "peakshapeinflate", config_group)
        self.noiseflag = read_attr(self.noiseflag, "noiseflag", config_group)
        self.linflag = read_attr(self.linflag, "linflag", config_group)
        self.cmap = read_attr(self.cmap, "cmap", config_group)
        self.peakcmap = read_attr(self.peakcmap, "peakcmap", config_group)
        self.spectracmap = read_attr(self.spectracmap, "spectracmap", config_group)
        self.publicationmode = read_attr(self.publicationmode, "publicationmode", config_group)
        self.isotopemode = read_attr(self.isotopemode, "isotopemode", config_group)
        self.peaknorm = read_attr(self.peaknorm, "peaknorm", config_group)
        self.baselineflag = read_attr(self.baselineflag, "baselineflag", config_group)
        self.orbimode = read_attr(self.orbimode, "orbimode", config_group)
        self.zout = read_attr(self.zout, "zout", config_group)

        self.pusher = read_attr(self.pusher, "pusher", config_group)
        self.mindt = read_attr(self.mindt, "mindt", config_group)
        self.maxdt = read_attr(self.maxdt, "maxdt", config_group)
        self.ccsub = read_attr(self.ccsub, "ccsub", config_group)
        self.ccslb = read_attr(self.ccslb, "ccslb", config_group)
        self.dtsig = read_attr(self.dtsig, "dtsig", config_group)
        self.csig = read_attr(self.csig, "csig", config_group)
        self.ccsbins = read_attr(self.ccsbins, "ccsbins", config_group)
        self.subbufdt = read_attr(self.subbufdt, "subbufdt", config_group)
        self.smoothdt = read_attr(self.smoothdt, "smoothdt", config_group)
        self.nativeccsub = read_attr(self.nativeccsub, "nativeccsub", config_group)
        self.nativeccslb = read_attr(self.nativeccslb, "nativeccslb", config_group)
        self.twaveflag = read_attr(self.twaveflag, "twaveflag", config_group)
        self.temp = read_attr(self.temp, "temp", config_group)
        self.pressure = read_attr(self.pressure, "pressure", config_group)
        self.volt = read_attr(self.volt, "volt", config_group)
        self.to = read_attr(self.to, "tnaught", config_group)
        self.driftlength = read_attr(self.driftlength, "driftlength", config_group)
        self.tcal1 = read_attr(self.tcal1, "tcal1", config_group)
        self.tcal2 = read_attr(self.tcal2, "tcal2", config_group)
        self.tcal3 = read_attr(self.tcal3, "tcal3", config_group)
        self.tcal4 = read_attr(self.tcal4, "tcal4", config_group)
        self.edc = read_attr(self.edc, "edc", config_group)
        self.gasmass = read_attr(self.gasmass, "gasmass", config_group)

        self.integratelb = read_attr(self.integratelb, "integratelb", config_group)
        self.integrateub = read_attr(self.integrateub, "integrateub", config_group)
        self.filterwidth = read_attr(self.filterwidth, "filterwidth", config_group)
        self.zerolog = read_attr(self.zerolog, "zerolog, config_group", config_group)
        self.mfileflag = read_attr(self.mfileflag, "mfileflag", config_group)
        self.manualfileflag = read_attr(self.manualfileflag, "manualfileflag", config_group)
        self.imflag = read_attr(self.imflag, "imflag", config_group)
        self.cdmsflag = read_attr(self.cdmsflag, "cdmsflag", config_group)

        self.exchoice = read_attr(self.exchoice, "exchoice", config_group)
        self.exchoicez = read_attr(self.exchoicez, "exchoicez", config_group)
        self.exthresh = read_attr(self.exthresh, "exthresh", config_group)
        self.exnorm = read_attr(self.exnorm, "exnorm", config_group)
        self.exnormz = read_attr(self.exnormz, "exnormz", config_group)
        self.datanorm = read_attr(self.datanorm, "datanorm", config_group)
        self.exwindow = read_attr(self.exwindow, "exwindow", config_group)

        self.masslist = get_dataset(config_group, "masslist")
        self.manuallist = get_dataset(config_group, "manuallist")
        self.oligomerlist = get_dataset(config_group, "oligomerlist").astype(np.unicode_)

        self.time_window = read_attr(self.time_window, "chrom_time_window", config_group)
        self.chrom_peak_width = read_attr(self.chrom_peak_width, "chrom_peak_width", config_group)
        self.sw_time_window = read_attr(self.sw_time_window, "sw_time_window", config_group)
        self.sw_scan_offset = read_attr(self.sw_scan_offset, "sw_scan_offset", config_group)

        self.time_start = read_attr(self.time_start, "time_start", config_group)
        self.time_end = read_attr(self.time_end, "time_end", config_group)

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
        self.smashfile = self.outfname + "_smashfile.dat"
        self.ofile = self.outfname + "_ofile.dat"
        self.matchfile = self.outfname + "_match.dat"
        self.peaksfile = self.outfname + "_peaks.dat"
        self.massdatfile = self.outfname + "_mass.txt"
        self.massgridfile = self.outfname + "_massgrid.bin"
        self.fitdatfile = self.outfname + "_fitdat.bin"
        self.errorfile = self.outfname + "_error.txt"
        self.deconfile = self.outfname + "_decon.txt"
        self.mzgridfile = self.outfname + "_grid.bin"
        self.cdrawextracts = self.outfname + "_rawdata.npz"
        self.cdchrom = self.outfname + "_chroms.txt"
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
        badest = 0
        warning = ""

        try:
            x = float(self.maxmz)
        except Exception:
            self.maxmz = 100000
        try:
            x = float(self.minmz)
        except Exception:
            self.minmz = 0

        # Check that mz min and max are not reversed
        if self.maxmz <= self.minmz:
            warning = "Max m/z is less or equal to Min m/z\nFix m/z range."
            badest = 1

        # Check that bin size is not too small
        if self.mzbins < 0.01 and self.linflag != 2:
            warning = "Bin size is really small!\nGo ahead, but know that it will be really slow."

        # Check Charge range
        if self.endz < self.startz:
            warning = "Max charge is less than Min charge\nFix charge range."
            badest = 1

        # Check Mass Range
        if abs(self.massub) <= abs(self.masslb):
            warning = "Max mass is less or equal to Min mass\nFix mass range."
            badest = 1

        # Check Native Z range
        if self.nativezub <= self.nativezlb:
            warning = "Max native z offset is less or equal to Min native z offset\nFix native z offset range."
            badest = 1

        # Check peak width
        # if self.mzsig == 0:
        #    warning = "Peak width is zero\nFix peak width to be positive number"
        #    badest = 1

        if self.imflag == 1:
            if self.nativeccsub <= self.nativeccslb:
                warning = "Max native ccs offset is less or equal to Min native ccs offset" \
                          "\nFix native ccs offset range."
                badest = 1

            if self.maxdt < self.mindt:
                warning = "Max Arrival Time is less or equal to Min Arrival Time\nFix Arrival Time range."
                badest = 1

            if self.ccsub < self.ccslb:
                warning = "CCS Upper Bound lower then CCS Lower Bound\nFix CCS range."
                badest = 1

            if self.twaveflag > 0 and (self.tcal1 == 0 or self.tcal2 == 0 or self.edc == 0):
                warning = "Note: One or more T-wave calibration parameters has been set to 0" \
                          "\nCheck to make sure this is correct"

        path = os.path.join(os.getcwd(), self.confname)
        if len(path) > 250:
            badest = 1
            warning = "ERROR: PATH LENGTH TOO LONG: Windows sets a 260 character path length." \
                      "\nPlease shorten your file and folder path lengths."

        return badest, warning

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
        pathtofile = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.defaultUnidecDir = os.path.join(pathtofile, 'bin')

        if self.system == 'Windows':
            self.defaultUnidecName = "unidec.exe"
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
            print(pathtofile)

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
                    self.defaultUnidecDir = os.path.split(pathtofile)[0]
                    self.UniDecPath = os.path.join(self.defaultUnidecDir, self.defaultUnidecName)
                    if not os.path.isfile(self.UniDecPath):
                        giveup()

        self.UniDecName = self.defaultUnidecName
        self.UniDecDir = self.defaultUnidecDir
        self.UniDecPath = os.path.join(self.UniDecDir, self.UniDecName)
        # self.rawreaderpath = os.path.join(self.UniDecDir, "rawreader.exe")
        self.cdcreaderpath = os.path.join(self.UniDecDir, "CDCreader.exe")
        self.defaultconfig = os.path.join(self.UniDecDir, "default_conf.dat")
        self.masstablefile = os.path.join(self.UniDecDir, "mass_table.csv")
        self.recentfile = os.path.join(self.UniDecDir, "recent.txt")
        self.recentfileCD = os.path.join(self.UniDecDir, "recentCD.txt")
        self.h5repackfile = os.path.join(self.UniDecDir, self.h5repackfile)
        self.presetdir = os.path.join(self.UniDecDir, "Presets")
        self.presetdirCD = os.path.join(self.UniDecDir, "Presets", "CDMS")
        self.exampledatadir = os.path.join(self.UniDecDir, "Example Data")
        self.exampledatadirCD = os.path.join(self.UniDecDir, "Example Data", "CDMS")
        self.exampledatadirUC = os.path.join(self.UniDecDir, "Example Data", "UniChrom")
        self.toplogofile = os.path.join(self.UniDecDir, "UniDecLogoMR.png")
        self.iconfile = os.path.join(self.UniDecDir, "logo.ico")
        self.mplstyledir = os.path.join(self.UniDecDir, "PlotStyles")
        self.mplstylefile = os.path.join(self.mplstyledir, "UniDec.mplstyle")

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
                except Exception:
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
        except Exception:
            pass

        # colormap = cm.get_cmap(cmap, n)
        colormap = mpl.colormaps[cmap].resampled(n)
        if colormap is None:
            # colormap = cm.get_cmap(u"rainbow", n)
            colormap = mpl.colormaps[u"rainbow"].resampled(n)
        colors = colormap(np.arange(n))
        return colors

    def get_preset_list(self):
        for dirpath, dirnames, files in os.walk(self.presetdir):
            print(files)


class OligomerContainer:
    def __init__(self):
        self.oligonames = np.array([])
        self.oligomernames = np.array([])
        self.oligomasslist = np.array([])
        self.oligomerlist = np.array([])

    def make_oligomers(self, isolated=False, oligomerlist=None, minsites=None, maxsites=None):
        print("Starting to make oligomers. Isolated=", isolated)
        stime = time.perf_counter()
        self.oligomerlist = oligomerlist
        if not isolated:
            self.oligomasslist, self.oligonames = ud.make_all_matches(oligomerlist)
        else:
            self.oligomasslist, self.oligonames = ud.make_isolated_match(oligomerlist)
        if minsites is not None:
            sums = np.sum(self.oligonames, axis=1)
            b1 = sums >= minsites
            self.oligomasslist = self.oligomasslist[b1]
            self.oligonames = self.oligonames[b1]
        if maxsites is not None:
            sums = np.sum(self.oligonames, axis=1)
            b1 = sums <= maxsites
            self.oligomasslist = self.oligomasslist[b1]
            self.oligonames = self.oligonames[b1]
        print("Oligomers Made in ", time.perf_counter() - stime, "s")

    def pair_glyco(self):
        self.oligomasslist, self.oligonames = ud.pair_glyco_matches(self.oligomasslist, self.oligonames,
                                                                    self.oligomerlist)

    def get_alts(self, pks, tolerance=10):
        altmasses = []
        altindexes = []
        matchcounts = []
        for p in pks.peaks:
            m = p.mass
            b1 = np.abs(self.oligomasslist - m) < tolerance
            altmasses.append(self.oligomasslist[b1])
            altindexes.append(self.oligonames[b1])
            matchcounts.append(np.sum(b1))
        return altmasses, altindexes, np.array(matchcounts)


class DataContainer:
    def __init__(self):
        """
        Initialize DataContainer with empty arrays.
        :return: None
        """
        self.fitdat = np.array([])  # Fit of m/z data from deconvolution
        self.baseline = np.array([])  # Baseline of the m/z data
        self.fitdat2d = np.array([])  # Fit of 2D data from deconvolution
        self.rawdata = np.array([])  # Raw m/z data
        self.rawdata3 = np.array([])  # Raw IMMS data
        self.data2 = np.array([])  # Processed m/z data
        self.data3 = np.array([])  # Processed IMMS data
        self.massdat = np.array([])  # Deconvolved 1D mass data
        self.mzgrid = np.array([])  # Deconvolved 2D m/z vs charge grid
        self.massgrid = np.array([])  # Deconvolved 2D mass vs charge grid
        self.mzmassgrid = np.array([])  # Deconvolved 2D m/z vs mass grid
        self.ztab = np.array([])  # Charge state table
        self.massccs = np.array([])  # Mass vs CCS table
        self.ccsz = np.array([])  # CCS vs charge table
        self.ccsdata = np.array([])  # CCS data in 1D
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


class Chromatogram:
    def __init__(self):
        self.chromdat = np.array([])
        self.decondat = np.array([])
        self.ccsdat = np.array([])
        self.label = ""
        self.color = "#000000"
        self.index = 0
        self.mzrange = [-1, -1]
        self.zrange = [-1, -1]
        self.sarray = [-1, -1, -1, -1]
        self.sum = -1
        self.ignore = 0
        self.ht = False
        self.massrange = [-1, -1]
        self.massmode = False
        self.dataobj = None
        self.swoopmode = False

    def to_row(self):
        out = [self.label, str(self.color), str(self.index), str(self.mzrange[0]), str(self.mzrange[1]),
               str(self.zrange[0]), str(self.zrange[1]), str(self.sum), str(self.ht), str(self.sarray[0])
               , str(self.sarray[1]), str(self.sarray[2]), str(self.sarray[3])]
        return out


class ChromatogramContainer:
    def __init__(self):
        self.chromatograms = []

    def add_chromatogram(self, data, decondat=None, ccsdat=None, label=None, color="#000000", index=None, mzrange=None,
                         zrange=None, massrange=None, massmode=False, mode="DM", sarray=None):
        chrom = Chromatogram()
        chrom.chromdat = data
        chrom.decondat = decondat
        chrom.ccsdat = ccsdat
        chrom.sum = np.sum(data[:, 1])

        if label is None:
            label = ""
            if mzrange is not None:
                label += "m/z: " + str(round(mzrange[0])) + "-" + str(round(mzrange[1]))
            if zrange is not None:
                label += " z: " + str(round(zrange[0])) + "-" + str(round(zrange[1]))
            if sarray is not None:
                label += "Swoop m/z: " + str(round(sarray[0])) + " z: " + str(round(sarray[1]))

        chrom.label = label

        chrom.color = color

        if index is not None:
            chrom.index = index
        else:
            chrom.index = len(self.chromatograms)

        if mzrange is not None:
            chrom.mzrange = mzrange
            chrom.massmode = False

        if zrange is not None:
            chrom.zrange = zrange

        if massrange is not None:
            chrom.massrange = massrange
            chrom.massmode = True

        if sarray is not None:
            chrom.sarray = sarray
            chrom.swoopmode = True
            chrom.massmode = False

        chrom.massmode = massmode

        # If label already exists, replace it
        if label in [x.label for x in self.chromatograms]:
            for i, x in enumerate(self.chromatograms):
                if x.label == label:
                    self.chromatograms[i] = chrom
                    break
        else:
            # Add it if new
            if type(self.chromatograms) is np.ndarray:
                self.chromatograms = np.concatenate((self.chromatograms, [chrom]))
            else:
                self.chromatograms.append(chrom)

        return chrom

    def clear(self):
        self.chromatograms = []

    def to_array(self):
        arr = []
        for chrom in self.chromatograms:
            arr.append(chrom.to_row())
        return np.array(arr).astype(str)


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
