from __future__ import unicode_literals
import string
import math
# import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import Normalize
import numpy as np
from unidec import tools as ud
import pandas as pd
from io import StringIO

__author__ = 'Michael.Marty'


class Peak:
    """
    Class for a single peak. Contains all key parameters for describing and plotting the peak.
    """

    def __init__(self):
        """
        Initialize all parameters for the peak to defaults
        """
        self.mass = 0
        self.height = 0
        self.ccs = 0
        self.centroid = 0
        self.area = 0
        self.color = [1, 1, 1]
        self.label = ""
        self.marker = "."
        self.textmarker = "."
        self.ignore = 0
        self.match = 0
        self.matcherror = 0
        self.altmatches = []
        self.altmatcherrors = []
        self.numberalts = 0
        self.integral = 0
        self.integralrange = []
        self.mztab = []
        self.mztab2 = []
        self.stickdat = []
        self.kendricknum = 0
        self.kendrickdefect = 0
        self.kmass = 0
        self.score = 0
        # self.corrint = 0
        # self.correrr = 0
        self.mztabi = []
        self.massavg = 0
        self.masserr = 0
        # self.tval = 0
        self.peakmasses = []
        # self.fitmassavg = 0
        # self.fitmasserr = 0
        # self.fitarea = 0
        # self.fitareaerr = 0
        self.diff = 0
        self.extracts = []
        self.errorFWHM = 0
        self.intervalFWHM = [0, 0]
        self.badFWHM = False
        self.errormean = -1
        self.errorreplicate = 0
        self.avgcharge = 0
        self.zstack = []
        self.mzstack = []
        self.mscore = 0
        self.uscore = 0
        self.cs_score = 0
        self.rsquared = 0
        self.fscore = 0
        self.dscore = 0
        self.lscore = 0
        self.mdist = None
        self.zdist = None
        self.estimatedarea = 0
        self.index = 0
        self.sdnum = 0
        self.sdval = 0
        self.filename = ""
        self.filenumber = -1
        self.zs = []

    def line_out(self, type="Full", rounding=3):
        if type == "Full":
            outputs = [self.textmarker, self.mass, self.centroid, self.height, self.integral, self.match,
                       self.matcherror, self.label,
                       self.area, self.diff, self.avgcharge, self.dscore, self.errorFWHM, self.intervalFWHM[0],
                       self.intervalFWHM[1], self.errormean, self.errorreplicate, self.numberalts, self.altmatches]
        elif type == "Basic":
            outputs = [self.mass, self.height, self.integral]
        elif type == "FullFiles":
            outputs = [self.textmarker, self.mass, self.centroid, self.height, self.integral, self.match,
                       self.matcherror, self.label,
                       self.area, self.diff, self.avgcharge, self.dscore, self.errorFWHM, self.intervalFWHM[0],
                       self.intervalFWHM[1], self.errormean, self.errorreplicate, self.numberalts, self.altmatches,
                       self.filename, self.filenumber]
        else:
            outputs = [self.mass, self.height]
        outstring = ""
        for o in outputs:
            # Check if o is a list or array
            if isinstance(o, (list, np.ndarray)):
                for oo in o:
                    if rounding is not None:
                        try:
                            oo = round(oo, rounding)
                        except:
                            pass
                    outstring += str(oo) + " "
                outstring += "\t"
            else:
                if rounding is not None:
                    try:
                        o=round(o, rounding)
                    except:
                        pass
                outstring += str(o) + "\t"
        return outstring


class Peaks:
    """
    Class containing all useful data about peaks.

    The peaks themselves are of the Peak class and contained within the self.peaks list.
    """

    def __init__(self):
        """
        Initialize Peaks class and set empty values
        :return: None
        """
        self.peaks = []
        self.plen = 0
        self.changed = 0
        self.masses = []
        self.heights = []
        self.centroids = []
        self.areas = []
        self.fwhms = []
        self.convolved = False
        self.composite = None
        self.peakcolors = []
        self.markers = []
        self.colormap = []
        self.textmarkers = []
        self.marklen = 0
        self.massbins = 0
        self.norm = 1
        self.uniscore = 0

    def add_peaks(self, parray, massbins=0, scores_included=False):
        """
        Create peak objects from an array
        :param parray: N x 2 array containing (mass, height) of each peak.
        :param massbins: Describes the precision of the mass inputs by describing the bin size on the mass axis.
        :return: None
        """
        for p in parray:
            newpeak = Peak()
            newpeak.mass = p[0]
            newpeak.height = p[1]
            if scores_included:
                newpeak.dscore = p[2]
            self.peaks.append(newpeak)
        self.masses = np.array([p.mass for p in self.peaks])
        self.heights = np.array([p.height for p in self.peaks])
        self.plen = len(self.peaks)
        self.convolved = False
        self.composite = None
        self.massbins = massbins

    def merge_in_peaks(self, pks, filename=None, filenumber=None):
        for p in pks.peaks:
            if filename is not None:
                p.filename = filename
            if filenumber is not None:
                p.filenumber = filenumber
            self.peaks.append(p)
        return self


    def merge_isodec_pks(self, pks, massbins=0, config=None):
        masspks = pks.masses
        # Sort pks
        pks = sorted(masspks, key=lambda x: x.monoiso)

        for p in pks:
            newpeak = Peak()
            newpeak.mass = p.monoiso
            newpeak.height = p.totalintensity
            newpeak.mztab = np.transpose([p.mzs, p.mzints])
            newpeak.zs = p.zs
            newpeak.stickdat = p.isodists

            self.peaks.append(newpeak)
        self.masses = np.array([p.mass for p in self.peaks])
        self.heights = np.array([p.height for p in self.peaks])
        self.plen = len(self.peaks)

        if config is not None:
            if config.peaknorm == 1:
                self.norm = np.amax(self.heights)
            elif config.peaknorm == 2:
                self.norm = np.sum(self.heights)

            if config.peaknorm > 0:
                for p in self.peaks:
                    p.height = p.height / self.norm
                self.heights = np.array([p.height for p in self.peaks])

        self.convolved = False
        self.composite = None
        self.massbins = massbins

        self.default_params()



    def default_params(self, cmap="rainbow"):
        """
        Set default parameters for peaks, such as color, label, and marker
        :param cmap: Colormap from matplotlib.cm
        :return: None
        """
        if cmap[:2] == "b'":
            cmap = cmap[2:-1]
        try:
            cmap = str(cmap, encoding="utf-8")
        except:
            pass

        # self.colormap = cm.get_cmap(cmap, len(self.peaks))
        self.colormap = mpl.colormaps[cmap].resampled(len(self.peaks))
        if self.colormap is None:
            # self.colormap = cm.get_cmap(u"rainbow", len(self.peaks))
            self.colormap = mpl.colormaps[u"rainbow"].resampled(len(self.peaks))
        self.peakcolors = self.colormap(np.arange(len(self.peaks)))
        self.markers = ['o', 'v', '^', '>', 's', 'd', '*']
        self.textmarkers = ['\u25CB', '\u25BD', '\u25B3', '\u25B7', '\u25A2', '\u2662', '\u2606']
        self.marklen = len(self.markers)
        for i in range(0, len(self.peaks)):
            self.peaks[i].marker = self.markers[i % self.marklen]
            self.peaks[i].textmarker = self.textmarkers[i % self.marklen]
            self.peaks[i].color = self.peakcolors[i]
            if i >= 26:
                self.peaks[i].label = string.ascii_uppercase[i % 26] + str(int(math.floor(i / 26) + 1))
            else:
                self.peaks[i].label = string.ascii_uppercase[i % 26]

    def color_by_score(self, e=0):
        scores = np.array([p.dscore for p in self.peaks])
        # colormap = cm.get_cmap('RdYlGn')
        colormap = mpl.colormaps['RdYlGn']
        norm = Normalize(vmin=0, vmax=1)
        colors = colormap(norm(scores))
        self.peakcolors = colors
        for i in range(0, len(self.peaks)):
            self.peaks[i].color = self.peakcolors[i]

    def get_mass_defects(self, kendrickmass, mode=0):
        """
        Get the mass defects and mass number for each peak
        :param kendrickmass: Kendrick reference mass
        :param mode: Select range of defects 1=(0,1), 0=(-0.5,0.5)
        :return: None
        """
        for p in self.peaks:
            p.kmass = p.mass / float(kendrickmass)
            if mode == 1:
                p.kendricknum = np.floor(p.kmass)
                p.kendrickdefect = p.kmass - np.floor(p.kmass)
            else:
                p.kendricknum = np.round(p.kmass)
                p.kendrickdefect = p.kmass - np.round(p.kmass)

    def get_bool(self):
        boo1 = []
        for p in self.peaks:
            if p.ignore == 0:
                boo1.append(True)
            else:
                boo1.append(False)
        return np.array(boo1)

    def diffs_from(self, target):
        for p in self.peaks:
            p.diff = p.mass - target
        return np.array([p.diff for p in self.peaks])

    def diffs_consecutive(self):
        b1 = self.get_bool()
        pmasses = np.array([p.mass for p in self.peaks])[b1]
        peakdiff = np.zeros(len(pmasses))
        peakdiff[1:] = np.diff(pmasses)
        for i, p in enumerate(self.peaks):
            p.diff = 0
        for i, p in enumerate(np.array(self.peaks)[b1]):
            p.diff = peakdiff[i]
        return np.array([p.diff for p in self.peaks])

    def integrate(self, data, lb=None, ub=None):
        self.areas = []
        for p in self.peaks:
            if lb is None:
                if len(p.intervalFWHM) == 2:
                    lb = (p.intervalFWHM[0] - p.mass) * 2
                else:
                    print("ERROR IN INTEGRATION. Need to calc FWHM first.")
                    lb = 0
            if ub is None:
                if len(p.intervalFWHM) == 2:
                    ub = (p.intervalFWHM[1] - p.mass) * 2
                else:
                    print("ERROR IN INTEGRATION. Need to calc FWHM first.")
                    ub = 0
            p.integralrange = [p.mass + lb, p.mass + ub]
            p.integral = ud.integrate(data, p.integralrange[0], p.integralrange[1])[0]
            self.areas.append(p.integral)
        self.areas = np.array(self.areas)

    def auto_format(self):

        colors = np.array([[[1, 0, 0], [1, 0, 0]], [[0, 0.7, 0], [1, 1, 0]], [[0, 0.5, 1], [1, 0, 1]]])

        newmarkers = [[0, 0], [1, 5], [4, 6]]

        # note = [[0, 1], [1, 1], [2, 1]]
        # self.markers = ['o', 'v', '^', '>', 's', 'd', '*']
        # self.textmarkers = [u'\u25CB', u'\u25BD', u'\u25B3', u'\u25B7', u'\u25A2', u'\u2662', u'\u2606']

        for p in self.peaks:
            n = p.label
            splits = n.split("[")

            try:
                n1 = int(splits[0])
            except:
                n1 = 0

            try:
                n2 = int(splits[1].split("]")[1])
            except:
                n2 = 0
            try:
                newcolor = colors[n1][n2]
                newmarker = self.markers[newmarkers[n1][n2]]
                newtextmarker = self.textmarkers[newmarkers[n1][n2]]
            except:
                newcolor = [0.5, 0.5, 1]
                newmarker = self.markers[6]
                newtextmarker = self.textmarkers[6]

            p.color = newcolor
            p.marker = newmarker
            p.textmarker = newtextmarker
            # print n1, n2, newcolor, newmarker, newtextmarker

    def copy(self, type="Full"):
        if type == "Full":
            outstring = "Symbol\tMass\tCentroid\tHeight\tIntegral\tMatch\tMatcherror\tLabel" \
                        "\tFit Area\tDiff\tAvgcharge\tDscore\tFWHM\tLowValFWHM\tHighValFWHM\tErrorMean\t" \
                        "ErrorReplicate\tNumMatches\tAltMatches\n"
        elif type == "Basic":
            outstring = "Mass\tHeight\tIntegral\n"
        elif type == "FullFiles":
            outstring = "Symbol\tMass\tCentroid\tHeight\tIntegral\tMatch\tMatcherror\tLabel" \
                        "\tFit Area\tDiff\tAvgcharge\tDscore\tFWHM\tLowValFWHM\tHighValFWHM\tErrorMean\t" \
                        "ErrorReplicate\tNumMatches\tAltMatches\tFileName\tFileNumber\n"
        # print("Columns:", outstring)

        for p in self.peaks:
            outstring += p.line_out(type=type) + "\n"
        return outstring

    def to_df(self, type="Full", drop_zeros=True):
        outstring = self.copy(type=type)
        # print(outstring)
        df = pd.read_csv(StringIO(outstring), sep="\t", index_col=False, na_values="", dtype=str)
        df.fillna("", inplace=True)
        if drop_zeros:
            df = df.loc[:, (df != "0").any(axis=0)]
            df = df.loc[:, (df != "-1").any(axis=0)]
            df = df.loc[:, (df != "").any(axis=0)]
        return df
