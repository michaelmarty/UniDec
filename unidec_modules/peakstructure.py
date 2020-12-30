from __future__ import unicode_literals
import string
import math
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np
from unidec_modules import unidectools as ud

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
        self.area = ""
        self.color = [1, 1, 1]
        self.label = ""
        self.marker = "."
        self.textmarker = "."
        self.ignore = 0
        self.match = 0
        self.matcherror = 0
        self.integral = 0
        self.integralrange = []
        self.mztab = []
        self.mztab2 = []
        self.stickdat = []
        self.kendricknum = 0
        self.kendrickdefect = 0
        self.kmass = 0
        self.score = 0
        #self.corrint = 0
        #self.correrr = 0
        self.mztabi = []
        self.massavg = 0
        self.masserr = 0
        #self.tval = 0
        self.peakmasses = []
        #self.fitmassavg = 0
        #self.fitmasserr = 0
        #self.fitarea = 0
        #self.fitareaerr = 0
        self.diff = 0
        self.extracts = []
        self.errorFWHM = 0
        self.intervalFWHM = [0, 0]
        self.badFWHM=False
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

    def line_out(self, type="Full"):
        if type == "Full":
            outputs = [self.mass, self.centroid, self.height, self.area, self.match, self.matcherror,
                       self.integral, self.diff, self.avgcharge, self.dscore, self.lscore]
        elif type == "Basic":
            outputs = [self.mass, self.height, self.integral]
        else:
            outputs = [self.mass, self.height]
        outstring = ""
        for o in outputs:
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
        self.plen = len(self.peaks)
        self.convolved = False
        self.composite = None
        self.massbins = massbins


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

        self.colormap = cm.get_cmap(cmap, len(self.peaks))
        if self.colormap is None:
            self.colormap = cm.get_cmap(u"rainbow", len(self.peaks))
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
        colormap = cm.get_cmap('RdYlGn')
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
            print("Columns: mass, centroid, height, area, match, matcherror, integral, diff, avgcharge, dscore, lscore")
        if type == "Basic":
            print("Columns: mass, height, integral")

        outstring = ""
        for p in self.peaks:
            outstring += p.line_out(type=type) + "\n"
        return outstring

    '''
    def score_peaks(self, thresh=0, ci=0.99):
        """
        For each peak, assign a score of the fractional number of charge states observed.
        :param thresh: Optional threshold to define a noise floor.
        :return: None
        """
        for p in self.peaks:
            boo1 = p.mztab[:, 1] > thresh
            boo2 = p.mztab2[:, 1] > 0
            booi = p.mztabi[:, 0] > 0
            boo3 = np.all([boo1, boo2, booi], axis=0)
            try:
                p.score = np.sum((np.ones_like(p.mztab[boo3, 1]) - np.clip(
                    np.abs((p.mztab[boo3, 1] - p.mztab2[boo3, 1])) / p.mztab2[boo3, 1], 0, 1)))
            except ZeroDivisionError:
                p.score = 0

            try:
                dof = p.score - 1
                if dof < 0:
                    tval = 0
                else:
                    tval = ud.get_tvalue(dof, ci=ci)
                p.massavg = np.average(p.peakmasses)
                p.masserr = tval * np.std(p.peakmasses, ddof=1) / np.sqrt(p.score)
                p.tval = tval
            except (ValueError, ZeroDivisionError, RuntimeWarning, RuntimeError):
                p.massavg = 0
                p.masserr = 0
                p.tval = 0
    '''

