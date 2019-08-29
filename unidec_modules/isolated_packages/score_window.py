import numpy as np
import os
import wx
import wx.grid
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import unidec_modules.unidectools as ud
from copy import deepcopy
import unidec

luminance_cutoff = 135
white_text = wx.Colour(250, 250, 250)
black_text = wx.Colour(0, 0, 0)


def scr(score):
    s = int(np.round(score * 100))
    return str(s)

def scr2(p):
    scores = [p.dscore, p.uni_score, p.pca_score, p.z_score, p.fwhm_score]
    strings = [scr(s) for s in scores]
    out = ""
    for s in strings:
        out+=s+","
    return out

def score_plots(eng):
    import matplotlib.pyplot as plt
    ztab = eng.data.ztab
    for k, p in enumerate(eng.pks.peaks):
        m = p.zstack[0]
        ints = p.zstack[1:]
        summ = np.sum(ints, axis=0)
        sumz = np.sum(ints, axis=1)
        sumz /= np.amax(sumz)

        plt.figure()
        plt.subplot(131)
        plt.title(
            "Combined Score: " + str(round(p.dscore * 100, 2)) + "\nCharge Score: " + str(round(p.z_score * 100, 2))
            )
        plt.plot(ztab, sumz)
        plt.plot(p.pca_zdist[:, 0], p.pca_zdist[:, 1], color="k")
        plt.subplot(132)
        for ival in ints:
            mdat = ival / np.amax(ints)
            plt.plot(m, mdat)
            plt.plot(p.pca_mdist[:, 0], p.pca_mdist[:, 1]/np.sum(p.pca_mdist[:,1])*np.sum(mdat), color="k")
        plt.title("Mass:" + str(p.mass) + "\nPCA Score: " + str(round(p.pca_score * 100, 2)))

        plt.subplot(133)
        for i, z in enumerate(ztab):
            v = p.mzstack[i]
            x = v[:, 0]
            y = v[:, 1]
            zvals = v[:, 2]
            plt.title(
                "Uniqueness : " + str(round(p.uni_score * 100, 2)) + "\nR Squared: " + str(round(p.rsquared * 100, 2))
                )
            if True:
                try:
                    x = x - np.amin(x)
                    x = x / np.amax(x)
                except:
                    pass
            plt.plot(x, y)  # , color="k")
            plt.plot(x, zvals, color="k")

    plt.show()

class ScoreListCtrl(wx.grid.Grid):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize, style=wx.LC_REPORT):
        wx.grid.Grid.__init__(self, parent, id_value, pos, size, style=style)
        self.parent = parent

    def populate(self, pks):
        self.ClearGrid()

        l = len(pks.peaks)

        collabels = ["Mass", "Intensity", "DScore", "Uniqueness/Fit", "Peak Shape",
                     "Charge Dist.", "FWHM Penalties"]

        cl = len(collabels)
        self.CreateGrid(l + 1, cl)
        w = 120
        for i, c in enumerate(collabels):
            self.SetCellValue(0, i, c)
            self.SetColSize(i, width=w)
            self.SetCellFont(0, i, wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        self.SetColSize(0, width=75)
        self.SetColSize(1, width=75)

        colormap = cm.get_cmap('RdYlGn')
        norm = Normalize(vmin=0, vmax=1)

        for index, p in enumerate(pks.peaks):
            i = index + 1
            self.SetCellValue(i, 0, str(p.mass))
            self.SetCellFont(i, 0, wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
            self.SetCellValue(i, 1, str(p.height))

            scores = [p.dscore, p.uni_score, p.pca_score, p.z_score, p.fwhm_score]
            cindex = 2
            for s in scores:
                self.SetCellValue(i, cindex, scr(s))

                color = colormap(norm(s))

                color = wx.Colour(int(round(color[0] * 255)), int(round(color[1] * 255)), int(round(color[2] * 255)),
                                  alpha=255)

                self.SetCellBackgroundColour(i, cindex, color)
                luminance = ud.get_luminance(color, type=2)
                if luminance < luminance_cutoff:
                    self.SetCellTextColour(i, cindex, white_text)
                else:
                    self.SetCellTextColour(i, cindex, black_text)

                cindex += 1

    def clear_list(self):
        self.ClearGrid()


class ScoreFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Score Window", size=(1000, 500))  # ,size=(-1,-1))
        self.parent = parent
        self.type = 1
        self.elevation = None
        self.azimuth = None

        self.panel = wx.Panel(self)

        self.listctrl = ScoreListCtrl(self, wx.ID_ANY, size=(1000, 500))

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.listctrl, 1, wx.LEFT | wx.TOP | wx.GROW)

        self.panel.SetSizer(self.sizer)
        self.panel.Fit()
        self.Show(True)

    def populate(self, pks):
        self.listctrl.populate(pks)

if __name__ == "__main__":

    if True:
        directory = "C:\\Data\\AlgorithmTesting"
        # directory = "C:\\Python\\UniDec3\\unidec_bin\\Example Data"
        os.chdir(directory)

        file = "ADH.txt"

        eng = unidec.UniDec()
        eng.open_file(file, directory)
        if True:
            eng.unidec_imports(everything=True)
        else:
            eng.process_data()
            eng.config.massub = 200000
            eng.config.masslb = 50000
            eng.run_unidec()
        eng.config.peakthresh = 0.025
        eng.pick_peaks()
        # eng.pks_pca()
        eng.dscore()

    app = wx.App()
    panel = ScoreFrame(None)
    panel.populate(eng.pks)
    app.MainLoop()
