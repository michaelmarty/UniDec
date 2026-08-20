import numpy as np
import os
import wx
import wx.grid
#import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.colors import Normalize
import unidec.tools as ud
from copy import deepcopy
from unidec import engine

luminance_cutoff = 135
white_text = wx.Colour(250, 250, 250)
black_text = wx.Colour(0, 0, 0)


def scr(score):
    s = int(np.round(score * 100))
    return str(s)


def scr2(p):
    scores = [p.dscore, p.uscore, p.mscore, p.cs_score, p.fscore]
    strings = [scr(s) for s in scores]
    out = ""
    for s in strings:
        out += s + ","
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
            "Combined Score: " + str(round(p.dscore * 100, 2)) + "\nCharge Score: " + str(round(p.cs_score * 100, 2))
        )
        plt.plot(ztab, sumz)
        plt.plot(p.zdist[:, 0], p.zdist[:, 1], color="k")
        plt.subplot(132)
        for ival in ints:
            mdat = ival / np.amax(ints)
            plt.plot(m, mdat)
            plt.plot(p.mdist[:, 0], p.mdist[:, 1] / np.sum(p.mdist[:, 1]) * np.sum(mdat), color="k")
        plt.title("Mass:" + str(p.mass) + "\nPeak Shape Score: " + str(round(p.mscore * 100, 2)))

        plt.subplot(133)
        for i, z in enumerate(ztab):
            v = p.mzstack[i]
            x = v[:, 0]
            y = v[:, 1]
            zvals = v[:, 2]
            plt.title(
                "Uniqueness : " + str(round(p.uscore * 100, 2)) + "\nFWHM: " + str(round(p.fscore * 100, 2))
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


def score_plots2(eng, save=False):
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams['ps.useafm'] = True
    rcParams['ps.fonttype'] = 42
    rcParams['pdf.fonttype'] = 42

    # colormap = cm.get_cmap('RdYlGn')
    colormap = mpl.colormaps["RdYlGn"]
    norm = Normalize(vmin=0, vmax=1)

    plt.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=[4, 3])
    p1 = plt.subplot(211)
    plt.plot(eng.data.data2[:, 0], eng.data.data2[:, 1] * 100, color="k")
    plt.xlabel(r"$\it{m/z}$")
    p1.get_yaxis().set_ticks([0, 50, 100])
    p1.get_yaxis().set_ticklabels(["0", '%', "100"])
    p1.spines['top'].set_visible(False)
    p1.spines['right'].set_visible(False)
    plt.ylim(0, 100)
    plt.xlim(np.amin(eng.data.data2[:, 0]), np.amax(eng.data.data2[:, 0]))
    p2 = plt.subplot(212)
    plt.plot(eng.data.massdat[:, 0] / 1000., eng.data.massdat[:, 1], color="k")
    plt.xlabel("Mass (kDa)")
    p2.get_yaxis().set_ticks([0, 50, 100])
    p2.get_yaxis().set_ticklabels(["0", '%', "100"])
    p2.spines['top'].set_visible(False)
    p2.spines['right'].set_visible(False)
    plt.ylim(0, 100)
    xmin = np.amin(eng.data.massdat[:, 0]) / 1000
    xmax = np.amax(eng.data.massdat[:, 0]) / 1000
    plt.xlim(xmin, xmax)
    offset = 0.2 * np.amax(eng.data.massdat[:, 1])
    for k, p in enumerate(eng.pks.peaks):
        text = str(int(round(p.dscore * 100)))
        color = colormap(norm(p.dscore))
        lum = ud.get_luminance(color)
        if lum > luminance_cutoff / 255.:
            #print(p.mass, lum, color)
            color = np.array(color) * 0.8
            color[3] = 1
        p2.text(p.mass / 1000., p.height + offset, text, horizontalalignment="center",
                verticalalignment="top", color=color)
        pass

    text = "UniScore=" + str(int(round(eng.pks.uniscore * 100)))
    color = colormap(norm(eng.pks.uniscore))
    p2.text(xmax, 100, text, horizontalalignment="right",
            verticalalignment="top", color=color)
    plt.tight_layout()

    if save:
        fig.savefig("Fig1.pdf")
        fig.savefig("Fig1.png")
    # plt.show()
    return fig


class ScoreListCtrl(wx.grid.Grid):
    def __init__(self, parent, id_value, pos=wx.DefaultPosition, size=wx.DefaultSize):  # , style=wx.LC_REPORT):
        wx.grid.Grid.__init__(self, parent, id_value, pos, size)  # , style=style)
        self.parent = parent
        self.Bind(wx.EVT_KEY_DOWN, self.OnKey)

    def OnKey(self, event):
        # If Ctrl+c is pressed...
        if event.ControlDown() and event.GetKeyCode() == 67:
            self.copy()

    def populate(self, pks):
        self.ClearGrid()

        l = len(pks.peaks)

        collabels = ["Mass", "Intensity", "DScore", "Uniqueness/Fit", "Peak Shape",
                     "Charge Dist.", "FWHM Penalties", "Limited Score (DScore without CSScore)"]

        cl = len(collabels)
        self.CreateGrid(l + 1, cl)
        w = 120
        for i, c in enumerate(collabels):
            self.SetCellValue(0, i, c)
            self.SetColSize(i, width=w)
            self.SetCellFont(0, i, wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        self.SetColSize(0, width=75)
        self.SetColSize(1, width=75)

        # colormap = cm.get_cmap('RdYlGn')
        colormap = mpl.colormaps['RdYlGn']
        norm = Normalize(vmin=0, vmax=1)

        for index, p in enumerate(pks.peaks):
            i = index + 1
            self.SetCellValue(i, 0, str(p.mass))
            self.SetCellFont(i, 0, wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
            self.SetCellValue(i, 1, str(p.height))

            scores = [p.dscore, p.uscore, p.mscore, p.cs_score, p.fscore, p.lscore]
            scores = deepcopy(scores)
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

    def copy(self):
        """Copies the current range of select cells to clipboard.
        """
        # Get number of copy rows and cols
        if self.GetSelectionBlockTopLeft() == []:
            rowstart = self.GetGridCursorRow()
            colstart = self.GetGridCursorCol()
            rowend = rowstart
            colend = colstart
        else:
            rowstart = self.GetSelectionBlockTopLeft()[0][0]
            colstart = self.GetSelectionBlockTopLeft()[0][1]
            rowend = self.GetSelectionBlockBottomRight()[0][0]
            colend = self.GetSelectionBlockBottomRight()[0][1]

        self.crows = rowend - rowstart + 1
        self.ccols = colend - colstart + 1

        # data variable contains text that must be set in the clipboard
        data = ''
        # For each cell in selected range append the cell value
        # in the data variable Tabs '\t' for cols and '\n' for rows
        for r in range(self.crows):
            for c in range(self.ccols):
                data += str(self.GetCellValue(rowstart + r, colstart + c))
                if c < self.ccols - 1:
                    data += '\t'
            data += '\n'

        # Create text data object
        clipboard = wx.TextDataObject()

        # Set data object value
        clipboard.SetText(data)

        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Error")


class ScoreFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Score Window", size=(1017, 500))  # ,size=(-1,-1))
        self.parent = parent
        self.type = 1
        self.elevation = None
        self.azimuth = None

        self.panel = wx.Panel(self)

        self.listctrl = ScoreListCtrl(self.panel, wx.ID_ANY, size=(1000, 500))

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.listctrl, 0, wx.LEFT | wx.TOP | wx.EXPAND)

        self.panel.SetSizer(self.sizer)
        self.panel.Fit()
        self.Show(True)

    def populate(self, pks):
        self.listctrl.populate(pks)


def score_window(eng):
    app = wx.App()
    panel = ScoreFrame(None)
    panel.populate(eng.pks)
    app.MainLoop()


if __name__ == "__main__":

    if True:
        directory = "C:\\Data\\AlgorithmTesting"
        # directory = "C:\\Python\\UniDec3\\bin\\Example Data"
        os.chdir(directory)

        file = "ADH.txt"
        file = "0.txt"

        eng = engine.UniDec()
        eng.open_file(file, directory)
        if True:
            eng.unidec_imports(everything=True)
        else:
            eng.process_data()
            eng.config.massub = 200000
            eng.config.masslb = 50000
            eng.run_unidec()
        # eng.config.peakthresh = 0.025
        eng.pick_peaks()
        # eng.pks_pca()
        eng.dscore()

    app = wx.App()
    panel = ScoreFrame(None)
    panel.populate(eng.pks)
    app.MainLoop()
