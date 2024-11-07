from matplotlib import colors as mcolors
import numpy as np
import wx
from unidec.metaunidec.mudeng import MetaUniDec
from unidec.modules import PlottingWindow, unidecstructure

__author__ = 'Michael.Marty'


def cc(arg):
    return mcolors.to_rgba(arg, alpha=0.6)


class PlatePlot(PlottingWindow.PlottingWindowBase):
    """
    Plotting class for 3D Waterfall plots.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize the plotting window. Specify axes. Specify that the plot does not have normal tick marks to repaint.
        :param args:
        :param kwargs:
        :return: Waterfall Object
        """
        if "frame" in kwargs:
            self.frame = kwargs['frame']
            del kwargs['frame']
        PlottingWindow.PlottingWindowBase.__init__(self, *args, **kwargs)
        self._axes = [0.05, 0.1, 0.9, 0.9]
        self.cids = []

    def make_plate(self, type=96):
        self.clear_plot("nopaint")
        if type == 96:
            self.xlen = 12
            self.xticks = np.arange(self.xlen)
            self.xlabels = [str(i + 1) for i in self.xticks]
            self.ylen = 8
            self.yticks = np.arange(self.ylen)
            self.ylabels = np.array([chr(65 + i) for i in self.yticks])[::-1]

            self.subplot1 = self.figure.add_axes(self._axes)

            self.subplot1.get_xaxis().set_ticks(self.xticks)
            self.subplot1.get_xaxis().set_ticks(np.arange(self.xlen + 1) - 0.5, minor=True)
            self.subplot1.get_xaxis().set_ticklabels(self.xlabels)

            self.subplot1.get_yaxis().set_ticks(self.yticks)
            self.subplot1.get_yaxis().set_ticks(np.arange(self.ylen + 1) - 0.5, minor=True)
            self.subplot1.get_yaxis().set_ticklabels(self.ylabels)


            self.subplot1.grid(which="minor")

            self.canvas = self.figure.canvas
            for cid in self.cids:
                self.canvas.mpl_disconnect(cid)
            self.cids.append(self.canvas.mpl_connect('button_release_event', self.onmove))
            try:
                self.repaint()
            except MemoryError:
                print("Memory Error: Not updating 2D plot")
            self.flag = True

    def populate(self, positions, colors=None):
        for i, p in enumerate(positions):
            letter = p[0]
            number = int(p[1:])
            xindex = number - 1
            yindex = -1
            for j, l in enumerate(self.ylabels):
                if letter == l:
                    yindex = j
            if colors is not None:
                color= colors[i]
            else:
                color="k"
            self.subplot1.plot(xindex, yindex, marker="o", markersize=20, color=color)
            self.subplot1.set_ylim(np.amin(self.yticks) - 0.5, np.amax(self.yticks) + 0.5)
            self.subplot1.set_xlim(np.amin(self.xticks) - 0.5, np.amax(self.xticks) + 0.5)
            self.repaint(setupzoom=False)

    def onmove(self, e=None):
        try:
            x = round(e.xdata)
            y = round(e.ydata)
            xval = self.xlabels[x]
            yval = self.ylabels[y]
            clickpos = yval + xval
            print(clickpos)
        except:
            print("Click. Not on Plot")
            clickpos = "A1"

        self.frame.onmove(clickpos)


class PlateFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Plate Plot", size=(1000, 600))  # ,size=(-1,-1))
        self.parent = parent
        self.type = 0
        self.position = "A1"
        self.data = None
        self.config = unidecstructure.UniDecConfig()
        self.config.initialize()

        self.panel = wx.Panel(self)
        self.plateplot = PlatePlot(self.panel, frame=self)
        self.plot1 = PlottingWindow.Plot1d(self.panel)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        s1 = wx.BoxSizer(wx.HORIZONTAL)
        s1.Add(self.plateplot, 1, wx.LEFT | wx.TOP | wx.GROW)
        s1.Add(self.plot1, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.sizer.Add(s1, 1)

        self.ctltype = wx.RadioBox(self.panel, label="Data Type", choices=["m/z", "Mass"])
        self.ctltype.SetSelection(self.type)
        self.replotbutton = wx.Button(self.panel, -1, "Replot")
        s2 = wx.BoxSizer(wx.HORIZONTAL)

        s2.Add(self.ctltype, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        '''
        self.ctlelevation = wx.TextCtrl(self.panel, value="")
        self.ctlazimuth = wx.TextCtrl(self.panel, value="")
        s2.Add(wx.StaticText(self.panel, label="  Elevation: "), flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(self.ctlelevation, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(wx.StaticText(self.panel, label="  Azimuth: "), flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(self.ctlazimuth, 0, flag=wx.ALIGN_CENTER_VERTICAL)'''

        s2.Add(self.replotbutton, 0, flag=wx.ALIGN_CENTER_VERTICAL)

        self.sizer.Add(s2, 0)

        self.Bind(wx.EVT_RADIOBOX, self.change_type, self.ctltype)
        self.Bind(wx.EVT_BUTTON, self.replot, self.replotbutton)

        self.panel.SetSizer(self.sizer)
        self.panel.Fit()
        self.Show(True)

        self.draw()

    def draw(self):
        self.plateplot.make_plate()

    def load_eng(self, eng):
        print(eng.data.var1)
        self.positions = eng.data.var1
        for i in range(len(self.positions)):
            pos = self.positions[i]
            digits = ''.join(filter(str.isdigit, pos))
            letters = ''.join([i for i in pos if not i.isdigit()])
            if len(letters) == 2:
                letters = letters[1]
            newpos = letters + digits
            print("Changing:", pos, "to", newpos)
            self.positions[i] = newpos
        print(self.positions)
        self.position = self.positions[0]
        self.data = eng.data
        self.config = eng.config
        self.colors = []
        for i, s in enumerate(self.data.spectra):
            self.colors.append(s.color)

        self.plateplot.populate(self.positions, self.colors)
        self.make_plot()



    def make_plot(self, data=None):
        if data is not None:
            self.data = data

        for i, s in enumerate(self.data.spectra):
            if self.positions[i] == self.position:
                if self.type == 1:
                    mdat = s.massdat
                    label = "Mass"
                elif self.type == 0:
                    mdat = s.data2
                    label = "m/z"
                self.plot1.plotrefreshtop(mdat[:, 0], mdat[:, 1], self.position, label,
                                          "Intensity", "", config=self.config, color=s.color)
                self.plot1.add_title(self.position)

    def replot(self, e=None):
        self.make_plot()

    def change_type(self, e=None):
        self.type = self.ctltype.GetSelection()
        self.make_plot()

    def onmove(self, e=None):
        self.position = e
        print("click", self.position)
        self.replot()
        pass


if __name__ == "__main__":
    app = wx.App()
    panel = PlateFrame(None)
    panel.draw()

    eng = MetaUniDec()
    path = "C:\\Data\\HTS_Sharon\\20220404-5_sequence_Shortlist.hdf5"
    eng.open(path)
    panel.load_eng(eng)

    app.MainLoop()
