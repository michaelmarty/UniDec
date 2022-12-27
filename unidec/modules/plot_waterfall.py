from matplotlib.collections import PolyCollection
from matplotlib import colors as mcolors
import numpy as np
import wx
from copy import deepcopy
from unidec.modules.PlottingWindow import PlottingWindowBase

__author__ = 'Michael.Marty'


def cc(arg):
    return mcolors.to_rgba(arg, alpha=0.6)


class Waterfall3DPlot(PlottingWindowBase):
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
        PlottingWindowBase.__init__(self, *args, **kwargs)
        self._axes = [0.05, 0.1, 0.9, 0.9]
        self.cids = []


    def waterfall(self, xgrid, ygrid, zarray, colors=None, xlabel="Mass", ylabel = ""):

        self.xlabel = xlabel
        if xlabel == "Mass":
            try:
                self.kda_test(xgrid[0])
                if self.kdnorm==1000:
                    xgrid2 = []
                    for x in xgrid:
                        x = x/self.kdnorm
                        xgrid2.append(x)
                    xgrid = xgrid2
            except:
                self.kdnorm = 1

        self.clear_plot("nopaint")
        self.subplot1 = self.figure.add_axes(self._axes, projection='3d', proj_type="ortho")
        verts = []
        for i, z in enumerate(zarray):
            verts.append(list(zip(xgrid[i], ygrid[i])))
        # ax = self.figure.gca(projection='3d')

        if colors is None:
            colors = cc('r')
            # wc = cc('w')

        poly = PolyCollection(verts, edgecolors=colors, facecolors=colors)
        #poly.set_facecolor((1,1,1,0))
        poly.set_alpha(0.7)
        labels=None
        try:
            float(zarray[0])
        except:
            labels=zarray
            zarray = np.arange(len(zarray))


        self.subplot1.add_collection3d(poly, zs=zarray, zdir='y')

        xgrid = np.hstack(xgrid)
        ygrid = np.hstack(ygrid)

        self.subplot1.set_xlabel(self.xlabel)
        self.subplot1.set_xlim3d(np.amin(xgrid), np.amax(xgrid))
        self.subplot1.set_ylabel(ylabel)
        self.subplot1.set_ylim3d(np.amax(zarray), np.amin(zarray))
        #self.subplot1.set_zlabel('Intensity')
        self.subplot1.set_zlim3d(np.amin(ygrid), np.amax(ygrid))
        self.subplot1.get_zaxis().set_ticks([0, np.amax(ygrid) / 2, np.amax(ygrid)])
        self.subplot1.get_zaxis().set_ticklabels(["0", '%', "100"])

        if labels is not None:
            self.subplot1.get_yaxis().set_ticks(zarray)
            self.subplot1.get_yaxis().set_ticklabels(labels)

        self.subplot1.xaxis.pane.fill = False
        self.subplot1.yaxis.pane.fill = False
        self.subplot1.zaxis.pane.fill = False

        # Now set color to white (or whatever is "invisible")
        self.subplot1.xaxis.pane.set_edgecolor('w')
        self.subplot1.yaxis.pane.set_edgecolor('w')
        self.subplot1.zaxis.pane.set_edgecolor('w')

        # Bonus: To get rid of the grid as well:
        self.subplot1.grid(False)

        self.canvas = self.figure.canvas
        for cid in self.cids:
            self.canvas.mpl_disconnect(cid)
        self.cids.append(self.canvas.mpl_connect('button_release_event', self.onmove))
        try:
            self.repaint()
        except MemoryError:
            print("Memory Error: Not updating 2D plot")
        self.flag = True

    def set_angle(self, elevation, azimuth):
        self.subplot1.view_init(elevation, azimuth)
        self.repaint()

    def get_angle(self):
        a = deepcopy(self.subplot1.azim)
        e = deepcopy(self.subplot1.elev)
        return a, e

    def onmove(self, e=None):
        self.frame.onmove()


class WaterfallFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Waterfall Plots", size=(1000, 500))  # ,size=(-1,-1))
        self.parent = parent
        self.type = 1
        self.elevation = None
        self.azimuth = None

        self.panel = wx.Panel(self)
        self.plot = Waterfall3DPlot(self.panel, frame=self)
        self.ctltype = wx.RadioBox(self.panel, label="Data Type", choices=["m/z", "Mass"])
        self.ctltype.SetSelection(self.type)
        self.ctlelevation = wx.TextCtrl(self.panel, value="")
        self.ctlazimuth = wx.TextCtrl(self.panel, value="")
        self.replotbutton = wx.Button(self.panel, -1, "Replot")
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.plot, 1, wx.LEFT | wx.TOP | wx.GROW)

        s2 = wx.BoxSizer(wx.HORIZONTAL)
        s2.Add(self.ctltype, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(wx.StaticText(self.panel, label="  Elevation: "), flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(self.ctlelevation, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(wx.StaticText(self.panel, label="  Azimuth: "), flag=wx.ALIGN_CENTER_VERTICAL)
        s2.Add(self.ctlazimuth, 0, flag=wx.ALIGN_CENTER_VERTICAL)

        s2.Add(self.replotbutton, 0, flag=wx.ALIGN_CENTER_VERTICAL)

        self.sizer.Add(s2, 0)

        #self.Bind(wx.EVT_TEXT, self.change, self.ctlelevation)
        #self.Bind(wx.EVT_TEXT, self.change, self.ctlazimuth)
        self.Bind(wx.EVT_RADIOBOX, self.change_type, self.ctltype)
        self.Bind(wx.EVT_BUTTON, self.change, self.replotbutton)

        self.panel.SetSizer(self.sizer)
        self.panel.Fit()
        self.Show(True)

    def draw(self):
        xs = []
        zs = [0.0, 1.0, 2.0, 3.0]
        ys = []
        x = np.arange(0, 10, 0.4)*10000
        for z in zs:
            y = np.random.rand(len(x))
            y[0], y[-1] = 0, 0
            ys.append(y)
            xs.append(x)
        self.plot.waterfall(xs, ys, zs)
        self.onmove()

    def make_plot(self, data=None):
        if data is not None:
            self.data = data

        xgrid = []
        ygrid = []
        colors = []
        zarray = self.data.var1
        label = "Mass"
        for s in self.data.spectra:
            if self.type == 1:
                mdat = s.massdat
            elif self.type == 0:
                mdat = s.data2
                label = "m/z"
            m = mdat[:, 0]
            i = mdat[:, 1]
            xgrid.append(m)
            ygrid.append(i)
            colors.append(s.color)

        ylabel= self.data.v1name

        self.plot.waterfall(xgrid, ygrid, zarray, colors=colors, xlabel=label, ylabel=ylabel)
        self.onmove()

    def change(self, e=None):
        elev = self.ctlelevation.GetValue()
        try:
            elev = float(elev)
        except:
            elev = None
        self.elevation = elev

        azi = self.ctlazimuth.GetValue()
        try:
            azi = float(azi)
        except:
            azi = None
        self.azimuth = azi

        # self.make_plot()
        self.plot.set_angle(self.elevation, self.azimuth)
        self.plot.get_angle()

    def change_type(self, e=None):
        self.type = self.ctltype.GetSelection()
        self.make_plot()

    def onmove(self, e=None):
        a, e = self.plot.get_angle()
        a = round(a)
        e = round(e)
        self.ctlazimuth.SetValue(str(a))
        self.ctlelevation.SetValue(str(e))
        pass


if __name__ == "__main__":
    app = wx.App()
    panel = WaterfallFrame(None)
    panel.draw()
    app.MainLoop()
