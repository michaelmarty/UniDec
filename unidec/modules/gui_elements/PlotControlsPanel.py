import matplotlib as mpl
import wx
from unidec.modules.PlottingWindow import Plot1d
from unidec.modules.unidecstructure import UniDecConfig

mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42

plot_params = {"font.size": 12, "xtick.labelsize": 12, "ytick.labelsize": 12, "lines.linewidth": 1.5,
               "lines.linestyle": "-", "font.family": "Arial", "axes.linewidth": 1.5, "axes.labelpad": 4.0,
               "axes.xmargin": 0.05, "axes.ymargin": 0.05,
               }


# Export plot_params to .mplstyle file
def export_mplstyle(params, filename="unidec.mplstyle"):
    with open(filename, "w") as f:
        for p in params:
            f.write(p + ": " + str(params[p]) + "\n")


# Plot Controls Panel
class PlotControlsPanel(wx.Panel):
    def __init__(self, parent, pres, config=None, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        if config is None:
            self.config = UniDecConfig()
        else:
            self.config = config
        # self.SetBackgroundColour("white")
        self.parent = parent
        self.pres = pres
        self.sizer = wx.GridBagSizer(wx.VERTICAL)
        self.SetSizer(self.sizer)

        # For p in plot_params, create a text box input
        for i, p in enumerate(plot_params):
            label = wx.StaticText(self, label=p)
            txt_ctrl = wx.TextCtrl(self, value=str(plot_params[p]), name=p)

            # Insert into sizer
            self.sizer.Add(label, pos=(i, 0), flag=wx.ALL, border=5)
            self.sizer.Add(txt_ctrl, pos=(i, 1), flag=wx.ALL, border=5)

        # Add in apply button
        self.apply_button = wx.Button(self, label="Apply")
        self.sizer.Add(self.apply_button, pos=(len(plot_params), 0), flag=wx.ALL, border=5)
        self.apply_button.Bind(wx.EVT_BUTTON, self.on_apply)

        # Finish
        self.sizer.Layout()
        self.Fit()

    def get_values_from_gui(self, event=None):
        for child in self.GetChildren():
            if isinstance(child, wx.TextCtrl):
                plot_params[child.GetName()] = child.GetValue()

    def set_values_to_gui(self, event=None):
        for child in self.GetChildren():
            if isinstance(child, wx.TextCtrl):
                child.SetValue(str(plot_params[child.GetName()]))

    def on_apply(self, event=None):
        self.get_values_from_gui()
        export_mplstyle(plot_params, filename=self.config.mplstylefile)
        self.pres.on_update_plot_params(stylefile=self.config.mplstylefile)


class PlotControlWindow(wx.Frame):
    '''
    Test window for PlotControlsPanel
    '''

    def __init__(self, parent, *args, **kwargs):
        wx.Frame.__init__(self, parent, *args, **kwargs)

        self.config = UniDecConfig()
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        plotpanel = wx.Panel(self)
        self.plot = Plot1d(plotpanel)
        sizer.Add(self.plot, 1, wx.EXPAND)

        self.panel = PlotControlsPanel(plotpanel, self)
        sizer.Add(self.panel, 0, wx.EXPAND)

        plotpanel.SetSizer(sizer)
        sizer.Fit(self)
        self.makeplot()
        self.Show()

    def makeplot(self):
        self.plot.plotrefreshtop([1, 2, 3], [1, 4, 9])
    def on_update_plot_params(self, event=None, stylefile=None):
        self.plot.update_style(stylefile)
        self.makeplot()



if __name__ == '__main__':
    app = wx.App(False)
    frame = PlotControlWindow(None, -1, "Plot Controls")
    app.MainLoop()
