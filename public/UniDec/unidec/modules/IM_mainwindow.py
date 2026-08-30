"""Dedicated main-window and plot layout for the ion mobility application."""

import wx

from unidec.modules.mainwindow import Mainwindow
from unidec.modules import miscwindows
from unidec.modules.plotting import ColorPlot, PlottingWindow, plot3d


class IMMainwindow(Mainwindow):
    """UniDec view configured permanently for IM-MS controls and plots."""

    mode_imflag = 1

    def __init__(self, parent, title, config, iconfile=None, tabbed=None):
        super().__init__(parent, title, config, iconfile=iconfile, tabbed=tabbed)

    def setup_mode_tabbed_plots(self, plotwindow, figsize):
        panels = {
            "data": wx.Panel(plotwindow),
            "fit": wx.Panel(plotwindow),
            "ccs": wx.Panel(plotwindow),
            "charges": wx.Panel(plotwindow),
            "mass_ccs": wx.Panel(plotwindow),
            "ccs_z": wx.Panel(plotwindow),
            "mz_cube": wx.Panel(plotwindow),
            "mass_cube": wx.Panel(plotwindow),
        }
        self.mode_tab_panels = panels
        self.plot1im = PlottingWindow.Plot2d(panels["data"], figsize=figsize)
        self.plot1fit = PlottingWindow.Plot2d(panels["fit"], figsize=figsize)
        self.plot2ccs = PlottingWindow.Plot1d(panels["ccs"], figsize=figsize)
        self.plot3color = ColorPlot.ColorPlot2D(panels["charges"], figsize=figsize)
        self.plot5mccs = PlottingWindow.Plot2d(panels["mass_ccs"], figsize=figsize)
        self.plot5ccsz = PlottingWindow.Plot2d(panels["ccs_z"], figsize=figsize)
        self.plot9 = plot3d.CubePlot(panels["mz_cube"], figsize=figsize)
        self.plot10 = plot3d.CubePlot(panels["mass_cube"], figsize=figsize)
        for panel, plot in (
                (panels["data"], self.plot1im),
                (panels["fit"], self.plot1fit),
                (panels["ccs"], self.plot2ccs),
                (panels["charges"], self.plot3color),
                (panels["mass_ccs"], self.plot5mccs),
                (panels["ccs_z"], self.plot5ccsz),
                (panels["mz_cube"], self.plot9),
                (panels["mass_cube"], self.plot10),
        ):
            miscwindows.setup_tab_box(panel, plot)

    def add_mode_tabs(self, plotwindow, section):
        panels = self.mode_tab_panels
        if section == "mz":
            plotwindow.AddPage(panels["data"], "IM-MS Data")
            plotwindow.AddPage(panels["fit"], "IM-MS Fit")
            plotwindow.AddPage(panels["charges"], "IM-MS Charges")
            plotwindow.AddPage(panels["mz_cube"], "m/z Cube")
        elif section == "ccs":
            plotwindow.AddPage(panels["ccs"], "CCS Distribution")
        elif section == "mass":
            plotwindow.AddPage(panels["mass_ccs"], "Mass vs. CCS ")
            plotwindow.AddPage(panels["ccs_z"], "CCS vs. Charge")
            plotwindow.AddPage(panels["mass_cube"], "Mass Cube")

    def setup_mode_scrolled_plots(self, plotwindow, figsize):
        self.plot1im = PlottingWindow.Plot2d(plotwindow, figsize=figsize)
        self.plot1fit = PlottingWindow.Plot2d(plotwindow, figsize=figsize)
        self.plot2ccs = PlottingWindow.Plot1d(plotwindow, figsize=figsize)
        self.plot5mccs = PlottingWindow.Plot2d(plotwindow, figsize=figsize)
        self.plot5ccsz = PlottingWindow.Plot2d(plotwindow, figsize=figsize)
        self.plot3color = ColorPlot.ColorPlot2D(plotwindow, figsize=figsize)
        self.plot9 = plot3d.CubePlot(plotwindow, figsize=figsize)
        self.plot10 = plot3d.CubePlot(plotwindow, figsize=figsize)

    def layout_scrolled_plots(self):
        placements = (
            (self.plot1, (0, 0)),
            (self.plot1im, (0, 1)),
            (self.plot3color, (1, 0)),
            (self.plot1fit, (1, 1)),
            (self.plot2, (2, 0)),
            (self.plot3, (2, 1)),
            (self.plot2ccs, (3, 0)),
            (self.plot5, (3, 1)),
            (self.plot4, (4, 0)),
            (self.plot5mccs, (4, 1)),
            (self.plot6, (5, 0)),
            (self.plot5ccsz, (5, 1)),
            (self.plot9, (6, 0)),
            (self.plot10, (6, 1)),
        )
        for plot, position in placements:
            self.sizerplot.Add(plot, position, span=(1, 1), flag=wx.EXPAND)

    def extend_mode_plot_metadata(self):
        self.plots.extend(
            [
                self.plot1im,
                self.plot1fit,
                self.plot2ccs,
                self.plot5mccs,
                self.plot5ccsz,
                self.plot3color,
                self.plot9,
                self.plot10,
            ]
        )
        self.plotnames.extend(
            [
                "Figure1im",
                "Figure1fit",
                "Figure2ccs",
                "Figure5massccs",
                "Figure5ccsz",
                "Figure3color",
                "mzCube",
                "massCube",
            ]
        )
