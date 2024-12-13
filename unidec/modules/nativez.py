import sys
import time
import wx
import numpy as np
# import matplotlib.cm as cm
import matplotlib as mpl
from scipy.interpolate import interp1d
from wx.lib.agw import ultimatelistctrl as ulc
from unidec import tools as ud
from modules.plotting import PlottingWindow
from unidec.modules import MassFitter

__author__ = 'Michael.Marty'


class Zoffset:
    def __init__(self):
        """
        A simple class for defining a charge offset species.
        :return: None
        """
        self.offset = 0
        self.intensity = 0
        self.index = 0
        self.color = [255, 255, 0, 255]
        self.marker = "."
        self.id = 0
        self.width = 0
        self.nstate = 0
        self.extractwidth = 1
        self.extract = []

    def make(self, offset, intensity, index, color, marker):
        """
        Add details to the charge offset species.
        :param offset: Charge offset value
        :param intensity: Intensity
        :param index: Index in the list
        :param color: Color
        :param marker: Marker for the plot
        :return: None
        """
        self.offset = offset
        self.intensity = intensity
        self.index = index
        self.color = color
        self.marker = marker


class NativeZ(wx.Dialog):
    def __init__(self, *args, **kwargs):
        """
        Create a dialog for examining the native charge/mass relationship and extraction specific CID states.
        Tries to create a big window, but will rescale to smaller if the screen is too small.
        :param args: Passed to wx.Dialog
        :param kwargs: Passed to wx.Dialog
        :return: None
        """
        wx.Dialog.__init__(self, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER, *args, **kwargs)
        defaultsize = [1400, 1000]
        displaysize = wx.GetDisplaySize()
        self.figsize = (4.75, 3)
        if defaultsize[0] > displaysize[0]:
            defaultsize[0] = np.round(displaysize[0] * 0.9)
            self.figsize = (3, 2)
        if defaultsize[1] > displaysize[1]:
            defaultsize[1] = np.round(displaysize[1] * 0.9)
            self.figsize = (3, 2)
        print(defaultsize, displaysize)
        self.SetSize(defaultsize)

        self.SetTitle("Native Charge Tools")

        self.config = None
        self.pks = None
        self.zoffs = []
        self.massoffset = 0
        self.massaxis = None
        self.chargeaxis = None
        self.igrid = None
        self.offset_totals = None
        self.offset_grid = None

        self.plot1 = None
        self.plot2 = None
        self.plot3 = None
        self.plot4 = None
        self.plot5 = None
        self.plot6 = None
        self.plot7 = None
        self.zlistctrl = None
        self.ctlmassoffset = None
        self.ctlfilt = None

    def initialize_interface(self, massaxis, chargeaxis, igrid, config, pks):
        """
        Initialize the parameters, setup the GUI, and plot the intial results.
        :param massaxis: Mass axis values
        :param chargeaxis: Charge axis value
        :param igrid: Intensities at each mass and charge point
        :param config: UniDecConfig object
        :param pks: Peaks object
        :return: None
        """
        # Initialize the parameters
        self.config = config
        self.pks = pks

        self.massaxis = np.array(massaxis)
        self.chargeaxis = np.array(chargeaxis)
        self.igrid = np.reshape(igrid, (len(massaxis), len(chargeaxis)))

        # Setup the GUI
        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        sb = wx.StaticBox(pnl, label='Set Parameters to Plot Native Z')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        self.plot1 = PlottingWindow.Plot1d(pnl, figsize=self.figsize)
        self.plot2 = PlottingWindow.Plot1d(pnl, figsize=self.figsize)
        self.plot3 = PlottingWindow.Plot1d(pnl, figsize=self.figsize)
        self.plot4 = PlottingWindow.Plot1d(pnl, figsize=self.figsize)
        self.plot5 = PlottingWindow.Plot1d(pnl, figsize=self.figsize)
        self.plot6 = PlottingWindow.Plot1d(pnl, figsize=self.figsize)
        self.plot7 = PlottingWindow.Plot2d(pnl, figsize=self.figsize)

        hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        hbox0.Add(self.plot1, 0, wx.EXPAND)
        hbox0.Add(self.plot3, 0, wx.EXPAND)
        hbox0.Add(self.plot4, 0, wx.EXPAND)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2.Add(self.plot2, 0, wx.EXPAND)
        hbox2.Add(self.plot5, 0, wx.EXPAND)
        hbox2.Add(self.plot6, 0, wx.EXPAND)

        sbs.Add(hbox0, 0, wx.EXPAND)
        sbs.Add(hbox2, 0, wx.EXPAND)

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        addbutton = wx.Button(pnl, label="Add Line")
        fitbutton = wx.Button(pnl, label="Fit")
        resetbutton = wx.Button(pnl, label="Reset to Default")
        extractbutton = wx.Button(pnl, label="Extract")
        self.ctlmassoffset = wx.TextCtrl(pnl, value=str(self.config.massoffset), size=(50, -1))
        self.ctlfilt = wx.RadioBox(pnl, label="Extract Shape", choices=["Box", "Gaussian"])
        self.ctlfilt.SetSelection(self.config.extractshape)
        savefigbutton = wx.Button(pnl, label="Save Figures")
        replotbutton = wx.Button(pnl, label="Replot")
        hbox1.Add(addbutton, 0)
        hbox1.Add(replotbutton, 0)
        hbox1.Add(resetbutton, 0)
        hbox1.Add(wx.StaticText(pnl, label="     "), 0, wx.ALIGN_CENTER_VERTICAL)

        hbox1.Add(fitbutton, 0)
        hbox1.Add(extractbutton, 0)
        hbox1.Add(wx.StaticText(pnl, label="     Monomer Mass: "), 0)  # , wx.ALIGN_CENTER_VERTICAL)
        hbox1.Add(self.ctlmassoffset, 0)
        hbox1.Add(self.ctlfilt, 0)

        hbox1.Add(savefigbutton, 0)
        sbs.Add(hbox1, 0, wx.EXPAND)

        hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.zlistctrl = ColorList(pnl)
        hbox3.Add(self.zlistctrl, 1, wx.EXPAND)
        hbox3.Add(self.plot7, 0, wx.EXPAND)
        sbs.Add(hbox3, 0, wx.EXPAND)

        pnl.SetSizer(sbs)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self, label='Ok')
        closebutton = wx.Button(self, label='Cancel')
        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, flag=wx.ALL | wx.EXPAND, border=5)
        vbox.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        self.Center()

        # Bind events
        self.zlistctrl.ultimateList.Bind(wx.EVT_BUTTON, self.on_delete)
        self.zlistctrl.ultimateList.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.update)

        fitbutton.Bind(wx.EVT_BUTTON, self.fit)
        addbutton.Bind(wx.EVT_BUTTON, self.on_add)
        resetbutton.Bind(wx.EVT_BUTTON, self.on_reset)
        extractbutton.Bind(wx.EVT_BUTTON, self.extract)
        savefigbutton.Bind(wx.EVT_BUTTON, self.save_figures)
        replotbutton.Bind(wx.EVT_BUTTON, self.on_replot)
        okbutton.Bind(wx.EVT_BUTTON, self.on_close)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)

        # Set Range Here and plot the initial results
        tstart = time.perf_counter()
        self.make_f_array(-50, 15)
        tend = time.perf_counter()
        print("F Array Time: %.2gs" % (tend - tstart))
        if self.config.zoffs == []:
            self.get_maxima()
        else:
            self.zoffs = self.config.zoffs
            self.plot_zoffs()
        self.populate_list(0)

    def update(self, e):
        """
        Update self.zoffs from the self.zlistctrl. Update the intensities for each zoffset value.
        :param e:
        :return:
        """
        self.zoffs = self.zlistctrl.return_data()
        zfunction = interp1d(self.offset_totals[:, 0], self.offset_totals[:, 1])
        for z in self.zoffs:
            if np.amin(self.offset_totals[:, 0]) <= z.offset <= np.amax(self.offset_totals[:, 0]):
                z.intensity = zfunction(z.offset)
            else:
                z.intensity = 0
            self.zlistctrl.ultimateList.SetStringItem(z.index, 1, str(z.intensity))

    def on_replot(self, e):
        """
        Update the parameters and replot the results.
        :param e: Unused event
        :return: None
        """
        self.update(e)
        self.plot_zoffs()
        self.update_list()

    def plot_zoffs(self):
        """
        Plot the total charge offsets along with fits. Call self.make_plot_7.
        :return: None
        """
        self.plot1.plotrefreshtop(self.offset_totals[:, 0], self.offset_totals[:, 1], title="Native Z Offset",
                                  xlabel="Offset",
                                  ylabel="Intensity", zoom="span")
        for i in range(0, len(self.zoffs)):
            self.plot1.plotadddot(self.zoffs[i].offset, self.zoffs[i].intensity,
                                  np.array(self.zoffs[i].color[:3]) / 255, self.zoffs[i].marker)
        self.plot1.repaint()
        self.make_plot_7(0)

    def make_plot_7(self, e):
        """
        Plot the 2D mass v. charge plot.
        Add colored charge offset bands on top of the plot.
        :param e: Unused event
        :return: None
        """
        self.plot7.contourplot(xvals=self.massaxis, yvals=self.chargeaxis, zgrid=self.igrid, config=self.config,
                               title="Mass vs. Charge", test_kda=True)
        try:
            for zoff in self.zoffs:
                eshape = self.ctlfilt.GetSelection()
                self.plot7.plot_native_z(zoff.offset, np.array(zoff.color) / 255, self.massaxis,
                                         width=zoff.extractwidth, alpha=0.5, shape=eshape)
        except Exception as e:
            print("Failed to plot", e)
            pass

    def on_delete(self, event):
        """
        Delete item from the list ctrl.
        :param event: wx.Button event
        :return: None
        """
        self.update(event)
        btn = event.GetId()
        ids = [int(p.id) for p in self.zoffs]
        i = [i for i, j in enumerate(ids) if j == btn][0]
        index = self.zoffs[i].index
        self.zlistctrl.ultimateList.DeleteItem(index)
        self.update(0)

    def fit(self, e):
        """
        Fit the offset totals to a series of overlapping Gaussians.
        Update the plots and listctrl from the result.
        :param e: Unused event
        :return: None
        """
        self.update(e)
        guess = []
        for z in self.zoffs:
            guess.append([z.offset, z.intensity * 1 / (0.5 * np.sqrt(2 * np.pi))])
        guess = np.array(guess)
        mf = MassFitter.MassFitter(self.offset_totals, guess, 0, "smallguess")
        fitspec, fitres = mf.perform_fit("nonorm", "smallguess")
        print("Output", fitres)

        for i in range(0, len(self.zoffs)):
            fit = fitres[i]
            self.zoffs[i].offset = fit[0]
            self.zoffs[i].width = fit[1]
            self.zoffs[i].intensity = fit[2] * 1 / (fit[1] * np.sqrt(2 * np.pi))
        self.update_list()
        self.plot_zoffs()
        for z in self.zoffs:
            yvals = ud.ndis_std(z.offset, self.offset_totals[:, 0], z.width, a=z.intensity, norm_area=False)
            self.plot1.plotadd(self.offset_totals[:, 0], yvals, np.array(z.color) / 255, "Fit")
        self.plot1.repaint()

    def on_reset(self, e):
        """
        Reset the listctrl and plots to the default states.
        :param e: Unused event
        :return: None
        """
        self.get_maxima()
        self.plot_zoffs()
        self.zlistctrl.super_delete()
        self.populate_list(e)

    def update_list(self):
        """
        Update self.zlistctrl from self.zoffs.
        :return: None
        """
        for z in self.zoffs:
            index = z.index
            self.zlistctrl.ultimateList.SetStringItem(index, 0, str(z.offset))
            self.zlistctrl.ultimateList.SetStringItem(index, 1, str(z.intensity))
            self.zlistctrl.ultimateList.SetStringItem(index, 2, str(z.width))

    def populate_list(self, e):
        """
        Add data from self.zoffs to self.zlistctrl.
        :param e: Unused event
        :return: None
        """
        for i in self.zoffs:
            self.zlistctrl.add_line(i)

    def on_add(self, e):
        """
        Add a blank line to self.zlistctrl.
        :param e: Unused event
        :return: None
        """
        self.zlistctrl.add_empty_line()

    def make_f_array(self, minval, maxval, zwidth=1):
        """
        Get the global offset parameters for the massaxis and chargeaxis.
        Create a grid of mass and charge value.
        Calculate the charge offset at each point and store it at self.offset_grid.
        Calculate the total intensity for each charge offset value +/- zwidth (default is 1).
        :param minval: Minimum charge offset value
        :param maxval: Maximum charge offset value
        :param zwidth: Tolerance window for the offset to be summed.
        :return: None
        """
        frange = np.arange(minval, maxval, 0.5)
        mgrid, zgrid = np.meshgrid(self.massaxis, self.chargeaxis, indexing='ij')
        self.offset_grid = ud.get_z_offset(mgrid, zgrid)
        ftot = []
        for f in frange:
            bool1 = self.offset_grid >= f - zwidth
            bool2 = self.offset_grid < f + zwidth
            bool3 = np.all([bool1, bool2], axis=0)
            ftot.append([f, np.sum(self.igrid[bool3])])
        self.offset_totals = np.array(ftot)

    def fast_extract(self, f, width, eshape):
        """
        Extract the intensity data for a specific charge offset from self.igrid. Mask all of the other data.
        :param f: Charge offset
        :param width: Width of the extraction (+/- from charge offset)
        :param eshape: Shape of the extraction
        0 = Box (+/- is a hard cutoff)
        1 = Gaussian (+/- is std deviation of Gaussian)
        :return: masked intensity values from self.igrid for a specific charge offset
        """
        if eshape == 1:
            width *= 3.
        bool1 = self.offset_grid >= f - width
        bool2 = self.offset_grid < f + width
        bool3 = np.all([bool1, bool2], axis=0)
        out = self.igrid * bool3
        if eshape == 1:
            wgrid = np.exp(-((self.offset_grid - f) ** 2.) / (2. * width * width))
            out = out * wgrid
        return out

    def get_maxima(self):
        """
        Detect peaks in self.offset_totals (the total extracted intensitites) and set default parameters.
        Plot the results.
        :return: None
        """
        defaultcolors = [[255, 0, 0, 255], [0, 0, 255, 255], [0, 255, 0, 255], [255, 0, 255, 255]]
        defaultmarkers = ['o', 'v', '^', '>', 's', 'd', '*']

        peaks = ud.peakdetect(self.offset_totals, window=2, threshold=0.1)
        print(peaks)
        self.zoffs = []
        for p in peaks:
            i = ud.nearest(self.offset_totals[:, 0], p[0])
            zoff = Zoffset()
            zoff.make(p[0], p[1], i, defaultcolors[len(self.zoffs) % len(defaultcolors)],
                      defaultmarkers[len(self.zoffs)])
            self.zoffs.append(zoff)
        self.plot_zoffs()

    def extract(self, e):
        """
        Extract the mass distribution at each offset value, plot it, and write the outputs.
        Call self.peak_extract
        :param e: Unused event
        :return: None
        """
        # Update the parameters and plots
        self.on_replot(e)
        self.plot2.clear_plot("nopaint")
        eshape = self.ctlfilt.GetSelection()
        self.massoffset = float(self.ctlmassoffset.GetValue())

        tstart = time.perf_counter()
        for i, z in enumerate(self.zoffs):
            # Extract mass v. intensity values for the offsets in self.zoffs
            intensity = self.fast_extract(z.offset, z.extractwidth, eshape)
            z.extract = np.transpose([self.massaxis + z.nstate * self.massoffset, np.sum(intensity, axis=1)])
            # Update plot2
            if not self.plot2.flag:
                self.plot2.plotrefreshtop(z.extract[:, 0], z.extract[:, 1], title="Extracted Intensities",
                                          xlabel="Mass (Da)", ylabel="Intensity", color=np.array(z.color) / 255,
                                          nticks=5, test_kda=True)
            else:
                self.plot2.plotadd(z.extract[:, 0], z.extract[:, 1], np.array(z.color) / 255, "Offset=" + str(z.offset))
            # Write each extract out
            fileout = self.config.outfname + "_extracted_intensities" + str(i) + ".txt"
            np.savetxt(fileout, z.extract)
        tend = time.perf_counter()
        print("Extraction Time: %.2gs" % (tend - tstart))
        self.plot2.repaint()

        # Extract Peaks
        self.peak_extract()

    def peak_extract(self):
        """
        Extract the values for local max (height) and area for peaks in pks from each zoff.extract.
        Plot the results.
        Write the results to files.
        :return: None
        """
        # TODO: Add extra peaks here to compensate for shifts in the peaks.
        peakextracts = np.zeros((len(self.zoffs), self.pks.plen))
        peakextractsarea = np.zeros((len(self.zoffs), self.pks.plen))

        try:
            xvals = self.pks.masses
        except AttributeError:
            print("No Peaks to Extract")
            return None

        if xvals is not None:
            # Extraction
            for i, z in enumerate(self.zoffs):
                for j, p in enumerate(self.pks.peaks):
                    x = z.extract[:, 0]
                    y = z.extract[:, 1]
                    index = ud.nearest(x, p.mass)
                    peakextracts[i, j] = ud.localmax2(y, index, int(self.config.peakwindow / self.config.massbins))
                    if not ud.isempty(p.integralrange):
                        integral, intdat = ud.integrate(z.extract, p.integralrange[0], p.integralrange[1])
                        peakextractsarea[i, j] = integral

            # Switch to subunit numbers
            if self.massoffset > 0:
                xvals = np.round(np.array(xvals) / self.massoffset)
                label = "Subunit Number"
            else:
                label = "Mass (Da)"

            # Plots
            # Create height extraction line plot
            sum1 = np.sum(peakextracts, axis=0)
            smax1 = np.amax(sum1)
            if np.amax(sum1) != 0:
                sum1 /= smax1
                self.plot3.plotrefreshtop(xvals, sum1, title="Extracted Heights", xlabel=label, ylabel="Intensity",
                                          integerticks=True, test_kda=True)
                for i, z in enumerate(self.zoffs):
                    self.plot3.plotadd(xvals, peakextracts[i] / smax1, np.array(z.color) / 255., str(i))
                self.plot3.repaint()

                # Create height extraction bar chart
                xvals2 = np.arange(0, len(sum1))
                # colormap = cm.get_cmap(self.config.peakcmap, len(xvals))
                colormap = mpl.colormaps[self.config.peakcmap].resampled(len(xvals))
                peakcolors = colormap(np.arange(len(xvals)))
                self.plot5.barplottop(xvals2, sum1, [int(i) for i in xvals], peakcolors, label,
                                      "Normalized Intensity", "Extracted Total Peak Heights")

            # Make Plots for Integrals
            sum2 = np.sum(peakextractsarea, axis=0)
            smax2 = np.amax(sum2)
            if np.amax(sum2) != 0:
                sum2 /= smax2
                # Make line plots
                self.plot4.plotrefreshtop(xvals, sum2, title="Extracted Areas", xlabel=label, ylabel="Area",
                                          integerticks=True, test_kda=True)
                for i, z in enumerate(self.zoffs):
                    self.plot4.plotadd(xvals, peakextractsarea[i] / smax2, np.array(z.color) / 255., str(i))
                self.plot4.repaint()
                # Make bar charts
                self.plot6.barplottop(xvals2, sum2, [int(i) for i in xvals], peakcolors, label,
                                      "Normalized Intensity", "Extracted Total Peak Areas")
            else:
                print("No integration provided")

            # Save total outputs for extracted peaks
            try:
                np.savetxt(self.config.outfname + "_extracted_heights.txt", peakextracts)
                np.savetxt(self.config.outfname + "_total_extracted_heights.txt",
                           np.transpose([self.pks.masses, xvals, sum1]))
                if np.amax(sum2) != 0:
                    np.savetxt(self.config.outfname + "_extracted_areas.txt", peakextractsarea)
                    np.savetxt(self.config.outfname + "_total_extracted_areas.txt",
                               np.transpose([self.pks.masses, xvals, sum2]))
            except Exception as e:
                print("Error saving files", e)

    def save_figures(self, e):
        """
        Save the figures as pdfs in the unidec folder.
        :param e: Unused event
        :return: None
        """
        extraheader = "_NativeZ"
        name1 = self.config.outfname + extraheader + "_Figure1.pdf"
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
        name2 = self.config.outfname + extraheader + "_Figure2.pdf"
        if self.plot2.flag:
            self.plot2.on_save_fig(e, name2)
        name3 = self.config.outfname + extraheader + "_Figure3.pdf"
        if self.plot3.flag:
            self.plot3.on_save_fig(e, name3)
        name4 = self.config.outfname + extraheader + "_Figure4.pdf"
        if self.plot4.flag:
            self.plot4.on_save_fig(e, name4)
        name5 = self.config.outfname + extraheader + "_Figure5.pdf"
        if self.plot5.flag:
            self.plot5.on_save_fig(e, name5)
        name6 = self.config.outfname + extraheader + "_Figure6.pdf"
        if self.plot6.flag:
            self.plot6.on_save_fig(e, name6)
        name7 = self.config.outfname + extraheader + "_Figure7.pdf"
        if self.plot7.flag:
            self.plot7.on_save_fig(e, name7)
        print("Saved to:", self.config.outfname + extraheader)

    def on_close(self, e):
        """
        Close the window and update the parameters
        :param e:
        :return:
        """
        self.update(e)
        self.config.zoffs = self.zoffs
        self.config.massoffset = self.massoffset
        self.config.extractshape = self.ctlfilt.GetSelection()
        self.Destroy()
        self.EndModal(0)

    def on_close_cancel(self, e):
        """
        Close the window but do not update the parameters.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        self.EndModal(1)


class ColorList(wx.Panel):
    def __init__(self, parent):
        """
        Create an ulitimape list control panel.
        :param parent: Parent passed to wx.Panel
        :return: None
        """
        wx.Panel.__init__(self, parent)
        self.ultimateList = ulc.UltimateListCtrl(self, size=(500, 200),
                                                 agwStyle=wx.LC_REPORT | wx.LC_VRULES | wx.LC_EDIT_LABELS | wx.LC_HRULES
                                                          | ulc.ULC_USER_ROW_HEIGHT | ulc.ULC_HAS_VARIABLE_ROW_HEIGHT
                                                          | ulc.ULC_EDIT_LABELS)

        info = ulc.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._image = []
        info._format = 0
        info._kind = 1
        info._text = "Native Z Offset"
        self.ultimateList.InsertColumnInfo(0, info)

        info = ulc.UltimateListItem()
        info._format = wx.LIST_FORMAT_RIGHT
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT | ulc.ULC_MASK_CHECK
        info._image = []
        info._text = "Intensity"
        self.ultimateList.InsertColumnInfo(1, info)

        info = ulc.UltimateListItem()
        info._format = wx.LIST_FORMAT_RIGHT
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT | ulc.ULC_MASK_CHECK
        info._image = []
        info._text = "Width"
        self.ultimateList.InsertColumnInfo(2, info)

        info = ulc.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = "Color"
        info._image = []
        self.ultimateList.InsertColumnInfo(3, info)

        info = ulc.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = "ID"
        info._image = []
        self.ultimateList.InsertColumnInfo(4, info)

        info = ulc.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = ""
        info._image = []
        self.ultimateList.InsertColumnInfo(5, info)

        info = ulc.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = "N"
        info._image = []
        self.ultimateList.InsertColumnInfo(6, info)

        info = ulc.UltimateListItem()
        info._mask = wx.LIST_MASK_TEXT | wx.LIST_MASK_IMAGE | wx.LIST_MASK_FORMAT
        info._format = 0
        info._text = "Z Width"
        info._image = []
        self.ultimateList.InsertColumnInfo(7, info)

        self.ultimateList.SetColumnWidth(0, 100)
        self.ultimateList.SetColumnWidth(1, 100)
        self.ultimateList.SetColumnWidth(2, 100)
        self.ultimateList.SetColumnWidth(3, 100)
        self.ultimateList.SetColumnWidth(4, 20)
        self.ultimateList.SetColumnWidth(5, 100)
        self.ultimateList.SetColumnWidth(6, 50)
        self.ultimateList.SetColumnWidth(7, 50)
        self.ultimateList.SetUserLineHeight(25)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.ultimateList, 0, wx.EXPAND)
        self.SetSizer(sizer)
        self.buttontot = 0

    def add_empty_line(self):
        """
        Add an empty line by adding a default Zoffset object.
        :return: None
        """
        self.add_line(Zoffset())

    def add_line(self, zoff):
        """
        Adds a new line from a Zoffset object.
        :param zoff: Zoffset object.
        :return:
        """
        index = self.ultimateList.InsertStringItem(sys.maxsize, str(zoff.offset))
        self.ultimateList.SetStringItem(index, 1, str(zoff.intensity))
        self.ultimateList.SetStringItem(index, 2, str(zoff.width))
        size = (50, 10)
        colorarray = list(zoff.color)
        if len(colorarray) < 4:
            colorarray.append(255)
        colorbox = wx.ColourPickerCtrl(self.ultimateList, size=size,
                                       colour=wx.Colour(int(colorarray[0]), int(colorarray[1]), int(colorarray[2]), alpha=int(colorarray[3])))
        self.ultimateList.SetItemWindow(index, col=3, wnd=colorbox, expand=True)
        self.ultimateList.SetStringItem(index, 4, str(self.buttontot))
        deletebutton = wx.Button(self.ultimateList, label="Delete", id=self.buttontot)
        self.ultimateList.SetItemWindow(index, col=5, wnd=deletebutton, expand=True)
        self.buttontot += 1
        textinput = wx.TextCtrl(self.ultimateList, value=str(zoff.nstate))
        self.ultimateList.SetItemWindow(index, col=6, wnd=textinput, expand=True)
        textinput2 = wx.TextCtrl(self.ultimateList, value=str(zoff.extractwidth))
        self.ultimateList.SetItemWindow(index, col=7, wnd=textinput2, expand=True)

    def return_data(self):
        """
        Make each line from the listctrl into a Zoffset object and return the list.
        :return: List of Zoffset objects
        """
        count = self.ultimateList.GetItemCount()
        zoffouts = []
        # print "Count",count
        defaultmarkers = ['o', 'v', '^', '>', 's', 'd', '*']
        for index in range(0, count):
            zout = Zoffset()
            zout.index = index
            # print index
            colorwindow = self.ultimateList.GetItemWindow(index, col=3)
            zout.color = colorwindow.GetColour()
            zout.offset = float(self.ultimateList.GetItem(index, col=0).GetText())
            zout.width = float(self.ultimateList.GetItem(index, col=2).GetText())
            # zout.intensity=float(self.ultimateList.GetItem(index,col=1).GetText())
            zout.id = int(self.ultimateList.GetItem(index, col=4).GetText())
            zout.marker = defaultmarkers[len(zoffouts)]
            inputbox = self.ultimateList.GetItemWindow(index, col=6)
            zout.nstate = int(inputbox.GetValue())
            inputbox = self.ultimateList.GetItemWindow(index, col=7)
            zout.extractwidth = int(inputbox.GetValue())
            zoffouts.append(zout)
        return zoffouts

    def super_delete(self):
        """
        Delete all items in the ultimate list control.
        :return: None
        """
        count = self.ultimateList.GetItemCount()
        topcount = count
        num = 0
        while count > 0 and num <= topcount:
            try:
                self.ultimateList.DeleteItemWindow(0, col=3)
                self.ultimateList.DeleteItem(0)
                count = self.ultimateList.GetItemCount()
                num += 1
            except Exception as e:
                num += 1
                print("Delete Failed:", e)
                pass
        self.ultimateList.DeleteAllItems()
