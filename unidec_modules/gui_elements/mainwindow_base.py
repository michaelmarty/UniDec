import wx
import platform, os
from pubsub import pub
import numpy as np
import unidec_modules.miscwindows as miscwindows


class MainwindowBase(wx.Frame):
    """
    Main UniDec GUI Window.
    """

    def __init__(self, parent, title, config, iconfile="logo.ico", tabbed=None):
        """
        initialize window and feed in links to presenter and config.

        :param parent: GUniDec Presenter -> self.pres
        :param title: Window title (string)
        :param config: UniDecConfig object ->self.config
        :return: None
        """
        wx.Frame.__init__(self, None, title=title)
        # Set presenter and config
        self.pres = parent
        self.config = config
        self.title = title
        self.controls = None
        self.plots = []
        self.plotnames = []

        self.version = self.pres.eng.version

        if iconfile is None:
            iconfile = "logo.ico"
        # Set Icon File
        try:
            if os.path.isfile(iconfile):
                favicon = wx.Icon(iconfile, wx.BITMAP_TYPE_ANY)
                wx.Frame.SetIcon(self, favicon)
                self.icon_path = os.path.abspath(iconfile)
            else:
                self.icon_path = None
        except Exception as e:
            print(e)
            self.icon_path = None

        # Get display size and intelligently reshape
        self.system = platform.system()
        self.displaysize = wx.GetDisplaySize()
        pub.subscribe(self.on_motion, 'newxy')

    def launch(self):
        self.Centre()
        self.Show(True)
        pass

    def setup_shortcuts(self, keys):
        """
        Setup shortcuts in GUI. Binds key combinations to functions in presenter (self.pres)
        :return: None
        """

        ids = [k[2].GetId() for k in keys]
        tab = []
        for i, k in enumerate(keys):
            self.Bind(wx.EVT_MENU, k[1], id=ids[i])
            tab.append((wx.ACCEL_CTRL, ord(k[0]), ids[i]))
        self.SetAcceleratorTable(wx.AcceleratorTable(tab))
        pass

    def import_config_to_gui(self):
        self.controls.import_config_to_gui()

    def export_gui_to_config(self):
        self.controls.export_gui_to_config()

    def on_defaults(self, e):
        self.menu.on_defaults(e)

    def on_motion(self, xpos, ypos):
        """
        Triggered by pubsub from plot windows. Reports text in Status Bar.
        :param xpos: x position fed from event
        :param ypos: y position fed from event
        :return: None
        """
        try:
            if xpos is not None and ypos is not None:
                self.SetStatusText("x=%.2f y=%.2f" % (xpos, ypos), number=6)
        except:
            pass

    def on_open_dir(self, e):
        save_dir = self.config.udir
        print("Opening directory:", save_dir)
        try:
            os.system(self.config.opencommand + "\"" + save_dir + "\"")
        except Exception as err:
            print("Error opening directory", err)

    def clear_all_plots(self, flag=0):
        """
        Clear All Plots
        :return: None
        """
        for plot in self.plots:
            try:
                plot.clear_plot()
            except AttributeError:
                pass

    def clear_plots(self, flag=0):
        self.clear_all_plots(flag)
        try:
            self.peakpanel.clear_list()
        except:
            pass

    # .......................................................
    #
    #  Saving Figures
    #
    # .......................................................

    def save_all_figures(self, extension, extension2='', e=0, header=None, **kwargs):
        """
        Save All of the Figures. Will name as header+extension2+_FigureX.+exetension
        :param extension: Figure type (pdf, eps, png). Anything accepted by matplotlib
        :param extension2: Additional text to include in the figure header.
        :param e: Dummy wx Event
        :param header: Option to add different header. Default of none yields self.outfname as the path header
        :param kwargs: Any keywards to pass to the matplotlib savefig command such as Transparent or DPI
        :return: figureflags, files (the figures that were successfully saved and the files that they were saved to)
        """
        self.SetStatusText("Saving Figures", number=5)
        figureflags = []
        files = []
        if header is None:
            header = self.config.outfname + extension2
        else:
            header += extension2

        for i, plot in enumerate(self.plots):
            name1 = header + "_" + self.plotnames[i] + "." + extension
            if plot.flag:
                plot.on_save_fig(e, name1, **kwargs)
                figureflags.append(i + 1)
                files.append([i + 1, name1])
        return figureflags, files

    def on_save_figure_eps(self, e):
        """
        Save all figures as EPS
        :param e: Dummy wx event
        :return: None
        """
        self.SetStatusText("Saving", number=5)
        figureflags, files = self.save_all_figures("eps")
        self.SetStatusText("Saved to .eps", number=5)
        return figureflags, files

    def on_save_figure_png(self, e, **kwargs):
        """
        Save all figures as PNG
        :param e: Dummy wx event
        :param kwargs: keywards to pass to matplotlib savefig
        :return: None
        """
        self.SetStatusText("Saving", number=5)
        flags, self.pngs = self.save_all_figures("png", **kwargs)
        self.SetStatusText("Saved to .png", number=5)
        pass

    def on_save_figure_pdf(self, e):
        """
        Saves all figures as PDF
        :param e: Dummy wx event
        :return: None
        """
        self.SetStatusText("Saving", number=5)
        figureflags, files = self.save_all_figures("pdf")
        self.SetStatusText("Saved to .pdf", number=5)
        return figureflags, files

    def on_save_figure_dialog(self, e):
        """
        Open dialog box to set the parameters for figure type, size, and path to save.
        :param e: Dummy wx event
        :return: None
        """
        figsize = self.plots[0].GetSize()
        dpi = wx.ScreenDC().GetPPI()
        defaultsize = np.array([figsize[0] / dpi[0], figsize[1] / dpi[1]])
        defaultrect = np.array([0.1, 0.1, 0.8, 0.8])
        figureflags = []
        files = []

        dlg = miscwindows.SaveFigureDialog(self)
        dlg.initialize_interface(self.config)
        code = dlg.ShowModal()
        if code == 0:
            directory = dlg.directory
            header = dlg.header
            extension = dlg.extension
            transparent = dlg.transparent
            dpi = dlg.dpi
            try:
                dpi = int(dpi)
            except Exception as e:
                print(e, dpi)
                dpi = None

            path = os.path.join(directory, header)
            self.figsize = np.array(dlg.figsize)
            self.rect = np.array(dlg.rect)


            if not np.all(defaultsize == self.figsize) or not np.all(defaultrect == self.rect):
                plots = self.shrink_all_figures()
                self.SetStatusText("Saving Figures", number=5)
                figureflags, files = self.save_all_figures(extension, extension2="", header=path,
                                                           transparent=transparent, dpi=dpi)
                for plot in plots:
                    plot.resize = 1
                    plot.size_handler()
            else:
                figureflags, files = self.save_all_figures(extension, extension2="", header=path, transparent=transparent, dpi=dpi)
            # print self.directory
            self.SetStatusText("Saved Figures", number=5)
            # TODO: Remember Values from previous
        return figureflags, files

    def shrink_figure(self, plot, figsize=None, rect=None):
        """
        Automatically shrinks the plot to a figure size in inches set in self.figsize.
        :param plot: Plot object to shrink
        :return: None
        """
        if figsize is None:
            figsize = self.figsize
        if rect is None:
            try:
                rect = self.rect
            except:
                rect = [0.1, 0.1, 0.8, 0.8]
                self.rect = rect
        print(self.rect)
        if plot.flag:
            dpi = wx.ScreenDC().GetPPI()
            figsize2 = (int(figsize[0] * dpi[0]), int(figsize[1] * dpi[1]))
            print(figsize2, figsize)
            plot.resize = 0
            plot.canvas.SetSize(figsize2)
            plot.canvas.draw()
            # plot.set_nticks(5)
            plot.subplot1.set_position(rect)
            if plot.cbar is not None:
                plot.cbar.ax.set_position([0.85, 0.2, 0.05, 0.7])
            plot.repaint()

    def shrink_all_figures(self, figsize=None, rect=None):
        """
        Shrinks all figures to the size specified in self.figsize
        :return: A list of plot objects that we shrunk
        """
        if figsize is None:
            figsize = self.figsize
        for plot in self.plots:
            self.shrink_figure(plot, figsize=figsize, rect=rect)
        return self.plots

    def on_save_figure_small(self, e):
        """
        Preset to shrink figures to 4.5 in by 3 in and save as PDF.
        :param e: Dummy wx event
        :return: None
        """
        self.figsize = (4.5, 3.0)
        self.rect = [0.2, 0.2, 0.6, 0.7]
        plots = self.shrink_all_figures()
        self.SetStatusText("Saving Figures", number=5)
        figureflags, files = self.save_all_figures("pdf", extension2="_Thumb")
        self.SetStatusText("Saved to Thumbnails", number=5)
        for plot in plots:
            plot.resize = 1
            plot.size_handler()
        pass

    def onFocus(self, event):
        self.plotpanel.SetFocus()

    def on_about(self, e):
        """
        Displays message about program
        :param e:
        :return:
        """
        dlg = wx.MessageDialog(self,
                               "UniDec GUI version " + self.version +
                               "\nPlease contact mtmarty@email.arizona.edu with any questions, bugs, or features to add.\n"
                               "The latest version may be found at https://github.com/michaelmarty/UniDec/releases.\n"
                               "RawFileReader reading tool. Copyright © 2016 by Thermo Fisher Scientific, Inc. All rights reserved.\n"
                               "If used in publication, please cite Marty et Al. Anal. Chem. 2015, DOI: 10.1021/acs.analchem.5b00140 ",
                               "About UniDec", wx.OK | wx.CENTER)
        dlg.ShowModal()
        dlg.Destroy()

    def on_exit(self, e=None):
        """
        Exit the Program
        :param e: Dummy wx event
        :return: None
        """
        pub.unsubAll()
        self.Close(True)


