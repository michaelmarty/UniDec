import numpy as np
from matplotlib.animation import FuncAnimation
import wx
from unidec.modules import PlottingWindow
from unidec.modules import unidecstructure
from unidec.modules import unidectools as ud
from unidec.modules.isolated_packages import FileDialogs
import os

__author__ = 'Michael.Marty'


class AnimationWindow(wx.Frame):
    def __init__(self, parent, data_list, config=None, yvals=None, mode="1D", pks=None, pksmode="mz", *args, **kwargs):
        """
        A simple window for animating mulitple 1D or 2D plots in a sequence.
        :param parent: Parent window. Passed to wx.Frame.
        :param data_list: List of data to be plotted
        :param config: UniDecConfig object
        :param yvals: Titles for the plots.
        :param mode: 1 = 1D plots, 2 = 2D plots
        :param args:
        :param kwargs: 
        :return: None
        """
        wx.Frame.__init__(self, parent, title="Plot Animations", size=(-1, -1))
        # Initialize parameters
        if mode == "2D":
            self.mode = 2
        else:
            self.mode = 1
        if config is None:
            self.config = unidecstructure.UniDecConfig()
            self.config.initialize()
        else:
            self.config = config

        self.datalist = data_list
        self.pks = pks
        self.pksmode = pksmode

        if self.pksmode == "mz":
            self.xlabel = "m/z (Th)"
            self.testkda = False
        elif self.pksmode == "CCS":
            self.xlabel = "CCS"
            self.testkda = False
        else:
            self.xlabel = "Mass (Da)"
            self.testkda = True

        self.yvals = yvals
        if self.yvals is None:
            self.yvals = list(range(0, len(data_list)))

        self.dim = 1
        self.pos = -1
        self.play = False

        self.animation = None

        # Create GUI

        # Make the menu
        filemenu = wx.Menu()
        menu_save_fig = filemenu.Append(wx.ID_ANY, "Save Figures",
                                        "Save all figures at selected path")
        self.Bind(wx.EVT_MENU, self.on_save_fig, menu_save_fig)

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        self.SetMenuBar(menu_bar)

        self.CreateStatusBar(2)
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)

        if self.mode == 1:
            self.plot = PlottingWindow.Plot1d(panel)
        else:
            self.plot = PlottingWindow.Plot2d(panel)
        sizer.Add(self.plot, 0, wx.EXPAND)

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)

        sb = wx.StaticBox(panel, label='Frame Rate (ms/frame)')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
        frmax = 2000

        frmin = 1
        self.frslider = wx.Slider(panel, wx.ID_ANY, 500, frmin, frmax, (30, 60), (250, -1),
                                  wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
        self.frslider.SetTickFreq(100)
        sbs.Add(self.frslider, 0, wx.EXPAND)
        self.Bind(wx.EVT_COMMAND_SCROLL_THUMBRELEASE, self.update_framerate, self.frslider)
        controlsizer.Add(sbs, 0, wx.EXPAND)

        self.playbutton = wx.ToggleButton(panel, label="Play")
        self.nextbutton = wx.Button(panel, label="Next")
        self.backbutton = wx.Button(panel, label="Back")

        controlsizer.Add(self.backbutton, 0, wx.EXPAND)
        controlsizer.Add(self.playbutton, 0, wx.EXPAND)
        controlsizer.Add(self.nextbutton, 0, wx.EXPAND)

        self.Bind(wx.EVT_TOGGLEBUTTON, self.on_play, self.playbutton)
        self.Bind(wx.EVT_BUTTON, self.on_next, self.nextbutton)
        self.Bind(wx.EVT_BUTTON, self.on_back, self.backbutton)

        self.ctlautoscale = wx.CheckBox(panel, label="Autoscale")
        controlsizer.Add(self.ctlautoscale, 0, wx.EXPAND)
        if self.mode == 2:
            self.ctlautoscale.SetValue(True)

        sizer.Add(controlsizer, 0, wx.EXPAND)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close, self)

        self.init()
        self.Centre()
        self.Show(True)

    def on_close(self, e):
        """
        Stop the animation and close the window.
        :param e: Unused event
        :return: None
        """
        print("Closing")
        try:
            self.animation._stop()
        except:
            pass
        self.Destroy()

    def update(self, frame_number):
        """
        Continues to increment the plot to the next value in the data_list. Will stop if self.play is False.
        :param frame_number: Required but unused. Filled by FuncAnimation.
        :return: 0
        """
        if self.play:
            self.pos += 1
            self.update_plot()
            return 0
        else:
            return 0

    def update_plot(self):
        """
        Increment to the next data set and update the plot with the new data.
        Tries to keep some of the old plotting parameters like the zoom the same.
        Stops the animation if an error occurs.
        :return: None
        """
        try:
            if self.pks is not None:
                self.refresh_plot()
                self.add_peaks()
                return
            self.pos %= len(self.datalist)
            newdata = self.datalist[self.pos]
            title = str(self.yvals[self.pos])

            if self.mode == 1:
                # 1D Plot
                line = self.plot.subplot1.lines[0]
                xlim = self.plot.subplot1.get_xlim()
                ylim = self.plot.subplot1.get_ylim()
                line.set_data(newdata[:, 0], newdata[:, 1])

                self.plot.subplot1.set_title(title)

                autoflag = self.ctlautoscale.GetValue()
                if autoflag:
                    self.plot.subplot1.set_autoscale_on(True)
                    self.plot.subplot1.relim()
                    self.plot.subplot1.autoscale_view(True, True, True)
                else:
                    self.plot.subplot1.set_xlim(xlim)
                    self.plot.subplot1.set_ylim(ylim)
                self.plot.repaint()

            else:
                # 2D plot
                xlim = self.plot.subplot1.get_xlim()
                ylim = self.plot.subplot1.get_ylim()
                self.plot.contourplot(newdata, self.config, xlab=self.xlabel, title=title, repaint=False)
                autoflag = self.ctlautoscale.GetValue()
                if not autoflag:
                    self.plot.subplot1.set_xlim(xlim)
                    self.plot.subplot1.set_ylim(ylim)
                self.plot.add_title(title)
                self.plot.repaint()
        except Exception as e:
            self.animation._stop()
            print(e)

    def init(self):
        """
        Create a fresh plot and start the animation.
        :return: None
        """
        self.pos = 0
        self.refresh_plot()
        self.animation = FuncAnimation(self.plot.figure, self.update, interval=500)
        self.animation._start()

    def on_play(self, e):
        """
        Toggles self.play on or off. Will break the while loop in self.update if self.play = False.
        :param e: Unused event
        :return: None
        """
        tog = self.playbutton.GetValue()
        if tog:
            self.play = True
        else:
            self.play = False
        pass

    def refresh_plot(self):
        """
        Create a fresh plot from the top.
        :return: None
        """
        self.pos %= len(self.datalist)
        newdata = self.datalist[self.pos]
        title = str(self.yvals[self.pos])

        if self.mode == 1:
            self.plot.plotrefreshtop(newdata[:, 0], newdata[:, 1], title, self.xlabel, "Intensity", "", self.config,
                                     test_kda=self.testkda)
            self.plot.add_title(title)
            if self.pks is not None:
                self.add_peaks()
        else:
            self.plot.contourplot(newdata, self.config, xlab=self.xlabel, title=title)
            self.plot.add_title(title)

    def on_next(self, e):
        """
        Plot the next data set in data_list.
        :param e: Unused event
        :return: None
        """
        self.pos += 1
        self.update_plot()
        pass

    def on_back(self, e):
        """
        Plot the previous data set in data_list.
        :param e: Unused event
        :return: None
        """
        self.pos -= 1
        self.update_plot()
        pass

    def update_framerate(self, e):
        """
        Change the frame rate. Restart the animation with fresh frame rate.
        :param e: Unused event
        :return: None
        """
        framerate = self.frslider.GetValue()
        # print "Updated framerate to:", framerate
        # self.animation._interval=framerate
        # self.animation.new_frame_seq()
        self.animation._stop()
        self.animation = FuncAnimation(self.plot.figure, self.update, interval=framerate)
        self.animation._start()

    def add_peaks(self, e=None):
        if self.pksmode == "mz":
            self.add_peaks_mz()
        else:
            self.add_peaks_mass()

    def add_peaks_mz(self, e=None):
        for p in self.pks.peaks:
            if p.ignore == 0:
                list1 = []
                list2 = []
                mztab = p.mztab[self.pos]
                mztab2 = p.mztab2[self.pos]
                if (not ud.isempty(mztab)) and (not ud.isempty(mztab2)):
                    mztab = np.array(mztab)
                    mztab2 = np.array(mztab2)
                    maxval = np.amax(p.mztab[:, :, 1])
                    for k in range(0, len(mztab)):
                        if mztab[k, 1] > self.config.peakplotthresh * maxval:
                            list1.append(mztab2[k, 0])
                            list2.append(mztab2[k, 1])
                            # print mztab[k]
                    self.plot.plotadddot(np.array(list1), np.array(list2), p.color, p.marker)
        self.plot.repaint()

    def add_peaks_mass(self, e=None):
        for p in self.pks.peaks:
            if p.ignore == 0:
                pos = ud.nearest(self.datalist[self.pos][:, 0], p.mass)
                data = self.datalist[self.pos][pos]
                if data[1] > self.config.peakplotthresh * np.amax(self.datalist[self.pos][:, 1]):
                    self.plot.plotadddot(data[0], data[1], p.color, p.marker)
        self.plot.repaint()

    def on_save_fig(self, e=None):
        path = FileDialogs.save_file_dialog()
        base, ext = os.path.splitext(path)
        #self.init()
        #self.save_fig(base,ext)
        for i in range(0,len(self.datalist)):
            self.on_next(None)
            self.save_fig(base, ext)

    def save_fig(self, base, ext):
        path=base+str(self.pos)+ext
        print(self.pos, path)
        self.plot.save_figure(path)



# Main App Execution
if __name__ == "__main__":
    x = np.arange(0.0, 10.0)
    datalist = np.array([np.transpose([x, x]), np.transpose([x, x * x]), np.transpose([x, x * x * x])])
    app = wx.App(False)
    frame = AnimationWindow(None, datalist)
    app.MainLoop()
