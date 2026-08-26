import wx
# noinspection PyUnresolvedReferences
import pandas as pd
# noinspection PyUnresolvedReferences
import numpy as np
# noinspection PyUnresolvedReferences
import os
from modules.plotting import PlottingWindow
import matplotlib.cm as cm
from pubsub import pub
import matplotlib.colors as colors
import unidec.tools as ud
from Scripts.LipiDec.SkyQuant import *

__author__ = 'Michael.Marty'

plotchoices = ["Automatic", "Class", "Tail Unsaturation", "Tail Length"]
normchoices = {"None": 0, "mTIC": 1, "STD": 2, "Both": 3}
luminance_cutoff = 135
white_text = wx.Colour(250, 250, 250)
black_text = wx.Colour(0, 0, 0)


class LipyDecWindow(wx.Frame):
    def __init__(self, parent, engine=None):
        wx.Frame.__init__(self, parent, title="LipyDec")  # ,size=(-1,-1))
        statusbar_log_silencer = wx.LogNull()
        self.CreateStatusBar(3)
        self.SetStatusWidths([-1, -3, -3])
        del statusbar_log_silencer
        pub.subscribe(self.on_motion, 'newxy')

        # Setup initial values
        try:
            self.directory = parent.directory
        except:
            self.directory = os.getcwd()

        self.df = None
        self.classes = []
        self.ngroups = 2
        self.nreps = 3
        self.normchoice = "None"

        self.outfname = "Figure"
        # Make the menu

        filemenu = wx.Menu()
        menu_open = filemenu.Append(wx.ID_ANY, "Open File", "Open File")
        self.Bind(wx.EVT_MENU, self.on_open, menu_open)

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        # menu_bar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(menu_bar)

        # Setup the GUI
        self.panel = wx.Panel(self)

        sizer = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.checkboxes = wx.CheckListBox(self.panel, size=(100, 400))
        self.Bind(wx.EVT_CHECKLISTBOX, self.update_checks, self.checkboxes)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.check_only, self.checkboxes)

        hbox.Add(self.checkboxes, 0, wx.EXPAND)

        self.plot1 = PlottingWindow.Plot1d(self.panel, figsize=(6, 4))

        plotsizer1 = wx.BoxSizer(wx.VERTICAL)
        plotsizer1.Add(self.plot1, 1, wx.EXPAND)
        hbox.Add(plotsizer1, 1, wx.EXPAND)

        self.plot2 = PlottingWindow.Plot1d(self.panel, figsize=(6, 4))

        plotsizer2 = wx.BoxSizer(wx.VERTICAL)
        plotsizer2.Add(self.plot2, 1, wx.EXPAND)
        hbox.Add(plotsizer2, 1, wx.EXPAND)

        self.plot4 = PlottingWindow.Plot1d(self.panel, figsize=(4, 4))

        plotsizer2 = wx.BoxSizer(wx.VERTICAL)
        plotsizer2.Add(self.plot4, 1, wx.EXPAND)
        hbox.Add(plotsizer2, 1, wx.EXPAND)

        sizer.Add(hbox, 1, wx.EXPAND)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)

        vbox = wx.BoxSizer(wx.VERTICAL)
        self.checkboxes2 = wx.CheckListBox(self.panel, size=(100, 200))
        self.Bind(wx.EVT_CHECKLISTBOX, self.update_checks, self.checkboxes2)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.check_only2, self.checkboxes2)
        vbox.Add(self.checkboxes2, 0, wx.EXPAND)

        self.checkboxes3 = wx.CheckListBox(self.panel, size=(100, 200))
        self.Bind(wx.EVT_CHECKLISTBOX, self.update_checks, self.checkboxes3)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.check_only3, self.checkboxes3)
        vbox.Add(self.checkboxes3, 0, wx.EXPAND)

        self.plotchoices = wx.ComboBox(self.panel, choices=plotchoices, style=wx.CB_READONLY)
        self.plotchoices.SetSelection(1)
        self.Bind(wx.EVT_COMBOBOX, self.update_checks, self.plotchoices)
        vbox.Add(self.plotchoices, 0, wx.EXPAND)

        self.normchoices = wx.ComboBox(self.panel, choices=list(normchoices.keys()), style=wx.CB_READONLY)
        self.normchoices.SetSelection(1)
        self.Bind(wx.EVT_COMBOBOX, self.update_checks, self.normchoices)
        vbox.Add(self.normchoices, 0, wx.EXPAND)

        hbox2.Add(vbox, 0, wx.EXPAND)

        self.plot3 = PlottingWindow.Plot1d(self.panel, figsize=(12, 5))
        plotsizer3 = wx.BoxSizer(wx.VERTICAL)
        plotsizer3.Add(self.plot3, 1, wx.EXPAND)
        hbox2.Add(plotsizer3, 1, wx.EXPAND)

        self.plot5 = PlottingWindow.Plot1d(self.panel, figsize=(4, 4))

        plotsizer2 = wx.BoxSizer(wx.VERTICAL)
        plotsizer2.Add(self.plot5, 1, wx.EXPAND)
        hbox2.Add(plotsizer2, 1, wx.EXPAND)


        sizer.Add(hbox2, 1, wx.EXPAND)

        self.panel.SetSizer(sizer)
        self.panel.SetDropTarget(LipDropTarget(self))
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close)

        self.Centre()
        # self.MakeModal(True)
        self.Show(True)
        self.Raise()

    def on_open(self, e=None):
        dlg = wx.FileDialog(self, "Choose a Skyline Output File", '', "", "*.csv")
        if dlg.ShowModal() == wx.ID_OK:
            self.SetStatusText("Opening", number=5)
            filename = dlg.GetPath()
            print("Opening: ", filename)
            self.on_open_file(filename)
        dlg.Destroy()

    def on_open_file(self, path):
        self.path = path
        self.outpath = os.path.splitext(path)[0] + "_filtered.csv"
        self.df = pd.read_csv(path)
        self.process_df()
        self.load_checkboxes()

    def process_df(self):
        self.df = full_processing(self.df)
        self.df = rel_quant(self.df, self.ngroups, self.nreps)
        self.df.to_csv(self.outpath)
        print("Saved to:", self.outpath)

    def load_checkboxes(self):
        self.checkboxes.Clear()
        self.classes = np.unique(self.df["Ontology"])
        print(self.classes)
        self.checkboxes.InsertItems(self.classes, 0)
        self.checkboxes.SetCheckedStrings(self.classes)
        self.update_checks()

    def load_check2(self, e=0, force_refresh=False):
        self.checked_classes2 = self.checkboxes2.GetCheckedStrings()
        self.checkboxes2.Clear()
        self.tailsat = np.unique(self.df2["Tail Unsaturation"]).astype(str)
        self.checkboxes2.InsertItems(self.tailsat, 0)
        self.checked_classes2 = np.intersect1d(self.checked_classes2, self.tailsat)
        if len(self.checked_classes2) != 0 and not force_refresh:
            self.checkboxes2.SetCheckedStrings(self.checked_classes2)
        else:
            self.checkboxes2.SetCheckedStrings(self.tailsat)

    def load_check3(self, e=0, force_refresh=False):
        self.checked_classes3 = self.checkboxes3.GetCheckedStrings()
        self.checkboxes3.Clear()
        self.taillen = np.unique(self.df2["Tail Lengths"]).astype(str)
        self.checkboxes3.InsertItems(self.taillen, 0)
        self.checked_classes3 = np.intersect1d(self.checked_classes3, self.taillen)
        if len(self.checked_classes3) != 0 and not force_refresh:
            self.checkboxes3.SetCheckedStrings(self.checked_classes3)
        else:
            self.checkboxes3.SetCheckedStrings(self.taillen)

    def update_checks(self, e=0, force_refresh2=False, force_refresh3=False):
        self.checked_classes = self.checkboxes.GetCheckedStrings()
        b1 = self.df["Ontology"].isin(self.checked_classes)
        self.df2 = self.df[b1]

        self.load_check2(force_refresh=force_refresh2)
        self.update_check2()
        self.load_check3(force_refresh=force_refresh3)
        self.update_check3()

        plotchoice = self.plotchoices.GetStringSelection()

        self.normchoice = self.normchoices.GetStringSelection()

        if plotchoice == "Class":
            self.colors = self.df2["colors"]
            self.color_class_checkbox()
        elif plotchoice == "Tail Unsaturation":
            if np.amax(self.df2["Tail Unsaturation"]) != 0:
                factor = np.amax(self.df2["Tail Unsaturation"])
            else:
                factor = 1
            self.colors = cm.rainbow(self.df2["Tail Unsaturation"].to_numpy() / factor)
            self.color_unsat_checkbox()
        elif plotchoice == "Tail Length":
            d = self.df2["Tail Lengths"].to_numpy()
            norm = d - np.amin(d)
            if np.amax(norm) != 0:
                norm = norm / np.amax(norm)
            self.colors = cm.rainbow(norm)
            self.color_length_checkbox()
        else:
            if len(np.unique(self.df2["Tail Unsaturation"])) > 1:
                if np.amax(self.df2["Tail Unsaturation"]) != 0:
                    factor = np.amax(self.df2["Tail Unsaturation"])
                else:
                    factor = 1
                self.colors = cm.rainbow(self.df2["Tail Unsaturation"].to_numpy() / factor)
                self.color_unsat_checkbox()
            elif len(np.unique(self.df2["Tail Lengths"])) > 1:
                d = self.df2["Tail Lengths"].to_numpy()
                norm = d - np.amin(d)
                norm = norm / np.amax(norm)
                self.colors = cm.rainbow(norm)
                self.color_length_checkbox()
            else:
                self.colors = self.df2["colors"]
                self.color_class_checkbox()

        self.make_plots()

    def color_class_checkbox(self, e=0, white=False):
        items = self.checkboxes.GetItems()
        for i, item in enumerate(items):
            self.checkboxes.SetItemBackgroundColour(i, (255, 255, 255))
            self.checkboxes.SetItemForegroundColour(i, black_text)
        if not white:
            items = self.checkboxes.GetCheckedItems()
            strings = self.checkboxes.GetCheckedStrings()
            for i, item in enumerate(items):
                s = strings[i]
                b1 = self.df2["Ontology"] == s
                try:
                    c = np.array(self.colors[b1])[0]
                    c = np.array(colors.to_rgba(c))
                    c = np.round(c * 255).astype(int)
                    luminance = ud.get_luminance(c, type=2)
                    if luminance < luminance_cutoff:
                        self.checkboxes.SetItemForegroundColour(i, white_text)
                    self.checkboxes.SetItemBackgroundColour(i, c)
                except:
                    pass
            self.color_length_checkbox(white=True)
            self.color_unsat_checkbox(white=True)
        self.checkboxes.SetFocus()
        pass

    def color_unsat_checkbox(self, e=0, white=False):
        items = self.checkboxes2.GetItems()
        for i, item in enumerate(items):
            self.checkboxes2.SetItemBackgroundColour(i, (255, 255, 255))
            self.checkboxes2.SetItemForegroundColour(i, black_text)
        if not white:
            items = self.checkboxes2.GetCheckedItems()
            strings = self.checkboxes2.GetCheckedStrings()
            for i, item in enumerate(items):
                s = strings[i]
                b1 = self.df2["Tail Unsaturation"] == int(s)
                c = self.colors[b1][0]
                c = np.round(c * 255).astype(int)
                luminance = ud.get_luminance(c, type=2)
                if luminance < luminance_cutoff:
                    self.checkboxes2.SetItemForegroundColour(item, white_text)
                self.checkboxes2.SetItemBackgroundColour(item, c)
            self.color_class_checkbox(white=True)
            self.color_length_checkbox(white=True)

    def color_length_checkbox(self, e=0, white=False):
        items = self.checkboxes3.GetItems()
        for i, item in enumerate(items):
            self.checkboxes3.SetItemBackgroundColour(i, (255, 255, 255))
            self.checkboxes3.SetItemForegroundColour(i, black_text)
        if not white:
            items = self.checkboxes3.GetCheckedItems()
            strings = self.checkboxes3.GetCheckedStrings()
            for i, item in enumerate(items):
                s = strings[i]
                b1 = self.df2["Tail Lengths"] == int(s)
                c = self.colors[b1][0]
                c = np.round(c * 255).astype(int)
                luminance = ud.get_luminance(c, type=2)
                if luminance < luminance_cutoff:
                    self.checkboxes3.SetItemForegroundColour(item, white_text)
                self.checkboxes3.SetItemBackgroundColour(item, c)
            self.color_class_checkbox(white=True)
            self.color_unsat_checkbox(white=True)

    def update_check2(self, e=0):
        self.checked_classes2 = np.array(self.checkboxes2.GetCheckedStrings()).astype(int)
        b1 = self.df2["Tail Unsaturation"].isin(self.checked_classes2)
        self.df2 = self.df2[b1]

    def update_check3(self, e=0):
        self.checked_classes3 = np.array(self.checkboxes3.GetCheckedStrings()).astype(int)
        b1 = self.df2["Tail Lengths"].isin(self.checked_classes3)
        self.df2 = self.df2[b1]

    def check_only(self, e=0):
        selection = self.checkboxes.GetSelection()
        checked = self.checkboxes.GetCheckedItems()

        if len(checked) == 1 and checked[0] == selection:
            for i, item in enumerate(self.checkboxes.GetItems()):
                self.checkboxes.Check(i, check=True)
        else:
            for i, item in enumerate(self.checkboxes.GetItems()):
                if i == selection:
                    self.checkboxes.Check(i, check=True)
                else:
                    self.checkboxes.Check(i, check=False)
        self.update_checks(force_refresh2=True, force_refresh3=True)

    def check_only2(self, e=0):
        selection = self.checkboxes2.GetSelection()
        checked = self.checkboxes2.GetCheckedItems()

        if len(checked) == 1 and checked[0] == selection:
            for i, item in enumerate(self.checkboxes2.GetItems()):
                self.checkboxes2.Check(i, check=True)
            force = True
        else:
            force = True
            for i, item in enumerate(self.checkboxes2.GetItems()):
                if i == selection:
                    self.checkboxes2.Check(i, check=True)
                else:
                    self.checkboxes2.Check(i, check=False)
        self.update_checks(force_refresh3=True)

    def check_only3(self, e=0):
        selection = self.checkboxes3.GetSelection()
        checked = self.checkboxes3.GetCheckedItems()

        if len(checked) == 1 and checked[0] == selection:
            for i, item in enumerate(self.checkboxes3.GetItems()):
                self.checkboxes3.Check(i, check=True)
            force = True
        else:
            for i, item in enumerate(self.checkboxes3.GetItems()):
                if i == selection:
                    self.checkboxes3.Check(i, check=True)
                else:
                    self.checkboxes3.Check(i, check=False)
            force = False
        self.update_checks(force_refresh2=True)

    def make_plots(self):
        self.normname = norms[normchoices[self.normchoice]]
        self.rt_plot()
        self.md_plot()
        self.fold_plot()
        self.make_plot_4()
        self.make_plot_5()

    def rt_plot(self, e=0):
        mz = self.df2["mass"].to_numpy()
        RT = self.df2["Average Rt(min)"].to_numpy()
        self.plot1.scatterplottop(mz, RT, xlabel="Mass (Da)", ylabel="Retention Time (min)", color=self.colors)

    def md_plot(self, e=0):
        mz = self.df2["mass"].to_numpy()
        md = self.df2["kmassd"].to_numpy()
        self.plot2.scatterplottop(mz, md, xlabel="Mass (Da)", ylabel="Mass Defect", color=self.colors)

    def fold_plot(self, e=0):
        labels = self.df2["Simp Name"].to_numpy()
        fold = self.df2["FoldChange"+self.normname].to_numpy()
        sd = self.df2["SDFold"+self.normname].to_numpy()
        # colors = self.df2["colors"].to_numpy()

        if len(labels) == 0:
            self.plot3.clear_plot()
        else:
            self.plot3._axes = [0.11, 0.3, 0.8, 0.65]
            self.plot3.barplottoperrors(np.arange(len(labels)), fold, labels, yerr=sd, colortab=self.colors,
                                        repaint=False,
                                        ylabel="Fold Change")
            self.plot3.plotadd([-0.5, len(labels)], [1, 1], linestyle="--", nopaint=False)

    def make_plot_4(self, e=0):
        """
        fold = self.df2["FoldChange"+self.normname].to_numpy()
        avga = self.df2["AvgA_0"].to_numpy()
        self.plot4.scatterplottop(avga, fold, color = self.colors)"""
        labels = self.df2["Simp Name"].to_numpy()
        avga = self.df2["AvgA" + self.normname + "_1"].to_numpy()
        sd = self.df2["SD" + self.normname + "_1"].to_numpy()
        # colors = self.df2["colors"].to_numpy()

        if len(labels) == 0:
            self.plot4.clear_plot()
        else:
            self.plot4._axes = [0.11, 0.3, 0.8, 0.65]
            self.plot4.barplottoperrors(np.arange(len(labels)), avga, labels, yerr=sd, colortab=self.colors,
                                        repaint=True,
                                        ylabel="Area")
            # self.plot5.plotadd([-0.5, len(labels)], [1, 1], linestyle="--", nopaint=False)
        pass

    def make_plot_5(self, e=0):
        labels = self.df2["Simp Name"].to_numpy()
        avga = self.df2["AvgA"+self.normname+"_0"].to_numpy()
        sd = self.df2["SD"+self.normname+"_0"].to_numpy()
        # colors = self.df2["colors"].to_numpy()

        if len(labels) == 0:
            self.plot5.clear_plot()
        else:
            self.plot5._axes = [0.11, 0.3, 0.8, 0.65]
            self.plot5.barplottoperrors(np.arange(len(labels)), avga, labels, yerr=sd, colortab=self.colors,
                                        repaint=True,
                                        ylabel="Area")
            #self.plot5.plotadd([-0.5, len(labels)], [1, 1], linestyle="--", nopaint=False)
        pass

    def on_motion(self, xpos, ypos):
        try:
            if xpos is not None and ypos is not None:
                self.SetStatusText("x=%.4f y=%.2f" % (xpos, ypos), number=2)
        except:
            pass
        try:
            if xpos is not None and ypos is not None:
                closest_thing = self.find_mouse_species(xpos, ypos)
                self.SetStatusText("Species: " + closest_thing, number=1)
        except:
            pass

    def find_mouse_species(self, xpos, ypos):
        active = self.find_active_plot()
        if xpos is not None and ypos is not None:
            if active == 0:
                mz = xpos
                rt = ypos
                mzarray = self.df2["mass"].to_numpy()
                rtarray = self.df2["Average Rt(min)"].to_numpy()

                dist = np.sqrt(np.square(mzarray - mz) + np.square(rtarray - rt))
                argmin = np.argmin(dist)
                species = self.df2.iloc[argmin]["Molecule Name"]
                return species
            elif active == 1:
                mz = xpos
                rt = ypos
                mzarray = self.df2["mass"].to_numpy()
                rtarray = self.df2["kmassd"].to_numpy()

                dist = np.sqrt(np.square(mzarray - mz) + np.square(rtarray - rt))
                argmin = np.argmin(dist)
                species = self.df2.iloc[argmin]["Molecule Name"]
                return species
            elif active == 2:
                argmin = round(xpos)
                species = self.df2.iloc[argmin]["Molecule Name"]
                return species
        return ""

    def find_active_plot(self):
        plotlist = [self.plot1, self.plot2, self.plot3]
        for i, p in enumerate(plotlist):
            if p.mouse_active:
                return i

    def on_close(self, e):
        """
        Close the window.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        # self.MakeModal(False)

    def on_save_fig(self, e):
        """
        Saves the figures in self.directory as PNGs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, self.outfname + "_1.png")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1

    def on_save_fig_pdf(self, e):
        """
        Saves the figures in self.directory as PDFs.
        :param e: Unused event
        :return: None
        """
        name1 = os.path.join(self.directory, self.outfname + "_1.pdf")
        if self.plot1.flag:
            self.plot1.on_save_fig(e, name1)
            # print name1


class LipDropTarget(wx.FileDropTarget):
    """"""

    def __init__(self, window):
        """Constructor"""
        wx.FileDropTarget.__init__(self)
        self.window = window

    def OnDropFiles(self, x, y, filenames):
        """
        When files are dropped, either open a single file or run in batch.
        """
        if len(filenames) == 1:
            path = filenames[0]
            self.window.on_open_file(path)
        elif len(filenames) > 1:
            for f in filenames:
                self.window.on_open_file(f)
        else:
            print("Error in drag and drop.")
        return 0


# Main App Execution
if __name__ == "__main__":
    app = wx.App(False)
    frame = LipyDecWindow(None)

    file = "C:\Data\Lipidomics\Exported_TLs\\220816_MS1_Ecoli_lipid_ND_res_outputs_1.csv"
    # file = "Z:\Group Share\Melanie Odenkirk\\20220906_EmptyND_Control\Exported_TLs_control_220906\C_comparison_NDminus_NDplus_220906_CVfilt_CL1_filtered.csv"
    # file = "Z:\Group Share\Melanie Odenkirk\\20220815_Ecoli_d7_lowerCone\Exported_TLs\\220822_MS1_avg_B1_ND_output.csv"
    file = "C:\Data\Lipidomics\\Exported_TLs_control_220906\C1_comparison_220915.csv"
    if os.path.isfile(file):
        frame.on_open_file(file)
        # frame.checked_classes=["PE"]
        # frame.make_plots()
    app.MainLoop()
