import wx
import os
from modules.plotting import PlottingWindow
from glycodec_eng import GlycoDecEngine, enzyme_options
from glycan_examples import *

__author__ = 'Michael.Marty'


class GlycoDecWindow(wx.Frame):
    def __init__(self, parent, engine=None):
        wx.Frame.__init__(self, parent, title="Glycan Builder")  # ,size=(-1,-1))
        statusbar_log_silencer = wx.LogNull()
        self.CreateStatusBar(3)
        self.SetStatusWidths([-1, -3, -3])
        del statusbar_log_silencer
        if engine is None:
            self.eng = GlycoDecEngine()
        else:
            self.eng = engine

        # Setup initial values
        try:
            self.directory = parent.directory
        except:
            self.directory = os.getcwd()

        self.outfname = "Glycofig"
        # Make the menu
        '''
        filemenu = wx.Menu()
        menu_oligomer_tools = filemenu.Append(wx.ID_ANY, "Generate from Masses", "Oligomer and Mass Tools")
        filemenu.AppendSeparator()
        menu_save_fig_png = filemenu.Append(wx.ID_ANY, "Save Figures as PNG",
                                            "Save all figures as PNG in central directory")
        menu_save_fig_pdf = filemenu.Append(wx.ID_ANY, "Save Figures as PDF",
                                            "Save all figures as PDF in central directory")
        self.Bind(wx.EVT_MENU, self.on_oligomer_tools, menu_oligomer_tools)
        self.Bind(wx.EVT_MENU, self.on_save_fig, menu_save_fig_png)
        self.Bind(wx.EVT_MENU, self.on_save_fig_pdf, menu_save_fig_pdf)

        self.plotmenu = wx.Menu()
        self.menuaddline = self.plotmenu.Append(wx.ID_ANY, "Add Horizontal Line",
                                                "Add Horizontal Line at Specific Y Value")
        self.plotmenu.AppendSeparator()
        self.menutotal = self.plotmenu.Append(wx.ID_ANY, "Plot Total Bound vs. Unbound",
                                              "Assumes the first value is the unbound and each other is 1, 2, 3, and so on bound ligands.")
        self.Bind(wx.EVT_MENU, self.on_add_line, self.menuaddline)
        self.Bind(wx.EVT_MENU, self.on_total, self.menutotal)

        menu_bar = wx.MenuBar()
        menu_bar.Append(filemenu, "&File")
        menu_bar.Append(self.plotmenu, "Plot")
        self.SetMenuBar(menu_bar)'''

        # Setup the GUI
        panel = wx.Panel(self)

        sizer = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        sb = wx.StaticBox(panel, label='GlycoCT')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)

        self.glytextbox = wx.TextCtrl(sb, style=wx.TE_MULTILINE, size=(210, 520))
        sbs.Add(self.glytextbox)

        '''
        importbutton = wx.Button(panel, label="Import from File")
        sbs.Add(importbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_import_masses, importbutton)
        
        addbutton = wx.Button(panel, label="Add Species")
        sbs.Add(addbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_add_mass, addbutton)
        
        self.masslistbox = masstools.MassListCtrlPanel(panel, columntitle="Mass Defect Value", size=(210, 320))
        sbs.Add(wx.StaticText(panel, label="Mass Defect List"))
        sbs.Add(self.masslistbox)
        '''
        updatebutton = wx.Button(sb, label="Update Plot")
        sbs.Add(updatebutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_update_button, updatebutton)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2.Add(wx.StaticText(sb, label="Mass (Da):"))
        self.masstextbox = wx.TextCtrl(sb, style=wx.TE_READONLY)
        hbox2.Add(self.masstextbox)
        sbs.Add(hbox2)
        # hbox.Add(sbs, 0, wx.EXPAND)
        self.plot1 = PlottingWindow.Plot1d(panel, figsize=(4, 4))

        plotsizer1 = wx.BoxSizer(wx.VERTICAL)
        plotsizer1.Add(self.plot1, 1, wx.EXPAND)

        hbox.Add(sbs, 0)
        hbox.Add(plotsizer1, 1, wx.EXPAND)

        sb2 = wx.StaticBox(panel, label='Digested Version')
        sbs2 = wx.StaticBoxSizer(sb2, orient=wx.VERTICAL)

        self.glytextbox2 = wx.TextCtrl(sb2, style=wx.TE_MULTILINE | wx.TE_READONLY, size=(210, 520))
        sbs2.Add(self.glytextbox2)

        controlsizer = wx.BoxSizer(wx.HORIZONTAL)
        controlsizer.Add(wx.StaticText(sb2, label="Enzyme Type: "), 0, wx.ALIGN_CENTER_VERTICAL)
        self.ctletype = wx.ComboBox(sb2, value=enzyme_options[0], choices=list(enzyme_options.values()),
                                    style=wx.CB_READONLY | wx.ALIGN_CENTER_VERTICAL)
        controlsizer.Add(self.ctletype, 0, wx.ALIGN_CENTER_VERTICAL)

        digestbutton = wx.Button(sb2, label="Digest")
        controlsizer.Add(digestbutton, 0, wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.on_digest, digestbutton)

        sbs2.Add(controlsizer)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        hbox2.Add(wx.StaticText(sb2, label="Mass (Da):"))
        self.masstextbox2 = wx.TextCtrl(sb2, style=wx.TE_READONLY)
        hbox2.Add(self.masstextbox2)
        sbs2.Add(hbox2)

        hbox.Add(sbs2, 0, wx.EXPAND)

        self.plot2 = PlottingWindow.Plot1d(panel, figsize=(4, 4))

        plotsizer2 = wx.BoxSizer(wx.VERTICAL)
        plotsizer2.Add(self.plot2, 1, wx.EXPAND)
        hbox.Add(plotsizer2, 1, wx.EXPAND)

        sizer.Add(hbox, 1, wx.EXPAND)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE, self.on_close)

        self.Centre()
        # self.MakeModal(True)
        self.Show(True)
        self.Raise()

        self.load_text(g3)
        self.on_update_button()

    def update(self):
        self.glytext = self.glytextbox.GetValue()
        # print(self.glytext)

    def load_text(self, text):
        self.glytextbox.SetValue(text)


    def load_text2(self, text):
        self.glytextbox2.SetValue(text)


    def on_update_button(self, e=None):
        self.update()
        self.eng.load_glycoct(self.glytext)
        self.glyco_plot()
        self.update_mass()

    def update_mass(self):
        mass = self.eng.get_mass()
        self.masstextbox.SetValue(str(mass))
        self.update_echoice()

    def update_mass2(self):
        mass = self.eng.get_mass()
        self.masstextbox2.SetValue(str(mass))
        self.update_echoice()

    def glyco_plot(self, e=None, plot=None):
        if plot is None:
            plot = self.plot1
        axis = plot.get_blank_axis(scale=[0,0,1,1])
        self.eng.plot_struct(axis)
        plot.repaint()
        #self.plot1.setup_zoom([axis], "box")
        pass

    def update_echoice(self, e=None):
        if self.eng.structure is not None:
            choices = self.eng.find_active_enzymes()
            self.ctletype.SetItems(choices)
            self.ctletype.SetSelection(0)
            if len(self.eng.digest_history)>0:
                dhist = "Parent"
                for s in self.eng.digest_history:
                    dhist = dhist + " -> " + s
                self.SetStatusText(dhist, 2)

    def on_digest(self, e=None):
        etype = self.ctletype.GetStringSelection()
        print("Digesting with:", etype)
        self.eng.enzyme(etype)
        self.load_text2(str(self.eng.structure))
        self.glyco_plot(plot=self.plot2)
        self.update_mass2()

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



# Main App Execution
if __name__ == "__main__":
    app = wx.App(False)
    frame = GlycoDecWindow(None)
    app.MainLoop()
