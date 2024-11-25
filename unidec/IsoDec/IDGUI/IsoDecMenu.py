import wx
import unidec.modules.isolated_packages.preset_manager as pm
import numpy as np
import os


class main_menu(wx.Menu):
    # noinspection PyMissingConstructor
    def __init__(self, parent, config, pres):
        super(wx.Menu, self).__init__()
        self.pres = pres
        self.config = config
        self.parent = parent

        self.filemenu = wx.Menu()
        self.toolsmenu = wx.Menu()
        self.analysismenu = wx.Menu()
        self.advancedmenu = wx.Menu()
        self.experimentalmenu = wx.Menu()
        self.menuOpenRecent = wx.Menu()

        # File Menu
        self.menuOpen = self.filemenu.Append(
            wx.ID_OPEN,
            "Open File (Text, mzML, or Thermo RAW)\tCtrl+O",
            " Open a Text File in x y text, mzML, or Thermo RAW format",
        )
        self.menuOpenRaw = self.filemenu.Append(

            wx.ID_ANY,
            "Open Waters or Agilent File",
            " Open a Waters .Raw or Agilent .D File",
        )
        self.filemenu.AppendSubMenu(self.menuOpenRecent, "Open Recent File")
        self.filemenu.AppendSeparator()

        self.menupastespectrum = self.filemenu.Append(
            wx.ID_ANY,
            "Get Spectrum From Clipboard\tCtrl+G",
            "Pastes the spectrum, formats, and loads",
        )
        self.filemenu.AppendSeparator()

        # Menu item to reset the config file to default
        self.menuResetConfig = self.filemenu.Append(
            wx.ID_ANY,
            "Set Settings to Default",
            "Resets the settings to default",
        )
        self.parent.Bind(wx.EVT_MENU, self.pres.on_init_config, self.menuResetConfig)
        self.filemenu.AppendSeparator()

        self.menufigdialog = self.filemenu.Append(
            wx.ID_ANY,
            "Save Figure As",
            "Dialog to define path and extension for figure saving",
        )
        # Figure Submenu
        self.figmenu = wx.Menu()

        self.menuSaveFigure0 = self.figmenu.Append(
            wx.ID_ANY, "Save Figures as .pdf", "Save all figures to .pdf format"
        )
        self.menuSaveFigure1s = self.figmenu.Append(
            wx.ID_ANY,
            "Save Figures as PDF thumbnail",
            "Save all figures to small PDF format",
        )
        self.menuSaveFigure1 = self.figmenu.Append(
            wx.ID_ANY, "Save Figures as .eps", "Save all figures to .eps format"
        )
        self.menuSaveFigure2 = self.figmenu.Append(
            wx.ID_ANY, "Save Figures as .png", "Save all figures to .png format"
        )

        self.filemenu.AppendSubMenu(self.figmenu, "Save Figure Presets")
        self.filemenu.AppendSeparator()

        # Example Data
        self.examplemenu, self.masterd2 = pm.make_preset_menu(
            os.path.join(self.config.exampledatadir, "IsoDec"),
            exclude_dir="_unidecfiles",
            topi=2500,
            exclude_ext="hdf5",
        )

        keys = []
        for i, d in enumerate(self.masterd2):
            if i < 10:
                keys.append([str(i + 1), self.pres.on_ex, d[2]])
        self.menukeys = keys

        self.filemenu.AppendSubMenu(self.examplemenu, "Load Example Data")
        for i, path, item in self.masterd2:
            # print(i, path, item)
            self.parent.Bind(wx.EVT_MENU, self.on_example_data, item)
        # ..............
        self.filemenu.AppendSeparator()
        self.menuAbout = self.filemenu.Append(
            wx.ID_ABOUT, "&About", " Information about this program"
        )
        self.menuExit = self.filemenu.Append(
            wx.ID_EXIT, "E&xit\tCtrl+Q", " Terminate the Program"
        )

        # Set Events for Menu Bar

        # File Menu
        self.parent.Bind(wx.EVT_MENU, self.pres.on_open, self.menuOpen)
        self.parent.Bind(wx.EVT_MENU, self.pres.on_raw_open, self.menuOpenRaw)

        self.parent.Bind(
            wx.EVT_MENU, self.pres.on_paste_spectrum, self.menupastespectrum
        )

        self.parent.Bind(
            wx.EVT_MENU, self.parent.on_save_figure_dialog, self.menufigdialog
        )
        self.parent.Bind(
            wx.EVT_MENU, self.parent.on_save_figure_pdf, self.menuSaveFigure0
        )
        self.parent.Bind(
            wx.EVT_MENU, self.parent.on_save_figure_eps, self.menuSaveFigure1
        )
        self.parent.Bind(
            wx.EVT_MENU, self.parent.on_save_figure_small, self.menuSaveFigure1s
        )
        self.parent.Bind(
            wx.EVT_MENU, self.parent.on_save_figure_png, self.menuSaveFigure2
        )
        self.parent.Bind(wx.EVT_MENU, self.parent.on_about, self.menuAbout)
        self.parent.Bind(wx.EVT_MENU, self.parent.on_exit, self.menuExit)

        # Add batch process to Tools Menu
        self.menuBatch = self.toolsmenu.Append(
            wx.ID_ANY, "Batch Process Full File", " Batch process all scans in a full file"
        )
        # Bind to self.parent.on_batch
        self.parent.Bind(wx.EVT_MENU, self.pres.on_batch, self.menuBatch)

        self.menuOpenDir = self.advancedmenu.Append(wx.ID_ANY, "Open Saved File Directory",
                                                    "Opens the save directory in the file explorer")
        self.parent.Bind(wx.EVT_MENU, self.parent.on_open_dir, self.menuOpenDir)

        # Setting Menu Bar
        self.menuBar = wx.MenuBar()
        self.menuBar.Append(self.filemenu, "&File")
        self.menuBar.Append(self.toolsmenu, "Tools")
        # self.menuBar.Append(self.analysismenu, "Analysis")
        self.menuBar.Append(self.advancedmenu, "Advanced")
        # self.menuBar.Append(self.experimentalmenu, "Experimental")
        self.parent.SetMenuBar(self.menuBar)

    def on_example_data(self, e):
        # print("Clicked", e)
        try:
            nid = e.GetId()
            ids = self.masterd2[:, 0].astype(float)
            pos = np.argmin(np.abs(ids - nid))
            self.load_example_data(pos)
        except Exception as error:
            try:
                nid = e
                ids = self.masterd2[:, 0].astype(float)
                pos = np.argmin(np.abs(ids - nid))
                self.load_example_data(pos)
            except Exception as error2:
                print(error)
                print(error2)

    def load_example_data(self, pos):
        path = self.masterd2[pos, 1]
        dir = os.path.dirname(path)
        file = os.path.split(path)[1]
        print("Opening Path:", path, file, dir)
        self.pres.on_open_file(file, dir)

    def update_recent(self):
        menu_list = self.menuOpenRecent.GetMenuItems()
        for i in range(len(menu_list) - 1, -1, -1):  # clear menu
            self.menuOpenRecent.DestroyItem(menu_list[i])

        max_items = 5  # can be changed to whatever
        added = 0
        for file_path in self.pres.recent_files:
            if added >= max_items:
                break
            filename = os.path.basename(file_path)
            if os.path.splitext(filename)[1] != ".hdf5":
                self.add_to_recent(file_path)
                added += 1

    def add_to_recent(self, file_path):
        # This needs to be separate from update_recent() for binding to work
        filename = os.path.basename(file_path)
        dirname = os.path.dirname(file_path)
        new_item = self.menuOpenRecent.Append(wx.ID_ANY, filename)
        self.parent.Bind(
            wx.EVT_MENU, lambda e: self.pres.on_open_file(filename, dirname), new_item
        )
