import wx
import pandas as pd
from unidec.modules.isolated_packages.spreadsheet import *
from unidec.batch import UniDecBatchProcessor as BPEngine


class UPPApp(wx.Frame):
    """"""

    def __init__(self, nrows=1, ncolumns=4, title="UniDec Pharma Pipeline"):
        """Constructor"""
        wx.Frame.__init__(self, parent=None, title=title, size=(1800, 600))
        self.use_decon = True
        self.use_converted = True
        self.bpeng = BPEngine()

        menu = wx.Menu()
        # Open File Menu
        open_file_menu_item = menu.Append(wx.ID_ANY, "Open File", "Open a CSV or Excel file")
        self.Bind(wx.EVT_MENU, self.on_load_file, open_file_menu_item)

        # Create the menubar
        menuBar = wx.MenuBar()
        menuBar.Append(menu, "&File")
        self.SetMenuBar(menuBar)

        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        # Insert a button and bind it with a handler called on_run
        btn = wx.Button(panel, label="Run")
        btn.Bind(wx.EVT_BUTTON, self.on_run)
        hsizer.Add(btn, 0)

        # Insert a button for Open All HTML Reports and bind to function
        btn = wx.Button(panel, label="Open All HTML Reports")
        btn.Bind(wx.EVT_BUTTON, self.on_open_all_html)
        hsizer.Add(btn, 0)

        # Insert a static text of directory
        # hsizer.Add(wx.StaticText(panel, label="   Data Directory:", style=wx.ALIGN_CENTER_VERTICAL))
        # Insert a text box to read out the directory
        # self.dirtxtbox = wx.TextCtrl(panel, size=(400, -1))
        # hsizer.Add(self.dirtxtbox, 0, wx.EXPAND)
        # Add a button to set the directory
        # btn = wx.Button(panel, label="...")

        # Insert a static text of tolerance
        #hsizer.Add(wx.StaticText(panel, label="   Tolerance:", style=wx.ALIGN_CENTER_VERTICAL))
        # Insert a text box to read out the directory
        #self.tolbox = wx.TextCtrl(panel, size=(50, -1))
        #self.tolbox.SetValue("50")
        #hsizer.Add(self.tolbox, 0, wx.EXPAND)
        #hsizer.Add(wx.StaticText(panel, label="Da   ", style=wx.ALIGN_CENTER_VERTICAL))

        # Insert a checkbox to select whether to use already converted data
        self.useconvbox = wx.CheckBox(panel, label="Use Converted Data")
        hsizer.Add(self.useconvbox, 0, wx.EXPAND)
        self.useconvbox.SetValue(self.use_converted)

        # Insert a checkbox to select whether to use already deconvolved data
        self.usedeconbox = wx.CheckBox(panel, label="Deconvolve Data")
        hsizer.Add(self.usedeconbox, 0, wx.EXPAND)
        self.usedeconbox.SetValue(self.use_decon)

        sizer.Add(hsizer, 0, wx.ALL | wx.EXPAND)

        self.ss = SpreadsheetPanel(self, panel, nrows, ncolumns).ss
        sizer.Add(self.ss, 1, wx.EXPAND)
        panel.SetSizer(sizer)
        self.Show()

    def on_run(self, event=None):
        print("Run button pressed")
        self.get_from_gui()
        self.bpeng.run_df(decon=self.use_decon, use_converted=self.use_converted)
        self.ss.set_df(self.bpeng.rundf)

    def load_file(self, filename):
        print("Loading File:", filename)
        self.bpeng.top_dir = os.path.dirname(filename)
        df = file_to_df(filename)
        self.ss.set_df(df)
        # dirname = os.path.dirname(filename)
        # self.set_dir_tet_box(dirname)

    def on_load_file(self, event):
        print("Load button pressed")
        # Create a file dialog
        with wx.FileDialog(self, "Open CSV or Excel File",
                           wildcard="CSV or Excel files (*.csv; *.xlsx; *.xls)|*.csv; *.xlsx; *.xls|"
                                    "CSV files (*.csv)|*.csv|"
                                    "Excel files (*.xlsx; *.xls)|*.xlsx; *.xls",
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:
            # Show the dialog and retrieve the user response. If it is the OK response,
            # process the data.
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return
            # Proceed loading the file chosen by the user
            pathname = fileDialog.GetPath()
            self.load_file(pathname)

    #def set_dir_tet_box(self, dirname):
    #    self.dirtxtbox.SetValue(dirname)

    def get_from_gui(self):
        self.use_converted = self.useconvbox.GetValue()
        self.use_decon = self.usedeconbox.GetValue()

        # dirname = self.dirtxtbox.GetValue()
        # tol = self.tolbox.GetValue()
        # self.bpeng.data_dir = dirname
        #try:
        #    self.bpeng.tolerance = float(tol)
        #except Exception as e:
        #    print("Error with Tolerance Value. Using default value of 50 Da", e)
        #    self.bpeng.tolerance = 10
        #    self.tolbox.SetValue("10")

        self.ss.remove_empty()
        ssdf = self.ss.get_df()
        self.bpeng.rundf = ssdf

    def on_open_all_html(self, event):
        print("Open All HTML Reports button pressed")
        self.bpeng.open_all_html()

    def open_unidec(self, row):
        print("Opening in UniDec:", row)

    def on_exit(self, event):
        self.Close()


if __name__ == "__main__":
    app = wx.App()
    frame = UPPApp(12, 8)
    frame.usedeconbox.SetValue(False)
    path = "C:\\Data\\Wilson_Genentech\\sequences_short.xlsx"
    frame.load_file(path)
    # frame.set_dir_tet_box("C:\\Data\\Wilson_Genentech\\Data")
    # print(df)
    frame.on_run()

    app.MainLoop()
