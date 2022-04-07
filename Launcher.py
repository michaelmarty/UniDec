import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

import wx
from unidec_modules.unidec_presbase import UniDecPres
from unidec_modules.unidec_enginebase import UniDecEngine
import multiprocessing
from CDPres import UniDecCDApp
from GUniDec import UniDecApp

import datacollector
from metaunidec import mudpres
from import_wizard import ImportWizard
from metaunidec.meta_import_wizard.meta_import_wizard import ImportWizard as HDF5Wizard
from metaunidec.ultrameta import DataCollector as UMDC
from UniChrom2 import ChromApp
import wx.py as py
import os
import sys

class UniDecLauncher(UniDecPres):
    """
    Main UniDec GUI Application.

    Presenter contains UniDec engine at self.eng and main GUI window at self.view
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize App
        :param args:
        :param kwargs:
        :return: UniDecApp object
        """
        UniDecPres.__init__(self, *args, **kwargs)
        self.init(*args, **kwargs)

    def init(self, *args, **kwargs):
        self.view = Lview(self)
        self.view.Bind(wx.EVT_CLOSE, self.on_close)

        if "--meta" in sys.argv[1:] or "-m" in sys.argv[1:]:
            print("Launching Meta")
            self.view.button4()

        if "--chrom" in sys.argv[1:] or "-c" in sys.argv[1:]:
            print("Launching UniChrom")
            self.view.button8()

        if "--unidec" in sys.argv[1:] or "-u" in sys.argv[1:]:
            print("Launching UniDec")
            self.view.button1()

        if len(sys.argv)>1:
            self.view.button8()

    def on_close(self, e=None):
        self.quit_application()
        self.view.Destroy()
        sys.exit()

class Lview(wx.Frame):
    def __init__(self, parent):
        self.eng = UniDecEngine()
        wx.Frame.__init__(self, None, title="UniDec Launcher Version " + str(self.eng.version))

        sizer = wx.GridBagSizer(wx.HORIZONTAL)
        panel = wx.Panel(self)
        button1 = wx.Button(panel, -1, "UniDec\n\nDeconvolve MS and IM-MS")
        button2 = wx.Button(panel, -1, "Data Collector\n\nVisualize multiple spectra\nExtract Trends\nFit Kd's")
        button3 = wx.Button(panel, -1, "Import Wizard\n\nBatch convert Waters Raw to Txt")
        button4 = wx.Button(panel, -1, "MetaUniDec\n\nBatch process and visualize MS spectra")
        button5 = wx.Button(panel, -1, "UniDec API Shell\n\nScript UniDec with console")
        button6 = wx.Button(panel, -1, "HDF5 Import Wizard\n\nImport Data into HDF5 for MetaUniDec")
        button7 = wx.Button(panel, -1, "UltraMeta Data Collector\n\nVisualize Multiple HDF5 Data Sets\nFit Trends")
        button8 = wx.Button(panel, -1, "UniChrom2\n\nDeconvolution of Chromatograms\nUniDec for LC/MS Data")
        button9 = wx.Button(panel, -1, "UniDecCD\n\nDeconvolution of Charge Detection MS\nUniDec for CD-MS Data")

        html = wx.html.HtmlWindow(panel, -1, size=(390,310))
        pathtofile = os.path.dirname(os.path.abspath(__file__))
        self.imagepath = os.path.join(pathtofile, "UniDecLogoMR.png")
        #print(self.imagepath)
        html.SetPage(
            "<html><body>"
            #"<h1>UniDec</h1>"
            "<img src=\"" + self.imagepath +"\" alt=\"PNG Icon\" height=\"250\" width=\"363\">"
            "<p>Please Cite: Marty et al. Anal. Chem. 2015. " 
            "DOI: 10.1021/acs.analchem.5b00140.</p>"
            "</body></html>"
        )

        sizer.Add(button1, (0, 0), flag=wx.EXPAND)
        sizer.Add(button2, (1, 0), flag=wx.EXPAND)
        sizer.Add(button3, (2, 0), flag=wx.EXPAND)
        sizer.Add(button4, (0, 1), flag=wx.EXPAND)
        sizer.Add(button6, (2, 1), flag=wx.EXPAND)
        sizer.Add(button7, (1, 1), flag=wx.EXPAND)
        sizer.Add(button9, (3, 0), flag=wx.EXPAND)
        sizer.Add(button5, (4, 0), span=(1,2), flag=wx.EXPAND)
        sizer.Add(button8, (3, 1), flag=wx.EXPAND)
        sizer.Add(html, (0, 2), span=(5, 2))

        self.Bind(wx.EVT_BUTTON, self.button1, button1)
        self.Bind(wx.EVT_BUTTON, self.button2, button2)
        self.Bind(wx.EVT_BUTTON, self.button3, button3)
        self.Bind(wx.EVT_BUTTON, self.button4, button4)
        self.Bind(wx.EVT_BUTTON, self.button5, button5)
        self.Bind(wx.EVT_BUTTON, self.button6, button6)
        self.Bind(wx.EVT_BUTTON, self.button7, button7)
        self.Bind(wx.EVT_BUTTON, self.button8, button8)
        self.Bind(wx.EVT_BUTTON, self.button9, button9)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Centre()
        self.Show(True)

    def button1(self, e=None):
        print("Launching UniDec")
        app = UniDecApp()
        app.start()

    def button2(self, e=None):
        print("Launching Data Collector")
        app = wx.App(False)
        frame = datacollector.DataCollector(None, "Collect Data")
        app.MainLoop()

    def button3(self, e=None):
        print("Launching Waters Converter Wizard")
        app = wx.App(False)
        frame = ImportWizard(None)
        frame.Show()
        app.MainLoop()

    def button4(self, e=None):
        print("Launching MetaUniDec")
        app = mudpres.UniDecApp()
        app.start()

    def button5(self, e=None):
        print("Launching Scripting Shell")
        app = Shell()
        app.start()

    def button6(self, e=None):
        print("Launching HDF5 Import Wizard")
        app = wx.App(False)
        frame = HDF5Wizard(None)
        frame.Show()
        app.MainLoop()

    def button7(self, e=None):
        print("Launching UltraMeta Data Collector")
        app = wx.App(False)
        frame = UMDC(None, "UltraMeta Data Collector")
        app.MainLoop()

    def button8(self, e=None):
        print("Launching UniChrom2")
        app = ChromApp()
        app.start()

    def button9(self, e=None):
        print("Launching UniDecCD")
        app = UniDecCDApp()
        app.start()

class Shell(object):
    def __init__(self, *args, **kwargs):
        self.__wx_app = wx.App(redirect=True)

        self.shell = py.shell.Shell(wx.Frame(None))

        self.shellwindow = py.shell.ShellFrame(self.shell, title="UniDecShell").Show()

        # self.shell.Execute('app=UniDecApp()')
        # self.shell.Execute('app.start()')
        # self.shellwindow.Center()
        # self.shell.setFocus()
        self.__wx_app.MainLoop()


if __name__ == '__main__':
    # app2 = Shell()
    multiprocessing.freeze_support()
    app = UniDecLauncher(sys.argv[1:])
    app.start()
