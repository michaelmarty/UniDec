import warnings
import os
import sys

warnings.simplefilter(action="ignore", category=FutureWarning)
# psims warns about optional Blosc support for mzMLb writing, which UniDec does not use.
warnings.filterwarnings("ignore", message="hdf5plugin is missing", category=UserWarning)

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import wx
import wx.html
import multiprocessing
import threading

from unidec.modules.unidec_presbase import UniDecPres
from unidec.modules.unidec_enginebase import UniDecEngine

from unidec.IsoDecGUI import IsoDecPres
from unidec.UniDecCD import UniDecCDApp
from unidec.GUniDec import UniDecApp
from unidec.UniDecIM import UniDecIMApp

import unidec.DataCollector as datacollector
from unidec import MetaUniDec as mudpres
from unidec.ImportWizard import ImportWizard
from unidec.metaunidec.meta_import_wizard.meta_import_wizard import ImportWizard as HDF5Wizard
from unidec.metaunidec.ultrameta import DataCollector as UMDC
from unidec.UniChrom import ChromApp
from unidec.UPP import UPPApp
from unidec.modules import unidecstructure
from unidec.UniChromCD import UniChromCDApp
import locale

locale.setlocale(locale.LC_ALL, 'C')


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
            return

        if "--chrom" in sys.argv[1:] or "-c" in sys.argv[1:]:
            print("Launching UniChrom")
            self.view.button8()
            return

        if "--unidec" in sys.argv[1:] or "-u" in sys.argv[1:]:
            print("Launching UniDec")
            self.view.button1()
            return

        if "--im" in sys.argv[1:] or "-i" in sys.argv[1:]:
            print("Launching UniDec IM")
            self.view.button5()
            return

        if "--ucd" in sys.argv[1:] or "-d" in sys.argv[1:]:
            print("Launching UniDecCD")
            self.view.button9()
            return

        if len(sys.argv) > 1:
            self.view.button8()

    def on_close(self, e=None):
        cancel_event = getattr(self.view, "test_cancel_event", None)
        if cancel_event is not None:
            cancel_event.set()
        self.quit_application()
        self.view.Destroy()
        sys.exit()


class Lview(wx.Frame):
    def __init__(self, parent):
        self.presenter = parent
        self.eng = UniDecEngine()
        wx.Frame.__init__(self, None, title="UniDec Launcher Version " + str(self.eng.version))

        self.test_thread = None
        self.test_cancel_event = None
        self.test_dialog = None
        self.test_timer = None
        self.test_completed = 0
        self.test_message = "Starting GUI smoke tests"

        menu_bar = wx.MenuBar()
        file_menu = wx.Menu()
        self.run_tests_item = file_menu.Append(wx.ID_ANY, "Run Tests")
        file_menu.AppendSeparator()
        about_item = file_menu.Append(wx.ID_ABOUT, "About")
        exit_item = file_menu.Append(wx.ID_EXIT, "Exit")
        menu_bar.Append(file_menu, "File")
        self.SetMenuBar(menu_bar)
        self.Bind(wx.EVT_MENU, self.on_run_tests, self.run_tests_item)
        self.Bind(wx.EVT_MENU, self.on_about, about_item)
        self.Bind(wx.EVT_MENU, self.on_exit, exit_item)

        try:
            if os.path.isfile(self.eng.config.iconfile):
                favicon = wx.Icon(self.eng.config.iconfile, wx.BITMAP_TYPE_ANY)
                wx.Frame.SetIcon(self, favicon)
                self.icon_path = os.path.abspath(self.eng.config.iconfile)
            else:
                self.icon_path = None
        except Exception as e:
            print(e)
            self.icon_path = None

        sizer = wx.GridBagSizer(wx.HORIZONTAL)
        panel = wx.Panel(self)
        button1 = wx.Button(panel, -1, "UniDec\n\nDeconvolve mass spectra")
        button2 = wx.Button(panel, -1, "Data Collector\n\nVisualize multiple spectra\nExtract Trends and Kd's")
        button3 = wx.Button(panel, -1, "Import Wizard\n\nBatch convert Waters Raw to Txt")
        button4 = wx.Button(panel, -1, "MetaUniDec\n\nBatch process and visualize MS spectra")
        button5 = wx.Button(panel, -1, "UniDec IM\n\nDeconvolve ion mobility-mass spectra")
        button6 = wx.Button(panel, -1, "HDF5 Import Wizard\n\nImport Data into HDF5 for MetaUniDec")
        button7 = wx.Button(panel, -1, "UltraMeta Data Collector\n\nVisualize Multiple HDF5 Data Sets\nFit Trends")
        button8 = wx.Button(panel, -1,
                            "UniChrom\n\nCoupled Deconvolution of Chromatograms\nUniDec for LC/MS Data")
        button9 = wx.Button(panel, -1, "UniDecCD\n\nDeconvolution of Charge Detection MS\nUniDec for CD-MS Data")
        button10 = wx.Button(panel, -1, "UniDec Processing Pipeline\n\nBatch Processing Workflow")
        button11 = wx.Button(panel, -1,
                             "UniChromCD\n\nDeconvolution of CD-MS Chromatograms")
        button12 = wx.Button(panel, -1, "IsoDec\n\nDeconvolution of Isotope Patterns")
        html = wx.html.HtmlWindow(panel, -1, size=(390, 310))
        pathtofile = os.path.dirname(os.path.abspath(__file__))
        self.imagepath = self.eng.config.toplogofile
        # print(self.imagepath)
        html.SetPage(
            "<html><body>"
            # "<h1>UniDec</h1>"
            "<img src=\"" + self.imagepath + "\" alt=\"PNG Icon\" height=\"250\" width=\"363\">"
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

        sizer.Add(button10, (4, 0), span=(1, 1), flag=wx.EXPAND)
        sizer.Add(button11, (4, 1), span=(1, 1), flag=wx.EXPAND)
        sizer.Add(button8, (3, 1), flag=wx.EXPAND)
        sizer.Add(button12, (5, 0), span=(1, 1), flag=wx.EXPAND)
        sizer.Add(button5, (5, 1), span=(1, 1), flag=wx.EXPAND)
        sizer.Add(html, (0, 2), span=(6, 2))

        self.Bind(wx.EVT_BUTTON, self.button1, button1)
        self.Bind(wx.EVT_BUTTON, self.button2, button2)
        self.Bind(wx.EVT_BUTTON, self.button3, button3)
        self.Bind(wx.EVT_BUTTON, self.button4, button4)
        self.Bind(wx.EVT_BUTTON, self.button5, button5)
        self.Bind(wx.EVT_BUTTON, self.button6, button6)
        self.Bind(wx.EVT_BUTTON, self.button7, button7)
        self.Bind(wx.EVT_BUTTON, self.button8, button8)
        self.Bind(wx.EVT_BUTTON, self.button9, button9)
        self.Bind(wx.EVT_BUTTON, self.button10, button10)
        self.Bind(wx.EVT_BUTTON, self.button11, button11)
        self.Bind(wx.EVT_BUTTON, self.button12, button12)

        panel.SetSizer(sizer)
        sizer.Fit(self)

        self.Centre()
        self.Show(True)

    def on_about(self, e=None):
        wx.MessageBox(
            "UniDec Launcher\nVersion " + str(self.eng.version),
            "About UniDec",
            wx.OK | wx.ICON_INFORMATION,
            self,
        )

    def on_exit(self, e=None):
        self.presenter.on_close()

    def on_run_tests(self, e=None):
        if self.test_thread is not None and self.test_thread.is_alive():
            return

        from unidec import gui_self_test

        self.run_tests_item.Enable(False)
        self.test_cancel_event = threading.Event()
        self.test_completed = 0
        self.test_message = "Starting GUI smoke tests"
        self.test_dialog = wx.ProgressDialog(
            "UniDec GUI Tests",
            self.test_message,
            maximum=len(gui_self_test.WORKERS),
            parent=self,
            style=wx.PD_CAN_ABORT | wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME,
        )
        self.test_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self._on_test_timer, self.test_timer)
        self.test_timer.Start(250)

        def progress(completed, result):
            wx.CallAfter(self._set_test_progress, completed, result)

        def run():
            try:
                report = gui_self_test.run_controller(
                    callback=progress,
                    cancel_event=self.test_cancel_event,
                )
                wx.CallAfter(self._finish_tests, report, None)
            except Exception as error:
                wx.CallAfter(self._finish_tests, None, str(error))

        self.test_thread = threading.Thread(target=run, name="UniDec GUI self-test", daemon=True)
        self.test_thread.start()

    def _set_test_progress(self, completed, result):
        self.test_completed = completed
        self.test_message = result["label"] + ": " + result["status"]

    def _on_test_timer(self, e=None):
        if self.test_dialog is None:
            return
        keep_going, _ = self.test_dialog.Update(self.test_completed, self.test_message)
        if not keep_going and self.test_cancel_event is not None:
            self.test_cancel_event.set()
            self.test_message = "Cancelling current test"

    def _finish_tests(self, report, error):
        if self.test_timer is not None:
            self.test_timer.Stop()
            self.test_timer.Destroy()
            self.test_timer = None
        if self.test_dialog is not None:
            self.test_dialog.Destroy()
            self.test_dialog = None
        self.run_tests_item.Enable(True)
        self.test_thread = None

        if error:
            message = "Unable to run GUI tests:\n" + error
            style = wx.OK | wx.ICON_ERROR
        else:
            summary = report["summary"]
            message = (
                f"Passed: {summary['pass']}\n"
                f"Skipped: {summary['skip']}\n"
                f"Failed: {summary['failure']}\n"
                f"Timed out: {summary['timeout']}\n"
                f"Cancelled: {summary['cancelled']}\n\n"
                f"Report: {report['report_path']}"
            )
            style = wx.OK | (wx.ICON_INFORMATION if report["success"] else wx.ICON_ERROR)
        wx.MessageBox(message, "UniDec GUI Tests", style, self)

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
        config = unidecstructure.UniDecConfig()
        config.initialize_system_paths()
        print("UniDec Directory: " + config.UniDecDir)
        app = wx.App(False)
        frame = ImportWizard(None, config.UniDecDir)
        frame.Show()
        app.MainLoop()

    def button4(self, e=None):
        print("Launching MetaUniDec")
        app = mudpres.UniDecApp()
        app.start()

    def button5(self, e=None):
        print("Launching UniDec IM")
        app = UniDecIMApp()
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
        print("Launching UniChrom")
        app = ChromApp()
        app.start()

    def button9(self, e=None):
        print("Launching UniDecCD")
        app = UniDecCDApp()
        app.start()

    def button10(self, e=None):
        print("Launching UPP")
        app = UPPApp()

    def button11(self, e=None):
        print("Launching UniChromCD")
        app = UniChromCDApp()
        app.start()

    def button12(self, e=None):
        print("Launching IsoDec")
        app = IsoDecPres()
        app.start()

def run_launcher(*args, **kwargs):
    multiprocessing.freeze_support()
    if "--gui-self-test" in sys.argv or "--gui-self-test-worker" in sys.argv:
        from unidec import gui_self_test
        return gui_self_test.main(sys.argv[1:])
    app = UniDecLauncher(sys.argv[1:])
    app.start()


if __name__ == '__main__':
    raise SystemExit(run_launcher())
