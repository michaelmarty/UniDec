import wx
import unidec_modules.isolated_packages.FileDialogs as FileDialogs
import os
import numpy as np
import unidec_modules.unidectools as ud
import unidec_modules.peakwidthtools as peakwidthtools
from unidec_modules import ManualSelectionWindow, AutocorrWindow


class UniDecPres(object):
    """
    Main UniDec GUI Application.

    Presenter contains UniDec engine at self.eng and main GUI window at self.view
    """

    def __init__(self, *args, **kwargs):
        self.wx_app = wx.App(redirect=False)
        self.eng = None
        self.view = None
        pass

    def start(self):
        """
        Launch view and start MainLoop
        :return:
        """
        self.wx_app.SetTopWindow(self.view)
        self.wx_app.SetAppName("GUniDec")
        self.wx_app.SetVendorName("Michael Marty - University of Arizona")

        self.view.Show()
        self.wx_app.MainLoop()

    def on_end_session(self):
        wx.CallAfter(self.quit_application, force=True)

    def quit_application(self):
        # self.wx_app.ProcessIdle()
        self.wx_app.ExitMainLoop()
        return True

    def import_config(self, file_name=None):
        """
        Import configuration from file to engine (if file_name is specified) and from engine to GUI.
        :param file_name: Path of file to import
        :return: None
        """
        if file_name is not None:
            extension = os.path.splitext(file_name)[1]
            if extension == ".hdf5":
                self.eng.config.read_hdf5(file_name)
            else:
                self.eng.config.config_import(file_name)
        self.view.import_config_to_gui()
        if self.eng.config.filetype == 1:
            self.eng.config.write_hdf5()
        self.eng.update_history()

    def export_config(self, file_name=None):
        """
        Get configuration from GUI and (if file_name is specified) write from engine to file_name
        :param file_name: Path of file to save config to
        :return: None
        """
        self.view.export_gui_to_config()
        if file_name is not None:
            extension = os.path.splitext(file_name)[1]
            if extension == ".hdf5":
                self.eng.config.write_hdf5(file_name)
            else:
                self.eng.config.config_export(file_name)
        else:
            if self.eng.config.filetype == 1:
                self.eng.config.write_hdf5()
        self.eng.update_history()

    def check_badness(self):
        """
        Check for any bad parameters and warn the user if something is off.
        :return: code (see unidecstructure.check_badness)
        """
        badness, warning = self.eng.config.check_badness()
        if warning is not "":
            self.warn(warning)
        return badness

    def warn(self, message, caption='Warning!'):
        """
        Send the user a message box.
        :param message: Message string
        :param caption: Caption string
        :return: None
        """
        dlg = wx.MessageDialog(self.view, message, caption, wx.OK | wx.ICON_WARNING)
        dlg.ShowModal()
        dlg.Destroy()

    def on_get_mzlimits(self):
        limits = self.view.plot1.subplot1.get_xlim()
        self.view.controls.ctlminmz.SetValue(str(limits[0]))
        self.view.controls.ctlmaxmz.SetValue(str(limits[1]))
        print "New m/z limits:", limits

    def on_load_conf_file(self, e=None):
        """
        Opens a file dialog and then imports a new _conf.dat config file
        :param e: unused space for event
        :return: None
        """
        cfilename = FileDialogs.open_file_dialog("Open Configuration File (_conf.dat)", file_types="*.*")
        if cfilename is not None:
            self.import_config(cfilename)
        pass

    def on_save_default(self, e=None):
        """
        Saves the default config file to self.eng.config.defaultconfig
        :param e: unused space for event
        :return: None
        """
        print "Saved: ", self.eng.config.defaultconfig
        self.export_config(self.eng.config.defaultconfig)
        pass

    def on_load_default(self, e=None):
        """
        Resets the configuration then loads the default config file from self.eng.config.defaultconfig.
        Ignores min and max data cutoffs.
        :param e: unused space for event
        :return: None
        """
        self.on_reset(e)
        if os.path.isfile(self.eng.config.defaultconfig):
            try:
                self.import_config(self.eng.config.defaultconfig)
                try:
                    if not ud.isempty(self.eng.data.rawdata):
                        self.view.controls.ctlminmz.SetValue(str(np.amin(self.eng.data.rawdata[:, 0])))
                        self.view.controls.ctlmaxmz.SetValue(str(np.amax(self.eng.data.rawdata[:, 0])))
                        if self.eng.config.imflag == 1:
                            self.view.controls.ctlmindt.SetValue(str(np.amin(self.eng.data.rawdata3[:, 1])))
                            self.view.controls.ctlmaxdt.SetValue(str(np.amax(self.eng.data.rawdata3[:, 1])))
                except:
                    pass
                print "Loaded: ", self.eng.config.defaultconfig
            except (ValueError, IndexError, TypeError):
                print "Failed to Load: ", self.eng.config.defaultconfig
        self.view.SetStatusText("Loaded Default", number=5)
        pass

    def on_reset(self, e=None):
        """
        Resets UniDecConfig to default and then loads to self.view.
        :param e: unused space for event
        :return:
        """
        self.eng.reset_config()
        self.import_config(None)
        self.view.SetStatusText("Reset", number=5)

    def on_auto_peak_width(self, e=None):
        self.export_config()
        if not ud.isempty(self.eng.data.data2):
            self.eng.get_auto_peak_width()
            self.import_config()
        else:
            print "Need to process data first"

    def on_peak_width_tool(self, e=None):
        """
        Open peak width tool window. After it has returned, update the GUI to reflect the new peak widths.
        :param e: unused event
        :return: None
        """
        self.export_config()
        if not ud.isempty(self.eng.data.data2):
            self.export_config(None)
            if self.eng.config.imflag == 0:
                dlg = peakwidthtools.PeakTools1d(self.view)
                dlg.initialize_interface(self.eng.config, self.eng.data.data2)
            else:
                dlg = peakwidthtools.PeakTools2d(self.view)
                dlg.initialize_interface(self.eng.data.data3, self.eng.data.data2, self.eng.config)
            dlg.ShowModal()
            self.import_config(None)
        else:
            print "Need to process data first"
        pass

    def on_manual(self, e=None):
        """
        Opens window for setting manual assignments. Window directly modifies self.eng.config.
        :param e: unused event
        :return: None
        """
        dlg = ManualSelectionWindow.ManualSelection(self.view)
        if self.eng.config.imflag == 0:
            dlg.initiate_dialog(self.eng.config, self.eng.data.data2)
        else:
            dlg.initiate_dialog(self.eng.config, self.eng.data.data3)
        dlg.ShowModal()

    def on_match(self, e=None):
        """
        Automatic matching to present oligomer list.

        Opens masstools window but doesn't make it visible.
        Matches results to any combination of oligomers.
        Hits ok and plots results.

        :param e:
        :return:
        """
        # TODO: Rewrite this so that it doesn't need to fake open the window
        self.on_mass_tools(0, show=False)

    def on_autocorr_window(self, e=None):
        """
        Opens the autocorrelation window. Feed the config and massdat from the engine.
        :param e: Unused event
        :return: None
        """
        dlg = AutocorrWindow.AutocorrWindow(self.view)
        dlg.initalize_dialog(self.eng.config, self.eng.data.massdat)
        dlg.ShowModal()