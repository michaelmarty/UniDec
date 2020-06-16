import wx
import unidec_modules.isolated_packages.FileDialogs as FileDialogs
import unidec_modules.isolated_packages.biopolymer_calculator as biopolymer_calculator
import os
import numpy as np
import unidec_modules.unidectools as ud
import unidec_modules.peakwidthtools as peakwidthtools
from unidec_modules import ManualSelectionWindow, AutocorrWindow, miscwindows, peakstructure
import time


class UniDecPres(object):
    """
    Main UniDec GUI Application.

    Presenter contains UniDec engine at self.eng and main GUI window at self.view
    """

    def __init__(self, *args, **kwargs):
        self.wx_app = wx.App(redirect=False)
        self.eng = None
        self.view = None
        self.recent_files = []
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

    # def on_end_session(self):
    #    wx.CallAfter(self.quit_application, force=True)

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
        # tstart = time.perf_counter()
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
                pass
            else:
                self.eng.config.config_export(file_name)
        else:
            if self.eng.config.filetype == 1:
                self.eng.config.write_hdf5()
                pass
        self.eng.update_history()

    def check_badness(self):
        """
        Check for any bad parameters and warn the user if something is off.
        :return: code (see unidecstructure.check_badness)
        """
        badness, warning = self.eng.config.check_badness()
        if warning != "":
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
        try:
            if not wx.GetKeyState(wx.WXK_CONTROL):
                limits = self.view.plot1.subplot1.get_xlim()
                self.view.controls.ctlminmz.SetValue(str(limits[0]))
                self.view.controls.ctlmaxmz.SetValue(str(limits[1]))
                print("New m/z limits:", limits)
                try:
                    self.on_dataprep_button()
                except:
                    pass
        except Exception as e:
            print(e)

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
        print("Saved: ", self.eng.config.defaultconfig)
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
                print("Loaded: ", self.eng.config.defaultconfig)
            except (ValueError, IndexError, TypeError):
                print("Failed to Load: ", self.eng.config.defaultconfig)
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

    def on_auto_peak_width(self, e=None, set=True):
        self.export_config()
        if not ud.isempty(self.eng.data.data2):
            self.eng.get_auto_peak_width(set=set)
            self.import_config()
        else:
            print("Need to process data first")

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
            print("Need to process data first")
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

    def on_biopolymer(self, e=0):
        bp = biopolymer_calculator.BiopolymerFrame(self.view)
        result = bp.Show()
        '''
        if result == 0:
            print("Calculated Mass from Sequence:", bp.mass)
            calc_mass = 0
            try:
                calc_mass = float(bp.mass)
                if calc_mass > 0:
                    label = bp.mass
                    mval = np.amax(self.eng.data.massdat[:, 1])
                    self.view.plot2.addtext(label, calc_mass, mval * 0.96, vlines=False)
                    pass
            except:
                pass
        bp.Destroy()'''

    def plot_theo_mass(self, e=0):
        dialog = miscwindows.SingleInputDialog(self.view)
        dialog.initialize_interface(title="Theoretical Mass", message="Set Theoretical Mass to Plot (Da): ",
                                    defaultvalue="")
        dialog.ShowModal()

        try:
            mass = float(dialog.value)
        except:
            print("Error with theoretical mass input:", dialog.value)
            mass = 0

        if self.eng.config.massbins >= 1:
            label = str(int(round(mass)))
        else:
            label = str(mass)
        mval = np.amax(self.eng.data.massdat[:, 1])
        self.view.plot2.addtext(label, mass, mval * 0.96, vlines=True)
        self.on_charge_states(mass=mass)

    def register(self, e=None):
        from unidec_modules.data_reader import register
        register()

    def on_linreg(self, e=0):
        fit, rsquared = self.eng.linear_regression_peaks()
        message = "Slope: " + str(fit[0]) + "\nIntercept: " + str(fit[1]) + "\nR-Squared: " + str(rsquared)
        self.warn(message, caption="Linear Regression Results")

    def fpop(self, e=0):
        output = self.eng.oxidation_analysis()
        outstring = ""
        for o in output:
            outstring += str(o) + "\n"
        # Create text data object
        clipboard = wx.TextDataObject()

        # Set data object value
        clipboard.SetText(outstring)

        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()

    def sub_div(self, e=0):
        try:
            masses = []
            for p in self.eng.pks.peaks:
                if p.ignore == 0:
                    masses.append(p.mass)
            defaultsub = np.amin(masses)
        except:
            defaultsub = 44088

        defaultdiv = self.eng.config.molig
        dlg = miscwindows.DoubleInputDialog(self.view)
        dlg.initialize_interface("Subtract and Divide", "Subtract:", str(defaultsub),
                                 "Divide:", str(defaultdiv))
        dlg.ShowModal()
        try:
            sub = float(dlg.value)
            div = float(dlg.value2)
        except:
            print("Error with Subtract and Divide Inputs:", dlg.value, dlg.value2)
            return 0

        try:
            sd_results = []
            message = "Average Masses:\n\n"
            outstring = ""
            for s in self.eng.data.spectra:
                pks = peakstructure.Peaks()
                peaks = ud.peakdetect(s.massdat, self.eng.config)
                pks.add_peaks(peaks, massbins=self.eng.config.massbins)
                sd_result = ud.subtract_and_divide(pks, sub, div)
                sd_results.append(sd_result)

                outstring += str(sd_result) + "\n"

            avg = np.mean(sd_results)
            message += outstring
            message += "\nOverall Average: " + str(avg)
            outstring += "\n" + str(avg)
            self.copy_to_clipboard(outstring)
            self.warn(message, caption="Subtract and Divide Results")

        except:
            sd_result = ud.subtract_and_divide(self.eng.pks, sub, div)
            outstring = str(sd_result)
            message = "Average Mass: " + outstring
            self.copy_to_clipboard(outstring)
            self.warn(message, caption="Subtract and Divide Results")

    def copy_to_clipboard(self, outstring):
        # Create text data object
        clipboard = wx.TextDataObject()

        # Set data object value
        clipboard.SetText(outstring)

        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()

    def write_to_recent(self, path=None):
        if path is None:
            path = os.path.join(self.eng.config.dirname, self.eng.config.filename)

        with open(self.eng.config.recentfile, 'a') as file:
            file.write(path + "\n")

        self.recent_files = self.read_recent()
        self.cleanup_recent_file(self.recent_files)

    def read_recent(self):
        if os.path.isfile(self.eng.config.recentfile):
            with open(self.eng.config.recentfile, 'r') as file:
                lines = []
                for l in file:
                    p = l.strip("\n")
                    if os.path.isfile(p):
                        lines.append(p)
        else:
            self.eng.config.recentfile = os.path.join(self.eng.config.UniDecDir, "recent.txt")
            with open(self.eng.config.recentfile, 'w') as file:
                pass  # Create empty file
            return []

        if not ud.isempty(lines):
            lines = np.array(lines)
            lines = lines[::-1][:10]
            unique, indexes = np.unique(lines, return_index=True)
            lines = np.array(unique[indexes.argsort()])
        else:
            lines = []
        return lines

    def cleanup_recent_file(self, recent_files):
        with open(self.eng.config.recentfile, "r+") as f:
            for p in recent_files[::-1]:
                l = p + '\n'
                f.write(l)
            f.truncate()
