import os
import threading
import multiprocessing
import wx
import numpy as np

import meta_import_wizard_grid, meta_import_wizard_treectrl, meta_data_importer
from unidec_modules.isolated_packages import FileDialogs
from import_wizard import ImportWizard as IW


class ImportWizard(wx.Frame):
    def __init__(self, parent, dir=None):
        wx.Frame.__init__(self, parent, size=(1000, 800))
        self.setup_frame()
        if dir is not None:
            self.exedir = dir
        else:
            self.exedir = os.path.dirname(os.path.abspath(__file__))

    def setup_frame(self):
        self.CreateStatusBar()
        self.Title = 'HDF5 File Import Wizard'

        panel = wx.Panel(self)

        menubar = wx.MenuBar()
        self.file_set = wx.Menu()
        self.file_set.Append(wx.ID_OPEN, 'Open CSV')
        self.file_set.Append(wx.ID_SAVE, 'Save CSV')
        self.file_set.AppendSeparator()
        self.iw_menu = self.file_set.Append(wx.ID_ANY, "Waters Conversion Wizard")
        self.file_set.AppendSeparator()
        self.file_set.Append(wx.ID_CLOSE, 'Close')
        menubar.Append(self.file_set, '&File')
        self.SetMenuBar(menubar)
        self.Bind(wx.EVT_MENU, self.close, id=wx.ID_CLOSE)
        self.Bind(wx.EVT_MENU, self.import_file, id=wx.ID_OPEN)
        self.Bind(wx.EVT_MENU, self.export_file, id=wx.ID_SAVE)
        self.Bind(wx.EVT_MENU, self.on_conversion_wizard, self.iw_menu)

        # add box sizers
        hb01 = wx.BoxSizer(wx.HORIZONTAL)
        hb02 = wx.BoxSizer(wx.HORIZONTAL)
        hb04 = wx.BoxSizer(wx.HORIZONTAL)
        hb05 = wx.BoxSizer(wx.HORIZONTAL)

        # Mode which to import
        self.rb = wx.RadioBox(panel, wx.ID_ANY, "Mode",
                              wx.DefaultPosition, wx.DefaultSize,
                              ['Time', 'Scans'], 0)
        self.rb.SetToolTip(wx.ToolTip("How you want to parse data"))

        # folder path stuff
        self.folder_path = wx.TextCtrl(panel, wx.ID_ANY, FileDialogs.default_dir, size=(450, -1))
        hb02.Add(self.folder_path, 0, wx.ALIGN_CENTRE_VERTICAL | wx.LEFT, border=5)
        hb02.Add(wx.Button(panel, 10, 'Browse'), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, border=5)
        box = wx.StaticBoxSizer(wx.StaticBox(panel, wx.ID_ANY, 'Path to Data Folder'), wx.VERTICAL)
        box.Add(hb02)

        hb01.Add(self.rb, 0, wx.ALL, border=5)
        hb01.Add(box, 0, wx.ALL, border=5)

        # add my grid
        panel2 = wx.Panel(panel)
        self.my_grid = meta_import_wizard_grid.WizardGrid(panel, self)
        self.my_tree = meta_import_wizard_treectrl.TreeCtrlPanel(panel, self)
        self.tree = self.my_tree.tree

        hb04.Add(self.my_tree)
        hb04.Add(wx.StaticText(panel, wx.ID_ANY, ''))
        # now make comobox for various files
        self.desc = wx.TextCtrl(panel, wx.ID_ANY, '', size=(400, 300), style=wx.TE_MULTILINE)
        hb04.Add(self.desc, wx.ALL, border=10)

        # make buttons to add, auto, clear all
        hb05.Add(wx.Button(panel, 2, 'Add'), 0, wx.ALL, border=5)
        # hb05.Add(wx.Button(panel, 3, 'Auto'), wx.ALL | wx.ALIGN_LEFT | wx.TOP, 10)
        hb05.Add(wx.Button(panel, 4, 'Clear All'), 0, wx.TOP | wx.BOTTOM | wx.RIGHT, border=5)

        # make a box around the tree ctrl
        box2 = wx.StaticBoxSizer(wx.StaticBox(panel, wx.ID_ANY, 'File(s) to Convert'), wx.VERTICAL)
        box2.Add(hb04)
        box2.Add(hb05)
        box2.Add(self.my_grid, 1, wx.EXPAND | wx.ALL, border=5)

        bottom_btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # bottom_btn_sizer.Add(wx.Button(panel, 7, 'Export Grid'), 0, wx.TOP | wx.BOTTOM | wx.RIGHT, border=5)
        bottom_btn_sizer.Add(wx.Button(panel, 8, 'Load All to HDF5'), 0, wx.TOP | wx.BOTTOM | wx.RIGHT, border=5)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(hb01, 0, wx.ALL, border=0)
        sizer.Add(box2, 1, wx.ALL, border=5)
        # sizer.Add(self.my_grid, 0, wx.EXPAND , border=5)
        # sizer.Add(wx.StaticText(panel, wx.ID_ANY, ''))
        sizer.Add(bottom_btn_sizer, 0, wx.ALIGN_RIGHT | wx.ALL, border=0)

        # bind events
        self.Bind(wx.EVT_RADIOBOX, self.my_grid.EvtDriftType, self.rb)
        self.Bind(wx.EVT_BUTTON, self.get_folder_path, id=10)
        self.Bind(wx.EVT_BUTTON, self.add_file, id=2)
        self.Bind(wx.EVT_BUTTON, self.my_grid.clear_all, id=4)
        # self.Bind(wx.EVT_BUTTON, self.export_file, id=7)
        self.Bind(wx.EVT_BUTTON, self.auto, id=8)

        self.folder_path.Bind(wx.EVT_KEY_UP, self.on_folder_path_change)

        panel.SetSizerAndFit(sizer)
        self.fill_text_box()
        self.SetFocus()
        self.Centre()

    # bound to browse button for folder path
    def get_folder_path(self, evt):
        '''
        Get path to folder and place in txtctrl
        '''
        path = FileDialogs.open_dir_dialog()
        if path != None:
            self.folder_path.SetValue(path)
            self.my_tree.populate_tree()

    # listens to changes in the folder path control
    def on_folder_path_change(self, evt):
        self.my_tree.populate_tree()

    def on_conversion_wizard(self, evt=None):
        print "Launching Waters Converter Wizard"
        app = wx.App(False)
        frame = IW(None)
        frame.Show()
        app.MainLoop()

    def add_file(self, evt):
        # scan through selected folders
        for item in self.tree.GetSelections():
            file_path = self.my_tree.tree.GetItemData(item)
            sel = self.rb.GetSelection()
            if sel == 0:
                exp_type = 'Time'
            elif sel == 1:
                exp_type = "Scans"

            output = meta_data_importer.parse_file(file_path, exp_type=exp_type, dir=self.exedir)
            print output
            if output is not None:
                self.my_grid.add_dataset(output)

    def auto(self, evt=None, file_path=None):
        if file_path is None:
            file_path = FileDialogs.save_file_dialog("Save File As", "*.hdf5")
            filtered_lines = []

        if file_path is not None:
            cols, data = self.my_grid.get_data()
            meta_data_importer.auto_from_wizard(data, file_path, self.rb.GetSelection())

    def import_file(self, evt=None, file_path=None):
        if file_path is None:
            file_path = FileDialogs.open_file_dialog("Open File", "*.csv")

        if file_path is not None:
            # read in the comma-delimited import file
            f = open(file_path, 'r')
            lines = f.readlines()
            f.close()

            type = lines[0].split(",")[1]
            if type.lower() == "scans":
                print "Type Scans:", type
                self.rb.SetSelection(1)
                self.my_grid.EvtDriftType(1)
            else:
                print "Type Time:", type
                self.rb.SetSelection(0)
                self.my_grid.EvtDriftType(0)

            data = np.genfromtxt(file_path, delimiter=",", skip_header=2, dtype=np.str)
            print "Data:", data

            self.my_grid.load_data(data)

    def export_file(self, evt):
        '''
        Export import file for later import
        '''
        file_path = FileDialogs.save_file_dialog(message="Save Import CSV File", file_types="CSV (*.csv)|*.csv")
        max_index = self.my_grid.next_free_row()
        if file_path != None and max_index > 0:

            if file_path[-17:] != '_import.csv' and file_path[-4:] != '.csv':
                file_path += '_import.csv'

            cols, data = self.my_grid.get_data()

            f = open(file_path, 'w')
            # Write Type at the Top
            if self.rb.GetSelection() == 0:
                f.writelines('MetaUniDec Import Wizard,Time,\n')
            elif self.rb.GetSelection() == 1:
                f.writelines("MetaUniDec Import Wizard,Scans,\n")

            # Write Column Values
            tmp = ''
            for c in cols:
                tmp += '%s,' % cols[c]
            f.writelines(tmp[:-1] + '\n')

            # Write Data
            for i, a in enumerate(data):
                for j, l in enumerate(a):
                    f.write(str(l) + ",")
                f.write('\n')

            # close the file
            f.close()

            return file_path

    def fill_text_box(self):
        self.desc.SetValue("First, select the directory. "
                           "Second, select the desired files. Third, you can set the Variable 1 and 2 metadata. "
                           "Fourth, you can select specific time or scans. "
                           "If you leave both as None, it will take all scans/times and sum them into a single spectrum. "
                           "Note: extraction of specific scans/times is only supported for Thermo RAW and mzML files. "
                           "However, opening of full files is supported for text files and Waters RAW files as well. "
                           "If you need to convert Waters files to text, use the File>Waters Conversion Wizard tool. "
                           "\n\nAfter you have entered all of the files and parameters, you can save the table as a CSV file, "
                           "which will allow you to load the same parameters again. You can also generate the CSV in Excel and import it. "
                           "Right clicking a cell will allow you to apply it to all rows. "
                           "\n\nWhen everything is ready, click Load All to HDF5 to convert. "
                           "The HDF5 file can then be opened directly in MetaUniDec."
                           )

    def close(self, evt):
        self.Close()


if __name__ == '__main__':
    multiprocessing.freeze_support()
    app = wx.App(False)
    frame = ImportWizard(None)

    #frame.folder_path.SetValue("C:\\Python\\UniDec\\TestSpectra")
    #frame.my_tree.populate_tree()
    # frame.import_file(file_path="C:\\Python\\UniDec\\TestSpectra\\test2.csv")
    # frame.auto(file_path="C:\\Python\\UniDec\\TestSpectra\\test_out2.hdf5")
    frame.Show()
    app.MainLoop()
