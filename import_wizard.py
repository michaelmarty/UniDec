import os
import threading
import multiprocessing
import wx

from unidec_modules.tims_import_wizard import import_wizard_grid, import_wizard_treectrl, ImporterProgressGauge, \
    data_importer
from unidec_modules.isolated_packages import FileDialogs


class ImportWizard(wx.Frame):
    def __init__(self, parent, dir=None):
        wx.Frame.__init__(self, parent, size=(1000, 900))
        self.setup_frame()
        if dir is not None:
            self.exedir = dir
        else:
            self.exedir = os.path.dirname(os.path.abspath(__file__))

    def setup_frame(self):
        self.CreateStatusBar()
        self.Title = 'Data Conversion Wizard - By Tim Allison with mods from MTM'

        panel = wx.Panel(self)

        menubar = wx.MenuBar()
        self.file_set = wx.Menu()
        self.file_set.Append(wx.ID_CLOSE, 'Close')
        menubar.Append(self.file_set, '&File')
        self.SetMenuBar(menubar)
        self.Bind(wx.EVT_MENU, self.close, id=wx.ID_CLOSE)

        # add box sizers
        hb01 = wx.BoxSizer(wx.HORIZONTAL)
        hb02 = wx.BoxSizer(wx.HORIZONTAL)
        hb04 = wx.BoxSizer(wx.HORIZONTAL)
        hb05 = wx.BoxSizer(wx.HORIZONTAL)

        # Mode which to import
        self.rb = wx.RadioBox(panel, wx.ID_ANY, "Type of Drift Cell",
                              wx.DefaultPosition, wx.DefaultSize,
                              ['Linear', 'T-wave', 'MS Only'], 2)
        defaultselection = 2
        self.rb.SetSelection(defaultselection)
        self.rb.SetToolTip(wx.ToolTip("Select type of drift cell in your ion mobility mass spectrometer"))

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
        self.my_grid = import_wizard_grid.WizardGrid(panel, self)
        self.my_grid.EvtDriftType(n=defaultselection)
        self.my_tree = import_wizard_treectrl.TreeCtrlPanel(panel, self)
        self.tree = self.my_tree.tree

        hb04.Add(self.my_tree)
        hb04.Add(wx.StaticText(panel, wx.ID_ANY, ''))
        # now make comobox for various files
        self.desc = wx.TextCtrl(panel, wx.ID_ANY, '', size=(600, 300), style=wx.TE_MULTILINE)
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
        bottom_btn_sizer.Add(wx.Button(panel, 7, 'Export Grid'), 0, wx.TOP | wx.BOTTOM | wx.RIGHT, border=5)
        bottom_btn_sizer.Add(wx.Button(panel, 8, 'Convert Raw to Txt'), 0, wx.TOP | wx.BOTTOM | wx.RIGHT, border=5)

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
        self.Bind(wx.EVT_BUTTON, self.export_file, id=7)
        self.Bind(wx.EVT_BUTTON, self.export_then_load, id=8)

        self.folder_path.Bind(wx.EVT_KEY_UP, self.on_folder_path_change)

        panel.SetSizerAndFit(sizer)
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

    def add_file(self, evt):
        # scan through selected folders
        for item in self.tree.GetSelections():
            file_path = self.my_tree.tree.GetItemData(item)
            sel = self.rb.GetSelection()
            if sel == 0:
                exp_type = 'linear'
            elif sel == 2:
                exp_type = "ms"
            else:
                exp_type = 't-wave'

            output = data_importer.parse_file(file_path, exp_type=exp_type, dir=self.exedir)
            if output is not None:
                self.my_grid.add_dataset(output)

    def export_then_load(self, evt):
        ''' save file then load it '''

        import_data = []
        if self.rb.GetSelection() == 0:
            import_data.append("ionMS import wizard,Linear,")
        elif self.rb.GetSelection() == 1:
            import_data.append("ionMS import wizard,T-wave,")
        elif self.rb.GetSelection() == 2:
            import_data.append("ionMS import wizard,MS,")
        else:
            print("No valid drift tube selected")

        # this really needs to be factorised

        max_index = self.my_grid.next_free_row()
        if max_index > 0:
            index = -1
            while True:
                cols = self.my_grid.column_labels()
                index += 1
                if index == max_index:
                    break

                if index == 0:
                    # write column labels
                    tmp = ''
                    for c in cols:
                        tmp += '%s,' % cols[c]
                    import_data.append(tmp[:-1])

                # write out user import list
                tmp = ''
                for c in cols:
                    tmp += '%s,' % self.my_grid.GetCellValue(index, c)
                import_data.append(tmp[:-1])

        self.Close()
        self.auto(evt=None, filtered_lines=import_data)
        # self.Close()

    def auto(self, evt, file_path=None, filtered_lines=None):
        if file_path is None and (filtered_lines is None or len(filtered_lines) <= 2):
            file_path = FileDialogs.open_file_dialog("Open File", "*.csv")
            filtered_lines = []

        if file_path is not None:
            # read in the comma-delimited import file
            f = open(file_path, 'r')
            lines = f.readlines()
            f.close()

            # strip lines beginning with # so can comment out lines
            # strip blank lines
            for line in lines:
                if line.strip() is '' or line[0] is '#':
                    continue

                filtered_lines.append(line)

        if filtered_lines is not None and len(filtered_lines) > 2:  # i.e. more than just hopefully the two header lines
            t = threading.Thread(target=data_importer.auto_from_wizard, args=(filtered_lines, self.exedir))
            t.start()
            # gauge=ImporterProgressGauge.ProgressDialog()
            # gauge.progress_dialog(self,"Importing raw file(s)","Importing file %s of %s, please wait..." % (1, len(filtered_lines) - 2),len(filtered_lines) - 2)
            ImporterProgressGauge.progress_dialog(self, "Importing raw file(s)",
                                                  "Importing file %s of %s, please wait..." % (
                                                      1, len(filtered_lines) - 2), len(filtered_lines) - 2)

    def export_file(self, evt):
        '''
        Export import file for later import
        '''
        file_path = FileDialogs.save_file_dialog(message="Save Import CSV File", file_types="CSV (*.csv)|*.csv")
        max_index = self.my_grid.next_free_row()
        if file_path != None and max_index > 0:

            if file_path[-17:] != '_ionMS_import.csv' and file_path[-4:] != '.csv':
                file_path += '_ionMS_import.csv'

            f = open(file_path, 'w')
            # get im type
            if self.rb.GetSelection() == 0:
                f.writelines('ionMS import wizard,Linear,\n')
            elif self.rb.GetSelection() == 2:
                f.writelines("ionMS import wizard,MS,\n")
            else:
                f.writelines('ionMS import wizard,T-wave,\n')

            index = -1
            while True:
                # get column labels
                cols = self.my_grid.column_labels()

                index += 1
                if index == max_index:
                    break

                if index == 0:
                    # write column labels
                    tmp = ''
                    for c in cols:
                        tmp += '%s,' % cols[c]

                    f.writelines(tmp[:-1] + '\n')

                # write out user import list
                tmp = ''
                for c in cols:
                    tmp += '%s,' % self.my_grid.GetCellValue(index, c)

                f.writelines(tmp[:-1] + '\n')

            # don't forget to close the file
            f.close()

            return file_path  # + '_ionMS_import.csv'

    def close(self, evt):
        self.Close()


if __name__ == '__main__':
    multiprocessing.freeze_support()
    app = wx.App(False)
    frame = ImportWizard(None)
    frame.Show()
    app.MainLoop()
