import wx
import wx.grid
from metaunidec.meta_import_wizard import MetaTagTypes as tt


class WizardGrid(wx.grid.Grid):
    '''
    Grid for data import wizard
    '''

    def __init__(self, parent, link, grid_size=(0, 6), labels=0):
        wx.grid.Grid.__init__(self, parent)
        self.CreateGrid(grid_size[0], grid_size[1])

        # events using
        self.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.showPopupMenu)

        # some grid specific dicts
        self.col_header = {0: ['File Name', 'Variable 1', 'Variable 2',
                               'Time Start', 'Time End',
                               'Full Path'],
                           1: ['File Name', 'Variable 1', 'Variable 2',
                               'Scan Start', 'Scan End',
                               'Full Path'],
                           }
        self.col_conv = {'Time': {tt.FILE_NAME: 0,
                                  tt.V1: 1,
                                  tt.V2: 2,
                                  tt.TIME_START: 3,
                                  tt.TIME_END: 4,
                                  tt.FILE_PATH: 5
                                  },
                         'Scans': {tt.FILE_NAME: 0,
                                   tt.V1: 1,
                                   tt.V2: 2,
                                   tt.SCAN_START: 3,
                                   tt.SCAN_END: 4,
                                   tt.FILE_PATH: 5},
                         }

        self.set_labels(labels)
        self.AutoSizeColumns()
        # self.sneaky_resize(600)

    def sneaky_resize(self, panel_width):

        last_col = self.GetNumberCols() - 1
        combined_width = 0
        for col in range(last_col):
            if col != 0:
                combined_width += self.GetColSize(col)

        minSize = panel_width - combined_width - wx.SYS_VSCROLL_X - 27

        # stretch the first column
        self.SetColMinimalWidth(0, minSize)
        self.SetColSize(0, minSize)

        # make the last column with file path invisible
        self.SetColMinimalWidth(last_col, 0)
        self.SetColSize(last_col, 0)

        self.SetMargins(0 - wx.SYS_VSCROLL_X, 0)
        self.ForceRefresh()

    def showPopupMenu(self, evt):
        """
        Create and display a popup menu on right-click event
        """
        # get position of right-click
        self.row, self.col = evt.GetRow(), evt.GetCol()

        if not hasattr(self, "popupID1"):
            self.popupID1 = wx.NewIdRef()
            # self.popupID2 = wx.NewIdRef()
            # self.popupID3 = wx.NewIdRef()
            # make a menu

        menu = wx.Menu()
        # Show how to put an icon in the menu
        item = wx.MenuItem(menu, self.popupID1, "Fill Down")
        menu.AppendItem(item)
        # menu.Append(self.popupID1, "Fill Down")
        # menu.Append(self.popupID1, "Two")
        # menu.Append(self.popupID3, "Three")
        menu.Bind(wx.EVT_MENU, self.fill_down, item)

        # Popup the menu.  If an item is selected then its handler
        # will be called before PopupMenu returns.
        self.PopupMenu(menu)
        menu.Destroy()

    def set_labels(self, mode):
        '''
        Set initial column headers based on mode
            0 - Linear
            1 - T-wave
        '''
        for col, value in enumerate(self.col_header[mode]):
            self.SetColLabelValue(col, value)

        self.SetRowLabelSize(35)
        self.AutoSizeColumns()
        # self.AutoSize()

    def EvtDriftType(self, evt):
        try:
            self.set_labels(evt.GetInt())
        except:
            self.set_labels(evt)
        self.ClearGrid()

    def clear_all(self, evt):
        self.ClearGrid()

    def add_dataset(self, out):
        index = self.next_free_row()
        if index != None:
            if out[tt.TYPE].lower() == 'time':
                conv = self.col_conv['Time']
            elif out[tt.TYPE].lower() == 'scans':
                conv = self.col_conv['Scans']
            for key in conv:
                self.SetCellValue(index, conv[key], str(out[key]))
            self.AutoSizeColumns()

    def next_free_row(self, column=0):
        index = -1;
        found = False
        while not found:
            index += 1
            if self.GetNumberRows() > index:  # check haven't gone past the end
                if self.GetCellValue(index, column) == '':
                    found = True
            if index == self.GetNumberRows():
                self.AppendRows()
                break
        return index

    def column_labels(self):
        col_key = {}
        for col in range(self.GetNumberCols()):
            col_key[col] = self.GetColLabelValue(col)
        return col_key

    def fill_down(self, evt):
        '''
        Try to fill down the columnn
        '''
        # get text to paste
        paste = self.GetCellValue(self.row, self.col)

        row = self.row
        while row < self.GetNumberRows() - 1:
            row += 1
            left = self.GetCellValue(row, 0)
            if left.strip() != '':  # filename col not empty
                current_value = '' if self.GetCellValue(row, self.col) == 'None' else self.GetCellValue(row,
                                                                                                        self.col).strip()
                if current_value == '':
                    self.SetCellValue(row, self.col, paste)
            else:  # filename empty, no more rows to fill for
                break

    def remove_row(self, evt):
        rows = self.GetSelectedRows()
        if rows != []:
            print(rows)

    def get_data(self, e=None):
        import_data = []
        max_index = self.next_free_row()
        if max_index > 0:
            index = -1
            while True:
                cols = self.column_labels()
                index += 1
                if index == max_index:
                    break

                tmp = []
                for c in cols:
                    tmp.append(self.GetCellValue(index, c))
                import_data.append(tmp)
        return cols, import_data

    def load_data(self, array):
        for i, a in enumerate(array):
            index = self.next_free_row()
            if index != None:
                for j, l in enumerate(a[:self.GetNumberCols()]):
                    self.SetCellValue(index, j, str(l))
        self.AutoSizeColumns()
