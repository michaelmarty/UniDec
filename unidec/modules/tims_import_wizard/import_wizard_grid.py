import wx
import wx.grid
from unidec.modules.tims_import_wizard import TagTypes as tt

class WizardGrid(wx.grid.Grid):
    '''
    Grid for data import wizard
    '''

    def __init__(self, parent, link, grid_size=(0, 14), labels=0):
        wx.grid.Grid.__init__(self, parent)
        self.CreateGrid(grid_size[0], grid_size[1])

        # events using
        self.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.showPopupMenu)

        # some grid specific dicts
        self.col_header = { 0  : ['Filename',
                                  #'Sample','Description',
                                  'Pusher', 'Drift Voltage',
                                  'Collision Voltage', 'Drift Pressure',
                                  'Drift Temperature', 'Atom',
                                  'm/z Start', 'm/z End', 'm/z Bin Size', 'Function',
                                  'Scan Start', 'Scan End', 'Full Path'],
                            1 :  ['Filename',
                                  #'Sample','Description',
                                  'Pusher', 'Calibration 1', 'Calibration 2', 'EDC',
                                  'Collision Voltage', 'Atom', 'm/z Start', 'm/z End', 'm/z Bin Size', 'Function',
                                  'Scan Start', 'Scan End', 'Full Path'],
                            2 :  ['Filename',
                                  #'Sample','Description',
                                  'Collision Voltage', 'm/z Start', 'm/z End', 'm/z Bin Size','Function',
                                  'Scan Start', 'Scan End', 'Full Path','','','','',''],
                            3: ['Filename',
                                # 'Sample','Description',
                                'Collision Voltage', 'm/z Start', 'm/z End', 'm/z Bin Size', 'Function',
                                'Time Start', 'Time End', 'Full Path', '', '', '', '', '']
                          }
        self.col_conv = { 'linear' : { tt.FILE_NAME: 0,
                                       #tt.SAMPLE: 1,
                                       #tt.DESCRIPTION: 2,
                                       tt.PUSHER: 1,
                                       tt.DRIFT_V: 2,
                                       tt.COLLISION_V: 3,
                                       tt.DRIFT_PRESSURE: 4,
                                       tt.TEMP: 5,
                                       tt.ATOM:6,
                                       tt.START: 7,
                                       tt.END: 8 ,
                                       tt.BIN: 9 ,
                                       tt.FUNCTION: 10,
                                       tt.SCAN_START: 11,
                                       tt.SCAN_END: 12,
                                       tt.FILE_PATH: 13},
                          't-wave' : { tt.FILE_NAME: 0,
                                       #tt.SAMPLE: 1,
                                       #tt.DESCRIPTION: 2,
                                       tt.PUSHER: 1,
                                       tt.TCAL1: 2,
                                       tt.TCAL2: 3,
                                       tt.EDC: 4,
                                       tt.COLLISION_V: 5,
                                       tt.ATOM: 6,
                                       tt.START: 7,
                                       tt.END: 8,
                                       tt.BIN: 9 ,
                                       tt.FUNCTION: 10,
                                       tt.SCAN_START: 11,
                                       tt.SCAN_END: 12,
                                       tt.FILE_PATH: 13},
                          'ms' : { tt.FILE_NAME: 0,
                                       #tt.SAMPLE: 1,
                                       #tt.DESCRIPTION: 2,
                                       tt.COLLISION_V: 1,
                                       tt.START: 2,
                                       tt.END: 3,
                                       tt.BIN: 4 ,
                                       tt.FUNCTION: 5,
                                       tt.SCAN_START: 6,
                                       tt.SCAN_END: 7,
                                       tt.FILE_PATH: 8},
                          'ms times': {tt.FILE_NAME: 0,
                                 # tt.SAMPLE: 1,
                                 # tt.DESCRIPTION: 2,
                                 tt.COLLISION_V: 1,
                                 tt.START: 2,
                                 tt.END: 3,
                                 tt.BIN: 4,
                                 tt.FUNCTION: 5,
                                 tt.TIME_START: 6,
                                 tt.TIME_END: 7,
                                 tt.FILE_PATH: 8}
                          }

        self.set_labels(labels)
        self.AutoSizeColumns()


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
            #self.popupID2 = wx.NewIdRef()
            #self.popupID3 = wx.NewIdRef()
            # make a menu

        menu = wx.Menu()
        # Show how to put an icon in the menu
        item = wx.MenuItem(menu, self.popupID1, "Fill Down")
        menu.Append(item)
        #menu.Append(self.popupID1, "Fill Down")
        #menu.Append(self.popupID1, "Two")
        #menu.Append(self.popupID3, "Three")
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
        #self.AutoSize()

    def EvtDriftType(self, evt=None, n=None):
        if n is None:
            self.set_labels(evt.GetInt())
        else:
            self.set_labels(n)
        self.ClearGrid()

    def clear_all(self, evt):
        self.ClearGrid()

    def add_dataset(self, out):
        index = self.next_free_row()
        if index != None:
            if out[tt.TYPE] == 'linear':
                conv = self.col_conv['linear']
            elif out[tt.TYPE].lower() == 'ms':
                conv = self.col_conv['ms']
            elif out[tt.TYPE].lower() == 'ms times':
                conv = self.col_conv['ms times']
            else:
                conv = self.col_conv['t-wave']
            for key in conv:
                self.SetCellValue(index, conv[key], str(out[key]))
            self.AutoSizeColumns()

    def next_free_row(self, column=0):
        index = -1; found = False
        while not found:
            index += 1
            if self.GetNumberRows() > index: # check haven't gone past the end
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
        while row < self.GetNumberRows()-1:
            row += 1
            left = self.GetCellValue(row, 0)
            if left.strip() != '': # filename col not empty
                current_value = '' if self.GetCellValue(row, self.col) == 'None' else self.GetCellValue(row, self.col).strip()
                if current_value == '':
                    self.SetCellValue(row, self.col, paste)
            else: # filename empty, no more rows to fill for
                break

    def remove_row(self, evt):
        rows = self.GetSelectedRows()
        if rows != []:
            print(rows)

