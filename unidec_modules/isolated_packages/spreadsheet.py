import wx.grid as gridlib
import wx


class Spreadsheet(wx.Frame):
    """"""

    def __init__(self, nrows, ncolumns, title=""):
        """Constructor"""
        wx.Frame.__init__(self, parent=None, title=title)
        panel = wx.Panel(self)

        self.ss = gridlib.Grid(panel)
        self.nrows = nrows
        self.ncolumns = ncolumns
        self.ss.CreateGrid(nrows, ncolumns)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.ss, 1, wx.EXPAND)
        panel.SetSizer(sizer)
        self.Show()


if __name__ == "__main__":
    app = wx.App()
    frame = Spreadsheet(12, 8)
    app.MainLoop()
