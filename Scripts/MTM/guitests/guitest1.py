import wx
import wx.lib.scrolledpanel as scrolled
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.figure import Figure

class PlottingWindow2(wx.Window):
    def __init__(self, *args, **kwargs):
        wx.Window.__init__(self, *args, **kwargs)
        self.figure = Figure(figsize=(6, 5))
        self.canvas = FigureCanvasWxAgg(self, -1, self.figure)
        self.Bind(wx.EVT_SIZE, self.size_handler)
        self.Bind(wx.EVT_MOUSEWHEEL, self.on_wheel)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)

    def size_handler(self, *args, **kwargs):
        self.canvas.SetSize(self.GetSize())

    def on_scroll(self, event):
        print("MPL Scroll")
        event.guiEvent.Skip()

    def on_wheel(self, event):
        print("Canvas Scroll")
        event.Skip()

class MyFrame(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=(400, 400))
        plotwindow = scrolled.ScrolledPanel(self, -1)
        sizerplot = wx.BoxSizer(wx.VERTICAL)
        self.plot1 = PlottingWindow2(plotwindow)
        sizerplot.Add(self.plot1, 0, flag=wx.EXPAND)
        test1 = wx.TextCtrl(plotwindow, size=(400, 400))
        sizerplot.Add(test1, 0, flag=wx.EXPAND)
        plotwindow.SetSizer(sizerplot)
        plotwindow.SetupScrolling()
        plotwindow.Bind(wx.EVT_MOUSEWHEEL, self.onWheel)
        self.Show(True)

    def onWheel(self, event):
        print("WX Scroll")
        event.Skip()

app = wx.App(False)
frame = MyFrame(None, "Test")
app.MainLoop()
