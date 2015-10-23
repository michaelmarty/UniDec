import sys
import string
from copy import deepcopy

import wx.lib.mixins.listctrl  as  listmix
import wx
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
from scipy.spatial.distance import euclidean
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve

from unidec_modules import plot1d, plot2d, peakstructure
from unidec_modules.isolated_packages import FileDialogs
import unidec_modules.unidectools as ud


class TestListCtrl(wx.ListCtrl,
                   listmix.ListCtrlAutoWidthMixin,
                   listmix.TextEditMixin):
    def __init__(self, parent, ID, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "Mass (Da)")
        self.SetColumnWidth(0,width=190)#, wx.LIST_AUTOSIZE)

    def Populate(self, listctrldata):
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, len(listctrldata)):
            index = self.InsertStringItem(sys.maxint, str(listctrldata[i]))
            self.SetItemData(index, i)
        self.currentItem = 0

    def Clear(self):
        self.DeleteAllItems()

    def AddLine(self):
        self.InsertStringItem(sys.maxint, str(0))

    def GetList(self):
        count = self.GetItemCount()
        list = []
        for i in range(0, count):
            list.append(float(self.GetItemText(i)))
        return list


class TestListCtrlPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)

        tID = wx.NewId()

        sizer = wx.BoxSizer(wx.VERTICAL)

        if wx.Platform == "__WXMAC__" and \
                hasattr(wx.GetApp().GetTopWindow(), "LoadDemo"):
            self.useNative = wx.CheckBox(self, -1, "Use native listctrl")
            self.useNative.SetValue(
                not wx.SystemOptions.GetOptionInt("mac.listctrl.always_use_generic"))
            self.Bind(wx.EVT_CHECKBOX, self.OnUseNative, self.useNative)
            sizer.Add(self.useNative, 0, wx.ALL | wx.ALIGN_RIGHT, 4)

        self.list = TestListCtrl(self, tID, size=(210, 380),style=wx.LC_REPORT )

        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)

    def OnUseNative(self, event):
        wx.SystemOptions.SetOptionInt("mac.listctrl.always_use_generic", not event.IsChecked())
        wx.GetApp().GetTopWindow().LoadDemo("ListCtrl_edit")


class TestListCtrl2(wx.ListCtrl,
                    listmix.ListCtrlAutoWidthMixin,
                    listmix.TextEditMixin):
    def __init__(self, parent, ID, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "Base Offset (Da)")
        self.InsertColumn(1, "Monomer Mass (Da)")
        self.InsertColumn(2, "Min # of Oligomers")
        self.InsertColumn(3, "Max # of Oligomers")
        self.InsertColumn(4, "Name")
        self.SetColumnWidth(0, 100)
        self.SetColumnWidth(1, 125)
        self.SetColumnWidth(2, 100)
        self.SetColumnWidth(3, 100)
        self.SetColumnWidth(4, 60)
        self.index=0

    def Clear(self):
        self.DeleteAllItems()
        self.index=0

    def AddLine(self):
        index=self.InsertStringItem(sys.maxint, str(0))
        self.SetStringItem(index, 1, str(0))
        self.SetStringItem(index, 2, str(0))
        self.SetStringItem(index, 3, str(1))
        self.SetStringItem(index, 4, string.uppercase[self.index])
        self.index=self.index+1


    def Populate(self, data):
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, len(data)):
            index=self.InsertStringItem(sys.maxint, str(data[i][0]))
            try:
                self.SetStringItem(index, 1, str(data[i][1]))
                self.SetStringItem(index, 2, str(data[i][2]))
                self.SetStringItem(index, 3, str(data[i][3]))
            except:
                self.SetStringItem(index, 1, "")
                self.SetStringItem(index, 2, "")
                self.SetStringItem(index, 3, "")
            try:
                self.SetStringItem(index, 4, str(data[i][4]))
            except:
                self.SetStringItem(index, 4, "")
            #self.SetItemData(index, i)


    def GetList(self):
        count = self.GetItemCount()
        list = []
        for i in range(0, count):
            sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),int(self.GetItem(i, col=2).GetText()), int(self.GetItem(i, col=3).GetText()),self.GetItem(i, col=4).GetText()]
            list.append(sublist)
        return list


class TestListCtrlPanel2(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)

        tID = wx.NewId()

        sizer = wx.BoxSizer(wx.VERTICAL)

        if wx.Platform == "__WXMAC__" and \
                hasattr(wx.GetApp().GetTopWindow(), "LoadDemo"):
            self.useNative = wx.CheckBox(self, -1, "Use native listctrl")
            self.useNative.SetValue(
                not wx.SystemOptions.GetOptionInt("mac.listctrl.always_use_generic"))
            self.Bind(wx.EVT_CHECKBOX, self.OnUseNative, self.useNative)
            sizer.Add(self.useNative, 0, wx.ALL | wx.ALIGN_RIGHT, 4)

        self.list = TestListCtrl2(self, tID, size=(500, 200),style=wx.LC_REPORT| wx.LC_SORT_ASCENDING)

        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


    def OnUseNative(self, event):
        wx.SystemOptions.SetOptionInt("mac.listctrl.always_use_generic", not event.IsChecked())
        wx.GetApp().GetTopWindow().LoadDemo("ListCtrl_edit")

class TestListCtrlPanelMatch(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)

        tID = wx.NewId()

        sizer = wx.BoxSizer(wx.VERTICAL)

        if wx.Platform == "__WXMAC__" and \
                hasattr(wx.GetApp().GetTopWindow(), "LoadDemo"):
            self.useNative = wx.CheckBox(self, -1, "Use native listctrl")
            self.useNative.SetValue(
                not wx.SystemOptions.GetOptionInt("mac.listctrl.always_use_generic"))
            self.Bind(wx.EVT_CHECKBOX, self.OnUseNative, self.useNative)
            sizer.Add(self.useNative, 0, wx.ALL | wx.ALIGN_RIGHT, 4)

        self.list = TestListCtrlMatch(self, tID, size=(500, 200),style=wx.LC_REPORT| wx.LC_SORT_ASCENDING)

        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)


    def OnUseNative(self, event):
        wx.SystemOptions.SetOptionInt("mac.listctrl.always_use_generic", not event.IsChecked())
        wx.GetApp().GetTopWindow().LoadDemo("ListCtrl_edit")


class TestListCtrlMatch(wx.ListCtrl,
                    listmix.ListCtrlAutoWidthMixin,
                    listmix.TextEditMixin):
    def __init__(self, parent, ID, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(0, "Peak Mass (Da)")
        self.InsertColumn(1, "Match")
        self.InsertColumn(2, "Error")
        self.InsertColumn(3, "Name")
        self.SetColumnWidth(0, 100)
        self.SetColumnWidth(1, 100)
        self.SetColumnWidth(2, 100)
        self.SetColumnWidth(3, 100)
        self.index=0

    def Clear(self):
        self.DeleteAllItems()
        self.index=0

    def Populate(self, data1,data2,data3,data4):
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, len(data1)):
            index=self.InsertStringItem(sys.maxint, str(data1[i]))
            self.SetStringItem(index, 1, str(data2[i]))
            self.SetStringItem(index, 2, str(data3[i]))
            self.SetStringItem(index, 3, str(data4[i]))
            #self.SetItemData(index, i)


    def GetList(self):
        count = self.GetItemCount()
        list = []
        for i in range(0, count):
            sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),float(self.GetItem(i, col=2).GetText()), self.GetItem(i, col=3).GetText()]
            list.append(sublist)
        return list


class CorrListCtrl(wx.ListCtrl,
                   listmix.ListCtrlAutoWidthMixin,
                   listmix.TextEditMixin):
    def __init__(self, parent, ID, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)
        self.InsertColumn(1, "Mass Diff (Da)")
        self.SetColumnWidth(1,width=100)#, wx.LIST_AUTOSIZE)
        self.InsertColumn(0, "")
        self.SetColumnWidth(0,width=25)#, wx.LIST_AUTOSIZE)

    def Populate(self, pks):
        self.DeleteAllItems()
        # for normal, simple columns, you can add them like this:
        for i in range(0, pks.plen):
            p=pks.peaks[i]
            index = self.InsertStringItem(i, p.textmarker)
            self.SetStringItem(i, 1, str(p.mass))
            color=wx.Colour(round(p.color[0]*255),round(p.color[1]*255),round(p.color[2]*255),alpha=255)
            self.SetItemBackgroundColour(i,col=color)

    def Clear(self):
        self.DeleteAllItems()

class CorrListCtrlPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)

        tID = wx.NewId()

        sizer = wx.BoxSizer(wx.VERTICAL)

        if wx.Platform == "__WXMAC__" and \
                hasattr(wx.GetApp().GetTopWindow(), "LoadDemo"):
            self.useNative = wx.CheckBox(self, -1, "Use native listctrl")
            self.useNative.SetValue(
                not wx.SystemOptions.GetOptionInt("mac.listctrl.always_use_generic"))
            self.Bind(wx.EVT_CHECKBOX, self.OnUseNative, self.useNative)
            sizer.Add(self.useNative, 0, wx.ALL | wx.ALIGN_RIGHT, 4)

        self.list = CorrListCtrl(self, tID, size=(200, 550),style=wx.LC_REPORT )

        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)

    def OnUseNative(self, event):
        wx.SystemOptions.SetOptionInt("mac.listctrl.always_use_generic", not event.IsChecked())
        wx.GetApp().GetTopWindow().LoadDemo("ListCtrl_edit")

class AutocorrWindow(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER,*args, **kwargs)
        self.SetSize((800, 600))
        self.SetTitle("Autocorrelation Plots")

    def InitUI(self, config,massdat):
        self.config=config
        self.massdat=massdat

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.sb = wx.StaticBox(self.pnl, label='Autocorrelation')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.plot1= plot1d.Plot1d(self.pnl)
        self.sbs.Add(self.plot1)
        self.hbox.Add(self.sbs)
        self.listpanel=CorrListCtrlPanel(self.pnl)
        self.hbox.Add(self.listpanel)
        self.pnl.SetSizer(self.hbox)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        hboxend.Add(okButton)
        hboxend.Add(closeButton, flag=wx.LEFT, border=5)
        self.vbox.Add(self.pnl, proportion=1,flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(hboxend,flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)
        self.SetSizer(self.vbox)
        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)
        self.CenterOnParent()
        self.Go(0)

    def Go(self,e):
        #corr=np.correlate(self.massdat[:,1],self.massdat[:,1],mode="same")
        corr=fftconvolve(self.massdat[:,1],self.massdat[:,1][::-1],mode='same')
        corr=corr/np.amax(corr)
        xdiff=self.massdat[1,0]-self.massdat[0,0]
        corrx=np.arange(0.0,len(corr))*xdiff
        maxpos=np.argmax(corr)
        corrx=corrx-corrx[maxpos]
        self.corr=corr
        self.corrx=corrx
        self.plot1.plotrefreshtop(corrx,corr,"Autocorrelation","Mass Difference","","",self.config)
        self.PickPeaks(0)

    def PickPeaks(self,e):
        self.pks2= peakstructure.Peaks()
        self.peaks=ud.peakdetect(np.transpose([self.corrx,self.corr]),self.config)
        self.pks2.add_peaks(self.peaks)
        self.pks2.default_params()
        if self.pks2.plen>0:
            for p in self.pks2.peaks:
                if p.ignore==0:
                    self.plot1.plotadddot(p.mass, p.height, p.color, p.marker)
            self.plot1.repaint()
        print self.peaks
        self.listpanel.list.Populate(self.pks2)

    def OnClose(self,e):
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self,e):
        self.Destroy()
        self.EndModal(0)

class MassSelection(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER,*args, **kwargs)
        self.SetSize((800, 700))
        self.SetTitle("Mass and Oligomer Tools")

    def InitUI(self, config,matchlistin,pks,massdat=None):
        massbins=config.massbins
        self.massdat=massdat
        self.config=config
        self.defaultmasslist = deepcopy(self.config.masslist)
        self.defaultoligolist = deepcopy(self.config.oligomerlist)
        self.defaultmatchlist=matchlistin
        self.newmasslist = deepcopy(self.defaultmasslist)
        self.oligos = deepcopy(self.defaultoligolist)
        self.matchlist=deepcopy(self.defaultmatchlist)
        self.newmatchlist=deepcopy(self.defaultmatchlist)
        self.pks=pks

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.sb = wx.StaticBox(self.pnl, label='Set the Mass List for Limited UniDec')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)
        self.peakpop = wx.Button(self.pnl, label="Populate from Peak List")
        self.Bind(wx.EVT_BUTTON, self.OnPopulateButton, self.peakpop)

        self.importbutton = wx.Button(self.pnl, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.OnImport, self.importbutton)

        self.oligopop = wx.Button(self.pnl, label="Populate from Isolated Oligomers")
        self.Bind(wx.EVT_BUTTON, self.OnPopulateButton2, self.oligopop)

        self.oligopop2 = wx.Button(self.pnl, label="Populate from All Possible Oligomers")
        self.Bind(wx.EVT_BUTTON, self.OnPopulateButton3, self.oligopop2)

        self.clearbutt = wx.Button(self.pnl, label="Clear List")
        self.Bind(wx.EVT_BUTTON, self.OnClear, self.clearbutt)

        self.addbutton = wx.Button(self.pnl, label="Manual Add Species")
        self.Bind(wx.EVT_BUTTON, self.OnAdd, self.addbutton)

        self.simbutton = wx.Button(self.pnl, label="Simulate These Masses")
        self.Bind(wx.EVT_BUTTON, self.OnSim, self.simbutton)

        self.sbs.Add(self.peakpop,0,wx.EXPAND)
        self.sbs.Add(self.importbutton,0,wx.EXPAND)
        self.sbs.Add(self.oligopop,0,wx.EXPAND)
        self.sbs.Add(self.oligopop2,0,wx.EXPAND)
        self.sbs.Add(self.addbutton,0,wx.EXPAND)
        self.sbs.Add(self.clearbutt,0,wx.EXPAND)
        self.masslistbox = TestListCtrlPanel(self.pnl)

        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick, self.masslistbox.list)
        self.sbs.Add(wx.StaticText(self.pnl, label="Mass List"))
        self.sbs.Add(self.masslistbox)
        self.sbs.Add(self.simbutton,0,wx.EXPAND)

        self.hbox.Add(self.sbs)

        self.sb2 = wx.StaticBox(self.pnl, label='Oligomer Maker')
        self.sbs2 = wx.StaticBoxSizer(self.sb2, orient=wx.VERTICAL)

        self.clearbutt2 = wx.Button(self.pnl, label="Clear Oligomer List")
        self.Bind(wx.EVT_BUTTON, self.OnClear2, self.clearbutt2)

        self.addbutton2 = wx.Button(self.pnl, label="Add Oligomer Species")
        self.Bind(wx.EVT_BUTTON, self.OnAdd2, self.addbutton2)

        self.importbutton2 = wx.Button(self.pnl, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.OnImport2, self.importbutton2)

        self.plotbutton = wx.Button(self.pnl, label="View Autocorrelation Plot")
        self.Bind(wx.EVT_BUTTON, self.OnPlot, self.plotbutton)
        self.buttonbox=wx.BoxSizer(wx.VERTICAL)
        self.hbox3=wx.BoxSizer(wx.HORIZONTAL)
        self.buttonbox.Add(self.importbutton2,0,wx.EXPAND)
        self.buttonbox.Add(self.addbutton2,0,wx.EXPAND)
        self.buttonbox.Add(self.clearbutt2,0,wx.EXPAND)
        self.buttonbox.Add(self.plotbutton,0,wx.EXPAND)
        self.hbox3.Add(self.buttonbox)
        text=wx.StaticText(self.pnl, label="  For i from Min # to Max #:\n      Mass(i)=Base Offset + Monomer Mass * i ")
        font=wx.Font(12,wx.FONTFAMILY_DEFAULT,wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_BOLD,False)
        text.SetFont(font)
        self.hbox3.Add(text)
        self.sbs2.Add(self.hbox3,0,wx.EXPAND)

        self.oligomerlistbox = TestListCtrlPanel2(self.pnl)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick2, self.oligomerlistbox.list)
        self.sbs2.Add(wx.StaticText(self.pnl, label="Oligomer List"))
        self.sbs2.Add(self.oligomerlistbox)
        self.sbs2.Add(wx.StaticText(self.pnl, label=""))

        self.sb4=wx.StaticBox(self.pnl,label="Match Peaks to Oligomers")
        self.sbs4 = wx.StaticBoxSizer(self.sb4, orient=wx.VERTICAL)
        self.matchIbutt = wx.Button(self.pnl, label="Match to Isolated Oligomers")
        self.matchIbutt.SetToolTip(wx.ToolTip("Match peaks to isolated oligomers from Oligomer Maker."))
        self.Bind(wx.EVT_BUTTON, self.OnMatchI, self.matchIbutt)
        self.matchAllbutt = wx.Button(self.pnl, label="Matched to Mixed Oligomers")
        self.matchAllbutt.SetToolTip(wx.ToolTip("Match peaks to any possible combination of oligomers from Oligomer Maker."))
        self.Bind(wx.EVT_BUTTON, self.OnMatchAll, self.matchAllbutt)
        self.matchlistbox=TestListCtrlPanelMatch(self.pnl)
        self.hbox2=wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.matchIbutt,1,wx.EXPAND)
        self.hbox2.Add(self.matchAllbutt,1,wx.EXPAND)
        self.sbs4.Add(self.hbox2,0,wx.EXPAND)
        self.sbs4.Add(self.matchlistbox)
        self.sbs2.AddStretchSpacer(prop=1)
        self.sbs2.Add(self.sbs4)
        self.hbox.Add(self.sbs2)

        self.pnl.SetSizer(self.hbox)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        hboxend.Add(okButton)
        hboxend.Add(closeButton, flag=wx.LEFT, border=5)

        self.vbox.Add(self.pnl, proportion=1,
                      flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(hboxend,
                      flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(self.vbox)

        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)

        self.masslistbox.list.Populate(self.defaultmasslist)
        self.oligomerlistbox.list.Populate(self.defaultoligolist)
        if(len(self.defaultmatchlist)==4):
            self.matchlistbox.list.Populate(self.defaultmatchlist[0],self.defaultmatchlist[1],self.defaultmatchlist[2],self.defaultmatchlist[3])

        self.CenterOnParent()

    def OnSim(self,e):
        self.newmasslist = self.masslistbox.list.GetList()
        self.newpeaks=[]
        f=interp1d(self.massdat[:,0],self.massdat[:,1],bounds_error=False,fill_value=0)
        for mass in self.newmasslist:
            intensity=f(mass)
            self.newpeaks.append([mass,intensity])
        self.newpeaks=np.array(self.newpeaks)
        self.pks.Peaks()
        self.pks.add_peaks(self.newpeaks)
        self.pks.default_params(cmap=self.config.peakcmap)
        self.pks.changed=1



    def OnPlot(self,e):
        if not ud.isempty(self.massdat):
            dlg=AutocorrWindow(self)
            dlg.InitUI(self.config,self.massdat)
            dlg.ShowModal()

    def OnRightClick(self,event):
            if not hasattr(self, "popupID1"):
                self.popupID1 = wx.NewId()
                self.Bind(wx.EVT_MENU, self.OnPopupOne, id=self.popupID1)
            menu = wx.Menu()
            menu.Append(self.popupID1, "Delete")
            self.PopupMenu(menu)
            menu.Destroy()

    def OnPopupOne(self, event):
        item=self.masslistbox.list.GetFirstSelected()
        num=self.masslistbox.list.GetSelectedItemCount()
        self.selection=[]
        self.selection.append(item)
        for i in range(1,num):
            item=self.masslistbox.list.GetNextSelected(item)
            self.selection.append(item)
        for i in range(0,num):
            self.masslistbox.list.DeleteItem(self.selection[num-i-1])

    def OnRightClick2(self,event):
            if not hasattr(self, "popupID1"):
                self.popupID2 = wx.NewId()
                self.Bind(wx.EVT_MENU, self.OnPopup2, id=self.popupID2)
            menu = wx.Menu()
            menu.Append(self.popupID2, "Delete")
            self.PopupMenu(menu)
            menu.Destroy()

    def OnPopup2(self, event):
        item=self.oligomerlistbox.list.GetFirstSelected()
        num=self.oligomerlistbox.list.GetSelectedItemCount()
        self.selection=[]
        self.selection.append(item)
        for i in range(1,num):
            item=self.oligomerlistbox.list.GetNextSelected(item)
            self.selection.append(item)
        for i in range(0,num):
            self.oligomerlistbox.list.DeleteItem(self.selection[num-i-1])

    def OnMatchI(self,e):
        self.oligos = self.oligomerlistbox.list.GetList()
        self.oligomasslist = []
        self.oligonames=[]
        for i in range(0, len(self.oligos)):
            start = int(self.oligos[i][2])
            end = int(self.oligos[i][3])
            for j in range(start, end + 1, 1):
                newmass = float(self.oligos[i][0]) + j * float(self.oligos[i][1])
                if (newmass > 0):
                    self.oligomasslist.append(newmass)
                    if j>0 or self.oligos[i][4]=="":
                        self.oligonames.append(str(j)+""+self.oligos[i][4])
                    else:
                        self.oligonames.append("")
                    #self.oligonames.append(str(j)+""+self.oligos[i][4])
        self.oligomasslist = np.array(self.oligomasslist)
        #self.oligomasslist = np.unique(self.oligomasslist)
        self.Match()

    def OnMatchAll(self,e):
        def lengths(list):
            top = []
            num = len(list)
            for i in range(0, num):
                start = int(list[i, 2])
                end = int(list[i, 3])
                top.append(end + 1 - start)
            return top

        def combine(array2):
            lens = lengths(array2)

            tup = tuple(lens)
            startindex = array2[:, 2]
            startindex=startindex.astype(np.int)
            basemass = array2[:, 0]
            basemass=basemass.astype(np.float)
            omass = array2[:, 1]
            omass=omass.astype(np.float)
            names= array2[:,4]

            finlist = []
            namelist=[]
            for index in np.ndindex(tup):
                name=""
                for i in range(0,len(index)):
                    val=index[i]+startindex[i]
                    if val>0:
                        name=name+str(val)+""+names[i]+" "
                    else:
                        pass
                sum = np.sum((index + startindex) * omass + basemass)
                if (sum > 0):
                    finlist.append(sum)
                    namelist.append(name)
                    #print index,name

            return finlist,namelist

        self.oligos = self.oligomerlistbox.list.GetList()
        #print self.oligos
        if(len(self.oligos)>1):
            self.oligos=np.array(self.oligos)
            self.oligomasslist,self.oligonames = combine(self.oligos)
        else:
            self.oligomasslist = []
            self.oligonames=[]
            for i in range(0, len(self.oligos)):
                start = int(self.oligos[i][2])
                end = int(self.oligos[i][3])
                for j in range(start, end + 1, 1):
                    newmass = self.oligos[i][0] + j * self.oligos[i][1]
                    if (newmass > 0):
                        self.oligomasslist.append(newmass)
                        if j>0 or self.oligos[i][4]=="":
                            self.oligonames.append(str(j)+""+self.oligos[i][4])
                        else:
                            self.oligonames.append("")

        self.oligomasslist = np.array(self.oligomasslist)
        #self.oligomasslist = np.unique(self.oligomasslist)
        self.Match()

    def Match(self):
        self.matches=[]
        self.error=[]
        self.peaks=[]
        self.names=[]
        for i in range(0,self.pks.plen):
            p=self.pks.peaks[i]
            target=p.mass
            nearpt=ud.nearestunsorted(self.oligomasslist,target)
            match=self.oligomasslist[nearpt]
            name=self.oligonames[nearpt]
            p.label=name
            p.match=match
            p.matcherror=target-match
            self.matches.append(match)
            self.error.append(target-match)
            self.peaks.append(target)
            self.names.append(name)
        self.matchlist=[self.peaks,self.matches,self.error,self.names]
        self.matchlistbox.list.Populate(self.matchlist[0],self.matchlist[1],self.matchlist[2],self.matchlist[3])

    def OnClose(self, e):
        self.newmasslist = self.masslistbox.list.GetList()
        #print self.newmasslist
        if not ud.isempty(self.newmasslist):
            self.newmasslist=np.array(self.newmasslist)
            self.newmasslist = self.newmasslist[self.newmasslist>0]
            #print self.newmasslist
            if len(self.newmasslist)>0:
                self.config.masslist=self.newmasslist
            else:
                self.config.masslist=[]
        self.oligos = self.oligomerlistbox.list.GetList()
        if not ud.isempty(self.oligos):
            self.oligos=np.array(self.oligos)
            self.oligoshort = self.oligos[:,:2]
            self.oligoshort=self.oligoshort.astype(np.float)
            self.oligos=self.oligos[np.any([self.oligoshort!=0],axis=2)[0],:]
            #self.oligos=self.oligos[::-1]
        self.config.oligomerlist=self.oligos

        if not ud.isempty(self.matchlist):
            self.newmatchlist=self.matchlist

        self.Destroy()
        try:
            self.EndModal(0)
        except:
            pass

    def OnCloseCancel(self, e):
        self.newmasslist = self.defaultmasslist
        self.oligos = self.defaultoligolist
        self.newmatchlist=self.defaultmatchlist
        self.Destroy()
        self.EndModal(1)

    def OnPopulateButton(self, e):
        try:
            self.masslistbox.list.Populate(self.pks.masses)
        except:
            print "Pick peaks first"

    def OnPopulateButton2(self, e):
        self.oligos = self.oligomerlistbox.list.GetList()
        self.oligomasslist = []
        for i in range(0, len(self.oligos)):
            start = int(self.oligos[i][2])
            end = int(self.oligos[i][3])
            for j in range(start, end + 1, 1):
                newmass = float(self.oligos[i][0]) + j * float(self.oligos[i][1])
                if (newmass > 0):
                    self.oligomasslist.append(newmass)
        self.oligomasslist = np.array(self.oligomasslist)
        self.oligomasslist = np.unique(self.oligomasslist)
        self.masslistbox.list.Populate(self.oligomasslist)



    def OnPopulateButton3(self, e):

        def lengths(list):
            top = []
            num = len(list)
            for i in range(0, num):
                start = int(list[i, 2])
                end = int(list[i, 3])
                top.append(end + 1 - start)
            return top

        def combine(array2):
            lens = lengths(array2)

            tup = tuple(lens)
            startindex = array2[:, 2]
            startindex=startindex.astype(np.int)
            basemass = array2[:, 0]
            basemass=basemass.astype(np.float)
            omass = array2[:, 1]
            omass=omass.astype(np.float)
            finlist = []
            for index in np.ndindex(tup):

                sum = np.sum((index + startindex) * omass + basemass)
                if (sum > 0):
                    finlist.append(sum)
            return finlist

        self.oligos = self.oligomerlistbox.list.GetList()

        if(len(self.oligos)>1):
            self.oligos=np.array(self.oligos)
            self.oligomasslist = combine(self.oligos)
        else:
            self.oligomasslist = []
            for i in range(0, len(self.oligos)):
                start = int(self.oligos[i][2])
                end = int(self.oligos[i][3])
                for j in range(start, end + 1, 1):
                    newmass = self.oligos[i][0] + j * self.oligos[i][1]
                    if (newmass > 0):
                        self.oligomasslist.append(newmass)
        self.oligomasslist = np.array(self.oligomasslist)
        self.oligomasslist = np.unique(self.oligomasslist)
        self.masslistbox.list.Populate(self.oligomasslist)

    def OnClear(self, e):
        self.masslistbox.list.Clear()

    def OnAdd(self, e):
        self.masslistbox.list.AddLine()

    def OnAdd2(self, e):
        self.oligomerlistbox.list.AddLine()

    def OnClear2(self, e):
        self.oligomerlistbox.list.Clear()

    def OnImport(self,e):
        self.mfilename= FileDialogs.open_file_dialog("Open Mass File (mfile)",file_types="*.*")
        if(self.mfilename!=None):
            self.importmass=np.loadtxt(self.mfilename)
            if(self.importmass[0]>0):
                self.masslistbox.list.Populate(self.importmass)

    def OnImport2(self,e):
        self.ofilename= FileDialogs.open_file_dialog("Open Oligomer File (ofile)",file_types="*.*")
        if(self.ofilename!=None):
            self.importolig=np.genfromtxt(self.ofilename,dtype='str')
            if(np.shape(self.importolig)==(5,) or np.shape(self.importolig)==(4,)):
                self.importolig=[self.importolig]
            self.oligomerlistbox.list.Populate(self.importolig)


class TestListCtrl3(wx.ListCtrl,
                    listmix.ListCtrlAutoWidthMixin,
                    listmix.TextEditMixin):
    def __init__(self, parent, ID, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0,imflag=0):
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        self.imflag=imflag
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        listmix.TextEditMixin.__init__(self)

        if self.imflag==0:
            self.InsertColumn(0, "m/z value (Th)")
            self.InsertColumn(1, "+/-")
            self.InsertColumn(2, "Charge")
            self.SetColumnWidth(0, 100)
            self.SetColumnWidth(1, 100)
            self.SetColumnWidth(2, 100)
        else:
            self.InsertColumn(0, "m/z value (Th)",wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(1, "+/-")
            self.InsertColumn(2, "Arrival Time (ms)",wx.LIST_FORMAT_RIGHT)
            self.InsertColumn(3, "+/-")
            self.InsertColumn(4, "Charge")
            self.SetColumnWidth(0, 90)
            self.SetColumnWidth(1, 50)
            self.SetColumnWidth(2, 105)
            self.SetColumnWidth(3, 50)
            self.SetColumnWidth(4, 75)

    def Clear(self):
        self.DeleteAllItems()

    def AddLine(self,line=[0,0,0,0]):
        index=self.InsertStringItem(sys.maxint, str(line[0]))
        self.SetStringItem(index, 1, str(line[1]))
        self.SetStringItem(index, 2, str(line[2]))
        if self.imflag==1:
            self.SetStringItem(index, 3, str(line[3]))
            self.SetStringItem(index, 4, str(0))


    def Populate(self, data,colors=None):
        self.DeleteAllItems()
        for i in range(0, len(data)):
            index = self.InsertStringItem(sys.maxint, str(data[i][0]))
            self.SetStringItem(index, 1, str(data[i][1]))
            self.SetStringItem(index, 2, str(data[i][2]))
            if self.imflag==1:
                self.SetStringItem(index, 3, str(data[i][3]))
                self.SetStringItem(index, 4, str(data[i][4]))
            if colors is not None:
                color=colors[i]
                color=wx.Colour(round(color[0]*255),round(color[1]*255),round(color[2]*255),alpha=255)
                self.SetItemBackgroundColour(index,col=color)

    def GetList(self):
        count = self.GetItemCount()
        list = []
        for i in range(0, count):
            if self.imflag==0:
                sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),float(self.GetItem(i, col=2).GetText())]
            else:
                sublist = [float(self.GetItem(i, col=0).GetText()), float(self.GetItem(i, col=1).GetText()),float(self.GetItem(i, col=2).GetText()),float(self.GetItem(i, col=3).GetText()),float(self.GetItem(i, col=4).GetText())]
            list.append(sublist)
        return list


class TestListCtrlPanel3(wx.Panel):
    def __init__(self, parent,imflag=0):
        wx.Panel.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        tID = wx.NewId()
        sizer = wx.BoxSizer(wx.VERTICAL)
        if wx.Platform == "__WXMAC__" and \
                hasattr(wx.GetApp().GetTopWindow(), "LoadDemo"):
            self.useNative = wx.CheckBox(self, -1, "Use native listctrl")
            self.useNative.SetValue(
                not wx.SystemOptions.GetOptionInt("mac.listctrl.always_use_generic"))
            self.Bind(wx.EVT_CHECKBOX, self.OnUseNative, self.useNative)
            sizer.Add(self.useNative, 0, wx.ALL | wx.ALIGN_RIGHT, 4)
        self.list = TestListCtrl3(self, tID, size=(300+100*imflag,150),style=wx.LC_REPORT,imflag=imflag)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnRightClick, self.list)
        sizer.Add(self.list, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.SetAutoLayout(True)

    def OnRightClick(self,e):
        if not hasattr(self, "popupID1"):
            self.popupID1 = wx.NewId()
            self.Bind(wx.EVT_MENU, self.OnPopupOne, id=self.popupID1)
        menu = wx.Menu()
        menu.Append(self.popupID1, "Delete")
        self.PopupMenu(menu)
        menu.Destroy()

    def OnPopupOne(self,e):
        #Delete
        item=self.list.GetFirstSelected()
        num=self.list.GetSelectedItemCount()
        self.selection=[]
        self.selection.append(item)
        for i in range(1,num):
            item=self.list.GetNextSelected(item)
            self.selection.append(item)
        for i in range(0,num):
            self.list.DeleteItem(self.selection[num-i-1])

    def OnUseNative(self, event):
        wx.SystemOptions.SetOptionInt("mac.listctrl.always_use_generic", not event.IsChecked())
        wx.GetApp().GetTopWindow().LoadDemo("ListCtrl_edit")

def closest(x,y,manlist):
    u=np.array([x,y])
    mindist=sys.maxint
    for i,l in enumerate(manlist):
        v=np.array([l[0],l[2]])
        dist=euclidean(u,v)
        if dist<mindist:
            mindist=dist
            out=l
            pos=i
    return out,pos
#NEED FIXING
def correctassignments(manlist,xdat,ydat):
    manlist2=[]
    for l in manlist:
        list=[l[0],l[0],l[2],l[2],l[4]]
        manlist2.append(list)
    manlist2=np.array(manlist2)

    for i,x in enumerate(xdat):
        for j,y in enumerate(ydat):
            l,pos=closest(x,y,manlist)
            if abs(x-l[0])<l[1] and abs(y-l[2])<l[3]:
                if x<manlist2[pos,0]:
                    manlist2[pos,0]=x
                if x>manlist2[pos,1]:
                    manlist2[pos,1]=x
                if y<manlist2[pos,2]:
                    manlist2[pos,2]=y
                if y>manlist2[pos,3]:
                    manlist2[pos,3]=y
    manlist3=[]
    for l in manlist2:
        xwin=(l[1]-l[0])/2.
        ywin=(l[3]-l[2])/2.
        list=[l[0]+xwin,xwin,l[2]+ywin,ywin,l[4]]
        manlist3.append(list)
    return np.array(manlist3)

def range_overlap(a_min, a_max, b_min, b_max):
    '''Neither range is completely greater than the other
    '''
    return not ((a_min > b_max) or (b_min > a_max))

def checkoverlap(l1,l2):
    return range_overlap(l1[0],l1[1],l2[0],l2[1]) and range_overlap(l1[2],l1[3],l2[2],l2[3])

def detectoverlap(manlist):
    manlist2=[]
    for l in manlist:
        list=[l[0]-l[1],l[0]+l[1],l[2]-l[3],l[2]+l[3],l[4]]
        manlist2.append(list)
    manlist2=np.array(manlist2)

    for l1 in manlist2:
        for l2 in manlist2:
            if np.all(l1 != l2) and checkoverlap(l1,l2):
                print "Overlap of rectangles detected"
                print l1,l2
                return True
    return False


class ManualSelection(wx.Dialog):
    def __init__(self,*args, **kwargs):
        wx.Dialog.__init__(self,style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER,*args, **kwargs)
        self.SetTitle("Manually Assign Masses")

    def InitUI(self, config,data):
        self.data=data
        self.SetSize((550+config.imflag*50, 600))
        self.config=config
        self.defaulttrunclist = deepcopy(self.config.manuallist)
        self.newtrunclist = deepcopy(self.defaulttrunclist)

        self.pnl = wx.Panel(self)

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        size=(7,4.5)
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        if self.config.imflag==0:
            self.plot1= plot1d.Plot1d(self.pnl,figsize=size)
        else:
            self.plot1= plot2d.Plot2d(self.pnl,figsize=size)

        self.vbox2.Add(self.plot1,0,wx.EXPAND)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.sb = wx.StaticBox(self.pnl, label='Manually Set Masses')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.importbutton = wx.Button(self.pnl, label="Import from File")
        self.Bind(wx.EVT_BUTTON, self.OnImport, self.importbutton)

        self.clearbutt = wx.Button(self.pnl, label="Clear List")
        self.Bind(wx.EVT_BUTTON, self.OnClear, self.clearbutt)

        self.addbutton = wx.Button(self.pnl, label="Manual Add Species")
        self.Bind(wx.EVT_BUTTON, self.OnAdd, self.addbutton)

        self.addbutton2 = wx.Button(self.pnl, label="Add from Plot Zoom Range")
        self.Bind(wx.EVT_BUTTON, self.OnAddFromPlot, self.addbutton2)

        self.plotbutton = wx.Button(self.pnl, label="Plot Manual Assignments")
        self.Bind(wx.EVT_BUTTON, self.OnPlot, self.plotbutton)

        self.sbs.Add(self.importbutton,0,wx.EXPAND)
        self.sbs.Add(self.addbutton,0,wx.EXPAND)
        self.sbs.Add(self.addbutton2,0,wx.EXPAND)
        self.sbs.Add(self.clearbutt,0,wx.EXPAND)
        self.sbs.Add(self.plotbutton,0,wx.EXPAND)
        self.hbox.Add(self.sbs)

        self.sb2 = wx.StaticBox(self.pnl, label='Manual List')
        self.sbs2 = wx.StaticBoxSizer(self.sb2, orient=wx.VERTICAL)
        self.masslistbox = TestListCtrlPanel3(self.pnl,imflag=self.config.imflag)
        #self.sbs2.Add(wx.StaticText(self.pnl, label="Manual List"))
        self.sbs2.Add(self.masslistbox)
        self.hbox.Add(self.sbs2)
        self.vbox2.Add(self.hbox,1)
        self.pnl.SetSizer(self.vbox2)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        hboxend.Add(okButton)
        hboxend.Add(closeButton, flag=wx.LEFT, border=5)

        self.vbox.Add(self.pnl, proportion=1,
                      flag=wx.ALL | wx.EXPAND, border=5)
        self.vbox.Add(hboxend,
                      flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.SetSizer(self.vbox)

        okButton.Bind(wx.EVT_BUTTON, self.OnClose)
        closeButton.Bind(wx.EVT_BUTTON, self.OnCloseCancel)

        self.masslistbox.list.Populate(self.defaulttrunclist)
        if self.config.imflag==0:
            self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1],"", "m/z (Th)","Normalized Intensity","Data",self.config)
        else:
            self.plot1.contourplot(self.data,self.config,xlab="m/z (Th)",ylab="Arrival Time (ms)",title="")
        self.CenterOnParent()

    def OnClose(self, e):
        self.newtrunclist = self.masslistbox.list.GetList()
        if not ud.isempty(self.newtrunclist):
            self.newtrunclist=np.array(self.newtrunclist)
            self.newtrunclist=self.newtrunclist[np.any([self.newtrunclist!=0],axis=2)[0],:]

        if (self.newtrunclist != ""):
            if (len(self.newtrunclist) > 0):
                self.config.manuallist = self.newtrunclist
            else:
                self.config.manuallist=[]
        else:
            self.config.manuallist=[]

        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.newtrunclist = self.defaulttrunclist
        self.Destroy()
        self.EndModal(1)

    def OnClear(self, e):
        self.masslistbox.list.Clear()

    def OnAdd(self, e):
        self.masslistbox.list.AddLine()

    def OnAddFromPlot(self,e):
        x0,x1=self.plot1.subplot1.get_xlim()
        xwin=(x1-x0)/2.
        if self.config.imflag==0:
            line=[x0+xwin,xwin,0,0]
        else:
            y0,y1=self.plot1.subplot1.get_ylim()
            ywin=(y1-y0)/2.
            line=[x0+xwin,xwin,y0+ywin,ywin]
        self.masslistbox.list.AddLine(line=line)
        pass

    def OnPlot(self,e):
        if self.config.imflag==0:
            self.plot1.plotrefreshtop(self.data[:, 0], self.data[:, 1],"", "m/z (Th)","Normalized Intensity","Data",self.config)
        else:
            self.plot1.contourplot(self.data,self.config,xlab="m/z (Th)",ylab="Arrival Time (ms)",title="")
        self.newtrunclist = self.masslistbox.list.GetList()
        if not ud.isempty(self.newtrunclist):
            self.newtrunclist=np.array(self.newtrunclist)
            self.newtrunclist=self.newtrunclist[np.any([self.newtrunclist!=0],axis=2)[0],:]
            colormap = cm.get_cmap('rainbow', len(self.newtrunclist))
            xcolors = colormap(np.arange(len(self.newtrunclist)))
            if self.config.imflag==1:
                try:
                    if detectoverlap(self.newtrunclist):
                        self.newtrunclist=correctassignments(self.newtrunclist,np.unique(self.data[:,0]),np.unique(self.data[:,1]))
                except:
                    print "Error with overlapping assignments. Try making sure regions don't intersect."
            for i,l in enumerate(self.newtrunclist):
                if self.config.imflag==0:
                    y0=np.amin(self.data[:,1])
                    ywidth=np.amax(self.data[:,1])-y0
                else:
                    y0=l[2]-l[3]
                    ywidth=l[3]*2.
                self.plot1.subplot1.add_patch(Rectangle((l[0]-l[1],y0),l[1]*2.,ywidth,alpha=0.5,facecolor=xcolors[i],edgecolor='black',fill=True))
            self.plot1.repaint()
            self.masslistbox.list.Populate(self.newtrunclist,colors=xcolors)
        pass

    def OnImport(self,e):
        self.truncfilename= FileDialogs.open_file_dialog("Open File",file_types="*.*")
        if(self.truncfilename!=None):
            self.importtrunc=np.loadtxt(self.truncfilename)
            if self.config.imflag==0:
                if(self.importtrunc.shape==(3,)):
                    self.importtrunc=[self.importtrunc]
                if(len(self.importtrunc[0])==3):
                    self.masslistbox.list.Populate(self.importtrunc)
                    #print self.importtrunc
            else:
                if(self.importtrunc.shape==(5,)):
                    self.importtrunc=[self.importtrunc]
                if(len(self.importtrunc[0])==5):
                    self.masslistbox.list.Populate(self.importtrunc)
                    #print self.importtrunc
