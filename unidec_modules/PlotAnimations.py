__author__ = 'Michael.Marty'

import numpy as np
from  matplotlib.animation import FuncAnimation
import wx

import plot2d
import plot1d
from unidec_modules import unidecstructure




#class Aniplot(wx.Window):
#    def __init__(self,*args,**kwargs):
#        wx.Window.__init__(self, *args, **kwargs)

x=np.arange(0.0,10.0)
datalist=np.array([np.transpose([x,x]),np.transpose([x,x*x]),np.transpose([x,x*x*x])])


class AnimationWindow(wx.Frame):
    def __init__(self, parent, datalist,config=None,yvals=None,mode="1D",*args,**kwargs):
        wx.Frame.__init__(self, parent, title="Plot Animations",size=(-1,-1))
        if mode=="2D":
            self.mode=2
        else:
            self.mode=1
        if config is None:
            self.config= unidecstructure.UniDecConfig()
            self.config.initialize()
        else:
            self.config=config

        self.datalist=datalist
        self.yvals=yvals

        self.dim=1
        self.pos=-1
        self.play=False

        self.CreateStatusBar(2)
        self.panel=wx.Panel(self)
        self.sizer = wx.BoxSizer(wx.VERTICAL)

        if self.mode==1:
            self.plot= plot1d.Plot1d(self.panel)
        else:
            self.plot= plot2d.Plot2d(self.panel)
        self.sizer.Add(self.plot,0,wx.EXPAND)

        self.controlsizer=wx.BoxSizer(wx.HORIZONTAL)

        self.sb = wx.StaticBox(self.panel, label='Frame Rate (ms/frame)')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)
        frmax=2000

        frmin=1
        self.frslider = wx.Slider(self.panel, wx.ID_ANY,500, frmin, frmax, (30, 60), (250, -1),
            wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS)
        self.frslider.SetTickFreq(100)
        self.sbs.Add(self.frslider,0,wx.EXPAND)
        self.Bind(wx.EVT_COMMAND_SCROLL_THUMBRELEASE,self.update_framerate,self.frslider)
        self.controlsizer.Add(self.sbs,0,wx.EXPAND)

        self.playbutton=wx.ToggleButton(self.panel,label="Play")
        self.nextbutton=wx.Button(self.panel,label="Next")
        self.backbutton=wx.Button(self.panel,label="Back")

        self.controlsizer.Add(self.backbutton,0,wx.EXPAND)
        self.controlsizer.Add(self.playbutton,0,wx.EXPAND)
        self.controlsizer.Add(self.nextbutton,0,wx.EXPAND)

        self.Bind(wx.EVT_TOGGLEBUTTON,self.on_play,self.playbutton)
        self.Bind(wx.EVT_BUTTON,self.on_next,self.nextbutton)
        self.Bind(wx.EVT_BUTTON,self.on_back,self.backbutton)

        self.ctlautoscale=wx.CheckBox(self.panel,label="Autoscale")
        self.controlsizer.Add(self.ctlautoscale,0,wx.EXPAND)

        self.sizer.Add(self.controlsizer,0,wx.EXPAND)

        self.panel.SetSizer(self.sizer)
        self.sizer.Fit(self)

        self.Bind(wx.EVT_CLOSE,self.on_close,self)

        self.init()
        self.Centre()
        self.Show(True)

    def on_close(self,e):
        print "Closing"
        self.animation._stop()
        self.Destroy()


    def update(self,frame_number):
        if self.play:
            self.pos+=1
            self.newplot()
            return 0
        else:
            return 0

    def newplot(self):
        try:
            self.pos=self.pos%len(self.datalist)
            self.data=self.datalist[self.pos]
            if self.yvals is None:
                title=""
            else:
                title=str(self.yvals[self.pos])
            if self.mode==1:
                line=self.plot.subplot1.lines[0]
                self.xlim=self.plot.subplot1.get_xlim()
                self.ylim=self.plot.subplot1.get_ylim()
                line.set_data(self.data[:,0],self.data[:,1])

                self.plot.subplot1.set_title(title)

                autoflag=self.ctlautoscale.GetValue()
                if autoflag:
                    self.plot.subplot1.set_autoscale_on(True)
                    self.plot.subplot1.relim()
                    self.plot.subplot1.autoscale_view(True,True,True)
                else:
                    self.plot.subplot1.set_xlim(self.xlim)
                    self.plot.subplot1.set_ylim(self.ylim)
                self.plot.repaint()
            else:
                self.xlim=self.plot.subplot1.get_xlim()
                self.ylim=self.plot.subplot1.get_ylim()
                self.plot.contourplot(self.data,self.config,xlab="",title=title,repaint=False)
                autoflag=self.ctlautoscale.GetValue()
                if not autoflag:
                    self.plot.subplot1.set_xlim(self.xlim)
                    self.plot.subplot1.set_ylim(self.ylim)
                self.plot.repaint()
        except:
            self.animation._stop()


    def init(self):
        self.pos=0
        self.RefreshPlot()
        self.animation = FuncAnimation(self.plot.figure, self.update,interval=500)
        self.animation._start()

    def on_play(self,e):
        tog=self.playbutton.GetValue()
        if tog:
            self.play=True
        else:
            self.play=False
        pass

    def RefreshPlot(self):
        self.pos=self.pos%len(self.datalist)
        self.data=self.datalist[self.pos]
        if self.yvals is None:
            title=""
        else:
            title=str(self.yvals[self.pos])
        if self.mode==1:
            self.plot.plotrefreshtop(self.data[:,0],self.data[:,1],title,"","","junk",self.config,test_kda=True)
        else:
            self.plot.contourplot(self.data,self.config,xlab="",title=title)

    def on_next(self,e):
        self.pos+=1
        self.newplot()
        pass

    def on_back(self,e):
        self.pos-=1
        self.newplot()
        pass

    def update_framerate(self,e):
        framerate=self.frslider.GetValue()
        #print "Updated framerate to:", framerate
        #self.animation._interval=framerate
        #self.animation.new_frame_seq()
        self.animation._stop()
        self.animation = FuncAnimation(self.plot.figure, self.update,interval=framerate)
        self.animation._start()




#Main App Execution
if __name__=="__main__":
    app = wx.App(False)
    frame = AnimationWindow(None, datalist)
    app.MainLoop()