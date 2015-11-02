from twython import Twython
import webbrowser
import wx
from wx.lib.statbmp import GenStaticBitmap as StaticBitmap

__author__ = 'Michael'


class TwitterWindow(wx.Dialog):
    def __init__(self, *args, **kwargs):
        if "pngs" in kwargs:
            self.pngs = kwargs["pngs"]
            del kwargs["pngs"]
        if "codes" in kwargs:
            self.codes = kwargs["codes"]
            del kwargs["codes"]
        if "imflag" in kwargs:
            self.imflag = kwargs["imflag"]
            del kwargs["imflag"]
        else:
            self.imflag = 0
        super(TwitterWindow, self).__init__(*args, **kwargs)

        self.previewsize = 400

        self.SetSize((420, 550))
        self.SetTitle("Twitter Extension - Bringing Mass Spectrometry to the Masses")
        self.APP_KEY = "tQoLvTjNPbeqZGl95ea8rqfvO"
        self.APP_SECRET = "6knaUv912Db37ZWMMSODuxZvmhjNOcxpHRV6YAyVNSvSfQHVz5"

        if self.imflag == 0:
            choices = ["None", "Data and Fit", "Zero-Charge Mass", "m/z Grid", "Individual Peaks", "Mass Grid",
                       "Bar Chart"]
        else:
            choices = [p[1] for p in self.pngs]
            choices = ["None"] + choices

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.sb = wx.StaticBox(self.pnl, label='Tweet a spectrum!')
        self.sbs = wx.StaticBoxSizer(self.sb, orient=wx.VERTICAL)

        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.loginbutton = wx.Button(self.pnl, label="Twitter Log In")
        self.hbox1.Add(wx.StaticText(self.pnl, label=''), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox1.Add(self.loginbutton, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.hbox1.Add(wx.StaticText(self.pnl, label=' User: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.userbox = wx.TextCtrl(self.pnl, value="", style=wx.TE_READONLY, size=(130, 25))
        self.hbox1.Add(self.userbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox1, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.sbs.AddSpacer(10)

        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.inputbox2 = wx.TextCtrl(self.pnl, value="#UniDec", size=(300, 50),
                                     style=wx.TE_CHARWRAP | wx.TE_MULTILINE | wx.TE_NO_VSCROLL)
        self.inputbox2.SetMaxLength(140)
        self.hbox2.Add(wx.StaticText(self.pnl, label='Tweet: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.inputbox2, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.countbox = wx.TextCtrl(self.pnl, value="7", style=wx.TE_READONLY | wx.TE_RIGHT, size=(30, 25))
        self.hbox2.Add(wx.StaticText(self.pnl, label=' '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(self.countbox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox2.Add(wx.StaticText(self.pnl, label=' Characters'), 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox2, 0, wx.ALIGN_CENTER_HORIZONTAL)
        self.sbs.AddSpacer(10)

        self.hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox3.Add(wx.StaticText(self.pnl, label='Image: '), 0, wx.ALIGN_CENTER_VERTICAL)
        self.imagechoice = wx.Choice(self.pnl, -1, (115, 50), choices=choices)
        self.imagechoice.SetSelection(0)
        self.previewbutton = wx.Button(self.pnl, label="Preview")
        self.hbox3.Add(self.imagechoice, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hbox3.Add(self.previewbutton, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox3, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.hbox4 = wx.BoxSizer(wx.HORIZONTAL)
        self.emptyimg = wx.EmptyImage(self.previewsize, self.previewsize)
        self.emptyimg.Replace(0, 0, 0, 255, 255, 255)
        self.imageCtrl = StaticBitmap(self.pnl, wx.ID_ANY, wx.BitmapFromImage(self.emptyimg))
        self.hbox4.Add(self.imageCtrl, 0)
        self.sbs.Add(self.hbox4, 1, wx.ALIGN_CENTER_HORIZONTAL)

        self.hbox5 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox5.Add(wx.StaticText(self.pnl, label=''), 0, wx.ALIGN_CENTER_VERTICAL)
        self.tweetbutton = wx.Button(self.pnl, label="Tweet it!")
        self.hbox5.Add(self.tweetbutton, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sbs.Add(self.hbox5, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.pnl.SetSizer(self.sbs)
        self.pnl.SetBackgroundColour("WHITE")

        self.vbox.Add(self.pnl, proportion=1,
                      flag=wx.ALL | wx.EXPAND, border=5)

        '''
        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Cancel')
        hboxend.Add(okButton)
        hboxend.Add(closeButton, flag=wx.LEFT, border=5)
        okButton.Bind(wx.EVT_BUTTON, self.on_close)
        closeButton.Bind(wx.EVT_BUTTON, self.on_close_cancel)
        self.vbox.Add(hboxend,
                      flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)
        '''

        self.SetSizer(self.vbox)

        self.loginbutton.Bind(wx.EVT_BUTTON, self.OnLaunchWeb)
        self.inputbox2.Bind(wx.EVT_TEXT, self.OnCharacterCount)
        self.previewbutton.Bind(wx.EVT_BUTTON, self.OnPreview)
        self.tweetbutton.Bind(wx.EVT_BUTTON, self.Tweet)
        self.Center()
        if self.codes is not None:
            self.LoadScreenName()

    def LoadScreenName(self):
        self.OAUTH_TOKEN = self.codes[0]
        self.OAUTH_TOKEN_SECRET = self.codes[1]
        twitter = Twython(self.APP_KEY, self.APP_SECRET,
                          self.OAUTH_TOKEN, self.OAUTH_TOKEN_SECRET)
        self.screen_name = twitter.verify_credentials()["screen_name"]
        print "Logged in successfully as: ", self.screen_name
        self.userbox.SetValue("@" + self.screen_name)

    def OnClose(self, e):
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.Destroy()
        self.EndModal(1)

    def OnPreview(self, e):
        choice = self.imagechoice.GetSelection()
        print choice
        self.imageFile = None
        self.imageCtrl.SetBitmap(wx.BitmapFromImage(self.emptyimg))
        if choice is not 0:
            for i in self.pngs:
                if i[0] == choice:
                    self.imageFile = i[1]
                    print self.imageFile
            if self.imageFile is not None:
                image = wx.Image(self.imageFile, wx.BITMAP_TYPE_ANY)
                W = image.GetWidth()
                H = image.GetHeight()
                if W > H:
                    NewW = self.previewsize
                    NewH = self.previewsize * H / W
                else:
                    NewH = self.previewsize
                    NewW = self.previewsize * W / H
                image = image.Scale(NewW, NewH)
                btm = wx.BitmapFromImage(image)
                self.imageCtrl.SetBitmap(btm)

        pass

    def OnCharacterCount(self, e):
        string = self.inputbox2.GetValue()
        count = len(string)
        self.countbox.SetValue(str(count))

    def OnLaunchWeb(self, e):

        twitter = Twython(self.APP_KEY, self.APP_SECRET)

        auth = twitter.get_authentication_tokens()

        self.OAUTH_TOKEN = auth['oauth_token']
        self.OAUTH_TOKEN_SECRET = auth['oauth_token_secret']

        self.url = auth["auth_url"]
        webbrowser.open(self.url, new=2)

        self.pinwindow = PinWindow(self)
        self.pinwindow.ShowModal()
        self.pin = self.pinwindow.pin

        if self.pin is not None:
            oauth_verifier = self.pin
            twitter = Twython(self.APP_KEY, self.APP_SECRET,
                              self.OAUTH_TOKEN, self.OAUTH_TOKEN_SECRET)

            final_step = twitter.get_authorized_tokens(oauth_verifier)

            self.OAUTH_TOKEN = final_step['oauth_token']
            self.OAUTH_TOKEN_SECRET = final_step['oauth_token_secret']
            self.codes = [self.OAUTH_TOKEN, self.OAUTH_TOKEN_SECRET]
            self.LoadScreenName()

        else:
            print "No pin provided"

    def Tweet(self, e):
        if self.codes is not None:
            self.OAUTH_TOKEN = self.codes[0]
            self.OAUTH_TOKEN_SECRET = self.codes[1]
            twitter = Twython(self.APP_KEY, self.APP_SECRET,
                              self.OAUTH_TOKEN, self.OAUTH_TOKEN_SECRET)
            tweet = self.inputbox2.GetValue()
            print "Tweeting: ", tweet
            choice = self.imagechoice.GetSelection()
            self.imageFile = None
            if choice is not 0:
                self.OnPreview(e)
                # self.imageFile=os.path.join(os.getcwd(),self.imageFile)
                print "\twith image: ", self.imageFile
                photo = open(self.imageFile, "rb")
                # result=twitter.upload_media(media=photo)
                # id=result['media_id']
                # twitter.update_status(status=tweet,media_ids=id)
                twitter.update_status_with_media(status=tweet, media=photo)
            else:
                twitter.update_status(status=tweet)
                pass
        else:
            print "Need to log in to Twitter"
        self.OnClose(e)


class PinWindow(wx.Dialog):
    def __init__(self, *args, **kwargs):
        super(PinWindow, self).__init__(*args, **kwargs)
        self.SetSize((175, 110))
        self.SetTitle("Enter Pin")

        self.pnl = wx.Panel(self)
        self.vbox = wx.BoxSizer(wx.VERTICAL)

        self.sbs = wx.BoxSizer(wx.VERTICAL)
        self.inputbox = wx.TextCtrl(self.pnl, value="")
        self.sbs.Add(self.inputbox, 1, wx.ALIGN_CENTER_VERTICAL)

        self.pnl.SetSizer(self.sbs)

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
        self.Center()

    def OnClose(self, e):
        try:
            self.pin = int(self.inputbox.GetValue())
        except ValueError:
            self.pin = None
        self.Destroy()
        self.EndModal(0)

    def OnCloseCancel(self, e):
        self.pin = None
        self.Destroy()
        self.EndModal(1)
