import numpy as np
import wx
from unidec_modules.biopolymertools import *


class BiopolymerFrame(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, title="Biopolymer Calculator", size=(500, 500))  # ,size=(-1,-1))
        self.parent = parent
        self.mass = ""
        self.type = 0  # Protein/Peptide = 0, RNA = 1, DNA = 2
        self.typerna = 1  # fivend parameter

        self.panel = wx.Panel(self)

        self.ctltype = wx.RadioBox(self.panel, label="Type", choices=["Protein/Peptide", "RNA", "ssDNA", "dsDNA"])
        self.ctltyperna = wx.RadioBox(self.panel, label="RNA or DNA 5' End", choices=["-OH", "MonoP", "TriP"])
        self.ctltype.SetSelection(self.type)
        self.ctltyperna.SetSelection(self.typerna)

        self.ctlseq = wx.TextCtrl(self.panel, value="", size=(450, 200), style=wx.TE_MULTILINE | wx.TE_CHARWRAP)
        self.ctlmass = wx.TextCtrl(self.panel, value="", )

        self.calcbutton = wx.Button(self.panel, -1, "Calculate", size=(450, 25))
        self.Bind(wx.EVT_BUTTON, self.calculate, self.calcbutton)

        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.sizer.Add(self.ctltype, 0)
        self.sizer.Add(self.ctltyperna, 0)

        self.sizer.Add(wx.StaticText(self.panel, label="  Sequence: "))
        self.sizer.Add(self.ctlseq, 0, flag=wx.ALIGN_CENTER_HORIZONTAL)

        self.sizer.Add(self.calcbutton, 0, flag=wx.ALIGN_CENTER_HORIZONTAL)

        self.sizer.Add(wx.StaticText(self.panel, label="\nCalculated Mass (Da): "))
        self.sizer.Add(self.ctlmass, 0)

        hboxend = wx.BoxSizer(wx.HORIZONTAL)
        okbutton = wx.Button(self.panel, label='Ok')
        closebutton = wx.Button(self.panel, label='Cancel')
        okbutton.Bind(wx.EVT_BUTTON, self.close_ok)
        closebutton.Bind(wx.EVT_BUTTON, self.on_close_cancel)

        hboxend.Add(okbutton)
        hboxend.Add(closebutton, flag=wx.LEFT, border=5)

        self.sizer.Add(hboxend, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10)

        self.panel.SetSizer(self.sizer)
        self.panel.Fit()

    def calculate(self, e=0):
        self.seq = self.ctlseq.GetValue()
        self.type = self.ctltype.GetSelection()
        self.typerna = self.ctltyperna.GetSelection()

        if self.typerna == 0:
            fiveend = "OH"
        elif self.typerna == 1:
            fiveend = "MP"
        elif self.typerna == 2:
            fiveend = "TP"
        if self.type == 0:
            self.mass = calc_pep_mass(self.seq)
        elif self.type == 1:
            self.mass = calc_rna_mass(self.seq, fiveend=fiveend)
        elif self.type == 2:
            self.mass = calc_dna_mass(self.seq, fiveend=fiveend)
        elif self.type == 3:
            self.mass = 2 * calc_dna_mass(self.seq, fiveend=fiveend)
        print("Mass:", self.mass)
        self.ctlmass.SetValue(str(self.mass))

    def close_ok(self, e=None):
        self.Destroy()
        try:
            self.EndModal(0)
        except Exception as e:
            pass

    def on_close_cancel(self, e):
        """
        Close the dialog but do not modify any of the values.
        :param e: Unused event
        :return: None
        """
        self.Destroy()
        try:
            self.EndModal(1)
        except Exception as e:
            pass

    # TODO: Add DNA feature
    # TODO: Cleanup buttons and positions
    # TODO: Add some color


if __name__ == "__main__":
    print(mass_HPO4 + mass_HPO4 - mass_O - mass_O + mass_H + mass_H)
    app = wx.App()
    panel = BiopolymerFrame(None)
    panel.Show()
    # panel.draw()
    app.MainLoop()
