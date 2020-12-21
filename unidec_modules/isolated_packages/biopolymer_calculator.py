import numpy as np
import wx

aa_masses = {'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766,
             'G': 57.0519, 'H': 137.1411, 'I': 113.1594, 'K': 128.1741, 'L': 113.1594,
             'M': 131.1926, 'N': 114.1038, 'P': 97.1167, 'Q': 128.1307, 'R': 156.1875,
             'S': 87.0782, 'T': 101.1051, 'V': 99.1326, 'W': 186.2132, 'Y': 163.1760}

rna_masses = {'A': 329.2, 'U': 306.2, 'C': 305.2, 'G': 345.2, 'T': 306.2}

dna_masses = {'A': 313.2, 'T': 304.2, 'C': 289.2, 'G': 329.2, 'U': 304.2, }

sequence = "testtest"
s2 = "MKTVVLAVAVLFLTGSQARHFWQRDDPQTPWDRVKDFATVYVDAVKDSGREYVSQFETSALGKQLNLNLLENWDTLGSTVGRLQEQLGPVTQEFWDNLEKETEWLRREMNKDLEEVKAKVQPYLDQFQTKWQEEVALYRQKMEPLGAELRDGARQKLQELQEKLTPLGEDLRDRMRHHVDALRTKMTPYSDQMRDRLAERLAQLKDSPTLAEYHTKAADHLKAFGEKAKPALEDLRQGLMPVFESFKTRIMSMVEEASKKLNAQ"

mass_water = 18.0153
mass_OH = 17.008
mass_O = 15.9994
mass_HPO4 = 95.9793
mass_H = 1.00794


def get_aa_mass(letter):
    try:
        return aa_masses[letter]
    except:
        print("Bad Amino Acid Code:", letter)
        return 0


def get_rna_mass(letter):
    if letter == "T":
        print("Assuming T means U")

    try:
        return rna_masses[letter]
    except:
        print("Bad RNA Code:", letter)
        return 0

def get_dna_mass(letter):
    try:
        return dna_masses[letter]
    except:
        print("Bad DNA Code:", letter)
        return 0


def calc_pep_mass(sequence):
    seq = sequence.upper()
    mass = np.sum([get_aa_mass(s) for s in seq]) + mass_water
    return np.round(mass, 2)


def calc_rna_mass(sequence, threeend="OH", fiveend="MP"):
    seq = sequence.upper()
    mass = np.sum([get_rna_mass(s) for s in seq])
    if threeend == "OH":
        mass += mass_OH

    if fiveend == "OH":
        mass -= mass_HPO4
        mass += mass_OH
    elif fiveend == "MP":
        mass += mass_H
    elif fiveend == "TP":
        mass += mass_HPO4 + mass_HPO4 - mass_O - mass_O + mass_H

    return round(mass, 2)

def calc_dna_mass(sequence, threeend="OH", fiveend="MP"):
    seq = sequence.upper()
    mass = np.sum([get_dna_mass(s) for s in seq])
    if threeend == "OH":
        mass += mass_OH

    if fiveend == "OH":
        mass -= mass_HPO4
        mass += mass_OH
    elif fiveend == "MP":
        mass += mass_H
    elif fiveend == "TP":
        mass += mass_HPO4 + mass_HPO4 - mass_O - mass_O + mass_H

    return round(mass, 2)

class BiopolymerFrame(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, title="Biopolymer Calculator", size=(500, 500))  # ,size=(-1,-1))
        self.parent = parent
        self.mass = ""
        self.type = 0 #Protein/Peptide = 0, RNA = 1, DNA = 2
        self.typerna = 1 #fivend parameter

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

        self.sizer.Add(self.ctltype, 0, flag=wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(self.ctltyperna, 0, flag=wx.ALIGN_CENTER_VERTICAL)

        self.sizer.Add(wx.StaticText(self.panel, label="  Sequence: "))
        self.sizer.Add(self.ctlseq, 0, flag=wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_CENTER_HORIZONTAL)

        self.sizer.Add(self.calcbutton, 0, flag=wx.ALIGN_CENTER_VERTICAL| wx.ALIGN_CENTER_HORIZONTAL)

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
            self.mass = 2*calc_dna_mass(self.seq, fiveend=fiveend)
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
