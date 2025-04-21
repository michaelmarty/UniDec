import numpy as np
import os
from unidec.IsoDec.IsoGen.isogen_tools import *

class UnimodModification():
    def __init__(self, lines):
        self.id = -1
        self.name = ""
        self.monoisotopic_mass = -1
        #Composition is ordered as: C, H, N, O, S, P, Na, K
        self.composition = {}
        self.sites = np.array([])
        for line in lines:
            #First split on ":"
            firstsplit = line.split(":")
            if firstsplit[0] == "id":
                self.id = int(firstsplit[2].strip())
            elif firstsplit[0] == "name":
                self.name = firstsplit[1].strip()
            elif firstsplit[0] == "xref":
                secondsplit = firstsplit[1].strip().split(" ")
                if "mono_mass" in secondsplit[0]:
                    self.monoisotopic_mass = float(secondsplit[1].strip().replace('"',''))
                if "site" in secondsplit[0]:
                    self.sites = np.append(self.sites, secondsplit[1].strip().replace('"',''))
                if "delta_composition" in secondsplit[0]:
                    print(secondsplit)
                    #Rejoin 1-end into a space-separated string
                    secondsplit[1] = " ".join(secondsplit[1:])
                    self.parse_composition(secondsplit[1].strip())


    def parse_composition(self, composition_string):
        #First remove quotes
        comp = composition_string.replace('"', '')
        elemental_comp = comp.split(" ")
        for i in range(len(elemental_comp)):
            element = "X"
            count = 1
            if ')' in elemental_comp[i]:
                element = elemental_comp[i].split("(")[0]
                count = int(elemental_comp[i].split("(")[1].split(")")[0])
            else:
                element = elemental_comp[i]

            if element not in elements:
                self.composition = None
                return
            if element is not None:
                self.composition[element] = count

            else:
                print("Unknown element:", element)
                self.composition = None

        return


if __name__ == "__main__":
    os.chdir("Z:\\Group Share\\JGP\\IsoGen\\Mods")
    if True:
        with open("unimod.obo.txt", "r") as f:
            lines = f.readlines()
            mods = []
            currmod = []
            for line in lines:
                if line == "[Term]\n" and len(currmod) > 0:
                    mod = UnimodModification(currmod)
                    if mod.composition is not None:
                        mods.append(mod)
                    currmod = []
                else:
                    currmod.append(line)

        np.savez_compressed("unimod_mods.npz", mods=mods)

    if True:
        mods  = np.load("unimod_mods.npz", allow_pickle=True)["mods"]
        elements = []
        for mod in mods:
            for element in mod.composition:
                if element not in elements and len(element) < 3:
                    elements.append(element)
        print(elements)
