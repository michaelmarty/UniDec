# import sys
# pathroot ="C:\\Python\\glypy"
# sys.path.append(pathroot)
import glypy
from glypy.io import glycoct, iupac
from glypy import plot as gplot
from glypy.enzyme import pathways
from glypy.algorithms import subtree_search, similarity
from glycan_examples import *
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

# enzyme_options = {0: "Sialidase", 1: "B1-4 Galactosidase"}
enzyme_options = {i: list(pathways.glycodigest_rules.keys())[i] for i in range(len(pathways.glycodigest_rules.keys()))}


class GlycoDecEngine:
    def __init__(self):
        self.structure = None
        self.text = None
        self.attrs = None
        self.mass = 0
        self.digest_history = []
        self.ncore = glypy.glycans["N-Linked Core"]

    def load_text(self, text, attrs=None):
        try:
            self.load_glycoct(text, attrs)
        except:
            try:
                self.load_iupac(text, attrs)
            except:
                print("Cannot Read Structure")
        return self.structure

    def load_glycoct(self, text, attrs=None):
        self.text = text
        self.attrs = attrs
        self.structure = glycoct.loads(text)
        self.topstruct = deepcopy(self.structure)
        self.digest_history = []
        return self.structure

    def load_iupac(self, text, attrs=None):
        self.attrs = attrs
        if text[-3:] == "-ol":
            self.text = text[:-3]
        elif text[-4:] == "(?1-" or text[-4:] == "(b1-" or text[-4:] == "(a1-":
            self.text = text[:-4]
        else:
            self.text = text
        self.structure = iupac.loads(self.text, dialect="simple")
        self.topstruct = deepcopy(self.structure)
        self.digest_history = []
        return self.structure

    def plot_struct(self, ax=None, label=True, border=False):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_axes([0, 0, 1, 1])
        else:
            fig = None
        # t = mpl.transforms.Affine2D().scale(1)
        gplot.plot(self.structure, ax=ax, center=True, label=label, layout="balanced", pad=0.5,
                   pad_factor=0)  # , layout_args={"transform":t})
        if not border:
            ax.axis("off")
        ax.set_aspect("equal")
        return fig, ax

    def get_mass(self, ndigits=3):
        self.mass = round(self.structure.mass(), ndigits)
        return self.mass

    def find_residues(self, rname):
        rids = []
        for i, r in enumerate(self.structure.indexed_traversal()):
            if str(r) == str(rname):
                rids.append(r.id)
        return rids

    def remove_monomers(self, ids, specificity="any", terminal=True):
        change_made = False
        for i, r in enumerate(self.structure.indexed_traversal()):
            if r.id in ids:
                plinks = r.parents(links=True)
                clinks = r.children(links=True)
                if len(plinks) == 1 and (not terminal or len(clinks) == 0):
                    pos = plinks[0][0]
                    link = plinks[0][1]
                    connection = str(link.child_position) + "+" + str(link.parent_position)
                    if specificity == "any" or connection == specificity:
                        r.drop_monosaccharide(pos)
                        change_made = True
                        print("Dropping", r)
        return change_made

    def enzyme(self, etype="ABS", overwrite_structure=True):
        if etype != "":
            enzyme = pathways.glycodigest_rules[etype]
            results = list(enzyme(self.structure, refund=False))
            while len(results) > 1:
                digested_structure = results[0][0]
                # print(digested_structure)
                results = list(enzyme(digested_structure, refund=False))
            if len(results) == 1:
                digested_structure = results[0][0]
                # print(digested_structure)
                if overwrite_structure:
                    self.structure = digested_structure
                    self.digest_history.append(etype)
                return True
            else:
                # print("Enzyme Not Active")
                return False
        else:
            return False

    def find_active_enzymes(self):
        active_enzymes = []
        for e in enzyme_options:
            if self.enzyme(enzyme_options[e], overwrite_structure=False):
                active_enzymes.append(enzyme_options[e])
        return active_enzymes

    def test_enzyme_sequence(self, seq):
        cleaved = []
        masses = [self.get_mass()]
        for s in seq:
            result = self.enzyme(s)
            if result:
                masses.append(self.get_mass())
            else:
                masses.append(self.mass)
            cleaved.append(result)
        allc = np.all(cleaved)
        anyc = np.any(cleaved)
        return allc, anyc, masses

    def test_equivalent_structures(self, other):
        return self.structure == other

    def test_subtree(self, other):
        results = subtree_search.subtree_of(other, self.structure, tolerance=1)
        return results

    def test_n_core(self):
        results = self.test_subtree(self.ncore)
        if results == 1:
            return True
        elif results > 1:
            print("Found more than one N-core?:", self.structure)
            return True
        else:
            return False

    def test_forbidden(self, r1name="Neu5Ac", r2name="GlcNAc", linkage="any", terminal=False):
        r1 = glypy.monosaccharides[r1name]
        r2 = glypy.monosaccharides[r2name]
        found = False
        for i, r in enumerate(self.structure.indexed_traversal()):
            s1 = similarity.monosaccharide_similarity(r, r1)
            s1 = s1[0] / s1[1]
            if s1 >= 0.9:
                pres = r.parents()
                cres = r.children()
                plinks = r.parents(links=True)
                clinks = r.children(links=True)

                for p in pres:
                    s2 = similarity.monosaccharide_similarity(p[1], r2)
                    s2 = s2[0] / s2[1]
                    if s2 >= 0.9:
                        found = True
                for p in cres:
                    s2 = similarity.monosaccharide_similarity(p[1], r2)
                    s2 = s2[0] / s2[1]
                    if s2 >= 0.9:
                        found = True
                '''
                if len(plinks) == 1 and (not terminal or len(clinks) == 0):
                    pos = plinks[0][0]
                    link = plinks[0][1]
                    connection = str(link.child_position) + "+" + str(link.parent_position)
                    if linkage == "any" or connection == linkage:
                        found = True'''
        return found


if __name__ == "__main__":
    eng = GlycoDecEngine()
    eng.load_text(iupac3)
    # print(eng.find_active_enzymes())
    print(eng.test_forbidden())
    # print(eng.ncore)
    # print(eng.structure)
    # eng.enzyme()
    # eng.enzyme(etype="sialidase")
    # eng.enzyme(etype="b1-4 galactosidase")
