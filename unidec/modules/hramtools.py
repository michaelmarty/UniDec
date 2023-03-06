import numpy as np
import pandas as pd
import time
import re

lipidsearchspace = [["1H", 54, 184], ["12C", 29, 93], ["14N", 0, 1], ["16O", 8, 17], ["31P", 0, 1]]#, ["13C", 0, 3]]
lipidsearchspace = [["1H", 0, 184], ["12C", 0, 93], ["14N", 0, 1], ["16O", 0, 17], ["31P", 0, 1]]#, ["13C", 0, 3]]
smallmolsearch = [["1H", 0, 100], ["12C", 0, 100], ["14N", 0, 5], ["16O", 0, 5], ["31P", 0, 1]]

class HRAMResult:
    def __init__(self):
        self.match_masses = None
        self.match_errors = None
        self.match_ppm = None
        self.match_comp = None
        self.match_formulas = None
        self.iso_keys = []
        self.elem_keys = []

    def to_df(self):
        self.df = pd.DataFrame()
        self.df["Mass"] = self.match_masses
        self.df["Error"] = self.match_errors
        self.df["PPM"] = self.match_ppm
        # self.df["Composition"] = self.match_comp
        self.df["Formula"] = self.match_formulas
        return self.df


class HRAMSearchSpace:
    def __init__(self, searchspace=None, efile=None):
        # Set up the search space, default to lipids
        if searchspace is None:
            searchspace = lipidsearchspace
        searchspace = np.array(searchspace, dtype=object)
        # Set up the element file
        if efile is None:
            efile = "C:\\Python\\UniDec3\\unidec\\bin\\IUPAC-atomic-masses.csv"
        self.elements = pd.read_csv(efile)
        # Set up the element dict
        self.edict = {}
        for n, m, u in self.elements.to_numpy():
            self.edict[n] = m
        # Set up a few key arrays, get the mass of each isotope
        self.isotopes = searchspace[:, 0]
        self.isomasses = np.array([self.edict[e] for e in searchspace[:, 0]])
        self.starts = searchspace[:, 1]
        self.ends = searchspace[:, 2]
        # Empty out positions for some results
        self.result = None
        self.indexes = None
        self.massvals = None

    def limit_search_space(self, mintarget, maxtarget=None, tolerance=100):
        if maxtarget is None:
            maxtarget = mintarget

        # Check for any things where going in a single dimension will go over the target
        startmass = np.sum(self.starts * self.isomasses)
        deltamzmin = tolerance * 1e-6 * mintarget
        deltamzmax = tolerance * 1e-6 * maxtarget
        print(self.ends)
        for i, m in enumerate(self.isomasses):
            counts = np.arange(0, self.ends[i] - self.starts[i] + 1)
            newmasses = counts * m + startmass
            b1 = newmasses < maxtarget + deltamzmax
            if not np.any(b1):
                continue
            maxcount = np.amax(counts[b1])
            self.ends[i] = self.starts[i] + maxcount
        print(self.ends)

        # Check for any things where going in a single dimension will go over the target
        endmass = np.sum(self.ends * self.isomasses)
        print(self.starts)
        for i, m in enumerate(self.isomasses):
            counts = np.arange(0, self.ends[i] - self.starts[i] + 1)
            newmasses = endmass - counts * m
            b1 = newmasses > mintarget + deltamzmin
            if not np.any(b1):
                continue
            maxcount = np.amax(counts[b1])
            self.starts[i] = self.starts[i] + len(newmasses) - maxcount - 1
        print(self.starts)

    def filter_indexes(self):
        pass

    def generate_masses(self):
        lens = self.ends - self.starts + 1
        print("Generating combinations:", lens, np.product(lens))
        self.indexes = np.array(list(np.ndindex(tuple(lens))))
        self.massvals = np.array([np.sum(self.isomasses * (i + self.starts)) for i in self.indexes])

    def match_to_target(self, target, tolerance=5, silent=True):
        self.result = HRAMResult()
        deltamz = tolerance * 1e-6 * target
        b1 = np.abs(self.massvals - target) < deltamz
        if np.any(b1):
            self.result.match_masses = self.massvals[b1]
            self.result.match_errors = self.result.match_masses - target
            self.result.match_ppm = self.result.match_errors / target * 1e6
            self.result.match_comp = self.indexes[b1] + self.starts
            self.result.match_formulas = np.array(["" for m in self.result.match_comp])
            for i, isoname in enumerate(self.isotopes):
                self.result.iso_keys.append(isoname)
                simpname = re.sub('\d', '', isoname)
                self.result.elem_keys.append(simpname)
                fpart = np.core.defchararray.add(" " + isoname + ":", self.result.match_comp[:, i].astype(str))
                self.result.match_formulas = np.core.defchararray.add(self.result.match_formulas, fpart)

            if not silent:
                print(self.result.to_df())
        else:
            if not silent:
                print("Match not found:", target)
        return self.result


def calc_fromula_from_mass(targets, Searcher=None, tolerance=5):
    st = time.perf_counter()
    if Searcher is None:
        Searcher = HRAMSearchSpace()
    Searcher.limit_search_space(np.amin(targets), np.amax(targets), tolerance=tolerance)
    Searcher.generate_masses()
    print("Generated Masses:", time.perf_counter() - st)
    results = [Searcher.match_to_target(target, tolerance=tolerance) for target in targets]
    print("Search Complete:", time.perf_counter() - st)
    return results



if __name__ == "__main__":
    target = 804.57434
    target = 807.5723599471314
    target = 729.2312
    #target = 18.034374
    target = 529.461
    searcher = HRAMSearchSpace()
    # target = 1500
    results = calc_fromula_from_mass([target], Searcher=searcher, tolerance=5)
    print(results[0].to_df())
    print(results[0].match_comp)
    print(results[0].elem_keys, results[0].iso_keys)


