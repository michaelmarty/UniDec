import pyteomics.mass as ms
import numpy as np
from collections import Counter
import molmass as mm
from numpy import fft as fftpack

from unidec.modules.isotopetools import isojim_rna

elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
edict = {}
for i, e in enumerate(elements):
    edict[e] = i

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

rnas = ['A', 'C', 'G', 'U']
# RNA Dictionary of Nucleotide Formulas
# This appears to be with one Phospho group and one hydroxyl group
rna_dict = {
    "A": "C10H12O6N5P",
    "C": "C9H12O7N3P",
    "G": "C10H12O7N5P",
    "U": "C9H11O8N2P"}

# RNA Dictionary of Separated Nucleotide Formulas
# Ordered by C,H,O,N,P
# These are the numbers of each element in the nucleotide
# in the context of an oligonucleotide to my best approximation.
# Does not take into account termininal groups (OH, PO4, etc)
rna_dict_sep = {
    "A": [10,12,5,5,1],
    "C": [9,12,6,3,1],
    "G": [10,12,6,5,1],
    "U": [9,11,7,2,1]}

rna_mass_dict = {
    "A": 347.221,
    "C": 323.209,
    "G": 363.221,
    "U": 324.203}


rnaveragine_comp_numerical = np.array([9.4999, 11.7038, 6.0010, 3.70037, 1])
rnaveragine_mass = 304.6166 #Da

rdict = {}
for i, r in enumerate(rna_dict):
    rdict[r] = i

dnas = ['A', 'C', 'G', 'T']
# DNA Dictionary of Nucleotide Formulas
dna_dict = {
    "A": "C10H12O5N5P",
    "C": "C9H12O6N3P",
    "G": "C10H12O6N5P",
    "T": "C10H13O7N2P"}

ddict = {}
for i, d in enumerate(dna_dict):
    ddict[d] = i


massavgine = 111.1254
avgine = np.array([4.9384, 7.7583, 1.3577, 1.4773, 0.0417])
rnaveragine_comp_numerical = np.array([9.50, 10.70, 6.09, 3.70, 1])
rnaveragine_mass = 321.29163925 #Da
isotopes = np.array(
    [
        [[12.0, 98.90], [13.0033548, 1.10], [0, 0], [0, 0]],
        [[1.0078250, 99.985], [2.0141018, 0.015], [0, 0], [0, 0]],
        [[14.0030740, 99.634], [15.0001089, 0.367], [0, 0], [0, 0]],
        [[15.9949146, 99.772], [16.9991315, 0.038], [17.9991605, 0.2], [0, 0]],
        [
            [31.9720707, 95.02],
            [32.9714585, 0.75],
            [33.9678668, 4.21],
            [35.9760809, 0.02],
        ],
    ]
)
def rna_mass_to_isolist(initial_mass, rna_mass = rnaveragine_mass, rna_comp = rnaveragine_comp_numerical):
    # Get the individual count of each of the C H N O P from the rna sequence based on the rnavergine_mass
    # This will be fed into isojim for the distribution
    value_per = initial_mass / rna_mass
    arr = np.zeros(5)
    for i in range(len(rna_comp)):
        arr[i] = value_per * rna_comp[i]
    print(arr)
    return arr

# @njit(fastmath=True)
def makemass(testmass):
    num = testmass / massavgine * avgine
    intnum = np.array([int(round(n)) for n in num])
    x = intnum * isotopes[:, 0, 0]
    minmassint = np.sum(x)
    formula = ""
    if intnum[0] != 0:
        formula = formula + "C" + str(intnum[0])
    if intnum[1] != 0:
        formula = formula + "H" + str(intnum[1])
    if intnum[2] != 0:
        formula = formula + "N" + str(intnum[2])
    if intnum[3] != 0:
        formula = formula + "O" + str(intnum[3])
    if intnum[4] != 0:
        formula = formula + "S" + str(intnum[4])
    return formula, minmassint, intnum


# @njit(fastmath=True)
def isojim(isolist, length=128):
    """Thanks to Jim Prell for Sketching this Code"""
    buffer = np.zeros(length)
    h = np.array([1, 0.00015, 0, 0])
    c = np.array([1, 0.011, 0, 0])
    n = np.array([1, 0.0037, 0, 0])
    o = np.array([1, 0.0004, 0.002, 0])
    s = np.array([1, 0.0079, 0.044, 0])
    h = np.append(h, buffer)
    c = np.append(c, buffer)
    n = np.append(n, buffer)
    o = np.append(o, buffer)
    s = np.append(s, buffer)

    dt = np.dtype(np.complex128)
    hft = fftpack.rfft(h).astype(dt)
    cft = fftpack.rfft(c).astype(dt)
    nft = fftpack.rfft(n).astype(dt)
    oft = fftpack.rfft(o).astype(dt)
    sft = fftpack.rfft(s).astype(dt)

    numc = isolist[0]
    numh = isolist[1]
    numn = isolist[2]
    numo = isolist[3]
    nums = isolist[4]

    allft = cft**numc * hft**numh * nft**numn * oft**numo * sft**nums

    # with nb.objmode(allift='float64[:]'):
    #    allift = fftpack.irfft(allft)
    allift = np.abs(fftpack.irfft(allft))
    # allift = np.abs(allift)
    allift = allift / np.amax(allift)
    return allift[:length]  # .astype(nb.float64)

def mass_to_dist(mass, isolength=128):
    _, minmassint, isolist = makemass(mass)
    intensities = isojim(isolist, length=isolength)
    #intensities /= np.sum(intensities)
    return intensities

def mass_to_vector(x):
    # x0 = x / 10000000.0
    x1 = x / 1000000.0
    x2 = x / 100000.0
    x3 = x / 10000.0
    x4 = x / 1000.0
    x5 = x / 100.0
    newx = np.array([x1, x2, x3, x4, x5])
    newx = np.clip(newx, 0, 1)
    return newx


# Calculate the isotopic distribution of the peptide
def peptide_to_dist(peptide):
    try:
        formula = ms.Composition(peptide)
        # Formula to string
        fstring = ""
        for k in formula:
            fstring += k + str(formula[k])
        dist = get_dist_from_formula1(fstring)
    except Exception as e:
        dist = None
    return dist

def peptide_to_mass(peptide):
    try:
        formula = ms.Composition(peptide)
        mass = formula.mass()
    except Exception:
        return None
    return mass
# Convert peptide sequence into vector
def peptide_to_vector(peptide):
    try:
        # Get counts of each amino acid
        counts = Counter(peptide)
        vec = np.zeros(len(aas))
        for i, aa in enumerate(aas):
            vec[i] = counts[aa]
    except Exception as e:
        vec = None
    return vec

def get_dist_from_formula1(formula, isolen=128, cutoff=0.001):
    mol = mm.Formula(formula)
    isolist = np.zeros(5)
    for a in mol._elements:
        if a == "C":
            isolist[0] = mol._elements[a][0]
        elif a == "H":
            isolist[1] = mol._elements[a][0]
        elif a == "N":
            isolist[2] = mol._elements[a][0]
        elif a == "O":
            isolist[3] = mol._elements[a][0]
        elif a == "S":
            isolist[4] = mol._elements[a][0]
        else:
            raise ValueError("Unknown atom found in formula")
    dist = isojim(isolist, length=isolen)
    return dist


def get_dist_from_formula(formula, isolen=128, cutoff=0.001):
    # Get the isotopic distribution from a chemical formula
    mol = mm.Formula(formula)
    spec = mol.spectrum().dataframe()["Intensity %"].to_numpy()
    # pad spec to isolen
    if (isolen - len(spec)) > 0:
        spec = np.pad(spec, (0, isolen - len(spec)))
    else:
        spec = spec[:isolen]
    spec /= np.sum(spec)
    spec[spec < cutoff] = 0
    return spec

def formula_to_vector(formula):
    """
    Convert a chemical formula to a vector of atom counts
    :param formula: Chemical formula string
    :return: Vector of atom counts
    """
    vector = np.zeros(len(elements))

    fvals = mm.Formula(formula)._elements

    for f in fvals:
        element = f
        num = fvals[f][0]
        vector[edict[element]] = num

    return vector

# RNA Sequence to Dist
def rnaseq_to_dist(rnaseq, isolen=128, cutoff=0.001):
    formula = ""
    for r in rnaseq:
        formula += rna_dict[r]
    # Get the isotopic distribution of the molecular weight
    dist = get_dist_from_formula(formula, isolen, cutoff)
    return dist

# DNA Sequence to Dist
def dnaseq_to_dist(dnaseq, isolen=128, cutoff=0.001):
    formula = ""
    for d in dnaseq:
        formula += dna_dict[d]
    # Get the isotopic distribution of the molecular weight
    dist = get_dist_from_formula(formula, isolen, cutoff)
    return dist

# RNA Sequence to Vector
def rnaseq_to_vector(rnaseq):
    # Get counts of each base
    counts = Counter(rnaseq)
    vec = np.zeros(len(rnas))
    for i, b in enumerate(rnas):
        vec[i] = counts[b]
    return vec

# DNA Sequence to Vector
def dnaseq_to_vector(dnaseq):
    # Get counts of each base
    counts = Counter(dnaseq)
    vec = np.zeros(len(dnas))
    for i, b in enumerate(dnas):
        vec[i] = counts[b]
    return vec

# RNA mass to dist
def rnamass_to_dist(mass, isolen=128, cutoff=0.001):
    # Get the isotopic distribution of the molecular weight
    dist = get_dist_from_formula(mass, isolen, cutoff)
    return dist


#Take the total mass, then divide by the total number of atoms in the averagine
#Then multiply by
if __name__ == "__main__":
    exit()


    rnaseq = "GUAC"
    dnaseq = "GATC"
    print(rnaseq_to_vector(rnaseq))
    print(dnaseq_to_vector(dnaseq))
    print(rnaseq_to_dist(rnaseq))
    print(dnaseq_to_dist(dnaseq))
    print(len(elements))

