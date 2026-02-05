import molmass
import pyteomics.mass as ms
import numpy as np
from collections import Counter
import molmass as mm
from unidec.modules.isotopetools import isojim
from numpy import fft as fftpack
import re

from unidec.modules.isotopetools import isojim_rna

elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
edict = {}
for i, e in enumerate(elements):
    edict[e] = i

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

'''
The format of this is a dictionary of elemental compositions of the 20 common amino acids.
The order of the compositions by element is C, H, N, O, S (as is taken by isojim)
'''
aa_elementalcomp_dict = {
    "A": np.array([3,5,1,1,0]),
    "R": np.array([6,12,4,1,0]),
    "N": np.array([4,6,2,2,0]),
    "D": np.array([4,5,1,3,0]),
    "C": np.array([3,5,1,1,1]),
    "Q": np.array([5,8,2,2,0]),
    "E": np.array([5,7,2,3,0]),
    "G": np.array([2,3,1,1,0]),
    "H": np.array([6,7,3,1,0]),
    "I": np.array([6,11,1,1,0]),
    "L": np.array([6,11,1,1,0]),
    "K": np.array([6,12,2,1,0]),
    "M": np.array([5,9,1,1,1]),
    "F": np.array([9,9,1,1,0]),
    "P": np.array([5,7,1,1,0]),
    "S": np.array([5,7,1,2,0]),
    "T": np.array([4,7,1,2,0]),
    "W": np.array([11,10,2,1,0]),
    "Y": np.array([9,9,1,2,0]),
    "V": np.array([5,9,1,1,0])
}

rnas = ['A', 'C', 'G', 'U']
# RNA Dictionary of Nucleotide Formulas
# This appears to be with one Phospho group and one hydroxyl group
rna_dict = {
    "A": "C10H12O6N5P",
    "C": "C9H12O7N3P",
    "G": "C10H12O7N5P",
    "U": "C9H11O8N2P"}

# RNA Dictionary of Separated Nucleotide Formulas
# Ordered by C,H,N,O,S
# These are the numbers of each element in the nucleotide
# in the context of an oligonucleotide to my best approximation.
# Does not take into account termininal groups (OH, PO4, etc)
rna_dict_sep = {
    "A": np.array([10, 11, 5, 2]),
    "C": np.array([9, 11, 3, 3]),
    "G": np.array([10, 11, 5, 3]),
    "U": np.array([9, 10, 2, 4])}

phospho = np.array([0, 0, 0, 4])

rnaformula_dict_keys = ["C", "H", "N", "O"]


rnaveragine_comp_numerical = np.array([9.50, 10.75, 3.75, 7.0, 0])
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
'''Ordered as C, H, N, O, S, P'''
rnaveragine_comp_forisolist = np.array([9.50, 10.70, 3.70, 7.00, 0])
rnaveragine_mass = 312.9514 #Da
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

formula_pattern = r'([A-Z][a-z]?\d*)'
fullseq_pattern = r'\[.*?\]'

# @njit(fastmath=True)
def pep_makemass(testmass):
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

def rna_makemass(mass):
    num = mass / rnaveragine_mass * rnaveragine_comp_forisolist
    intnum = [int(round(n)) for n in num]
    return intnum


def pepmass_to_dist(mass, isolen=128):
    _, minmassint, isolist = pep_makemass(mass)
    intensities = isojim(isolist, length=isolen)
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


def parse_chemical_formula(formula):
    '''
    Takes a chemical formula string and returns a dictionary of element counts.
    '''
    matches = re.findall(formula_pattern, formula)
    elements = {}
    for match in matches:
        # Split the match into the element and its count (if available)
        element = ''.join([char for char in match if char.isalpha()])  # Extract the element part
        count = ''.join([char for char in match if char.isdigit()])  # Extract the numeric part
        count = int(count) if count else 1
        if element not in elements:
            elements[element] = count
        else:
            elements[element] += count
    return elements


# Calculate the isotopic distribution of the peptide
def peptide_to_dist(peptide, isolen=128):
    try:
        isolist = np.array([0,2,0,1,0,0,0,0,0,0,0])

        mod_matches = re.findall(fullseq_pattern, peptide)
        if len(mod_matches) > 0:
            mod_matches = [match.strip('[]') for match in mod_matches]
            mod_chemical_formula = {}
            for mod in mod_matches:
                current_mod = parse_chemical_formula(mod)
                for element in current_mod:
                    if element not in mod_chemical_formula:
                        mod_chemical_formula[element] = current_mod[element]
                    else:
                        mod_chemical_formula[element] += current_mod[element]

            isolist[0] += mod_chemical_formula.get("C", 0)
            isolist[1] += mod_chemical_formula.get("H", 0)
            isolist[2] += mod_chemical_formula.get("N", 0)
            isolist[3] += mod_chemical_formula.get("O", 0)
            isolist[4] += mod_chemical_formula.get("S", 0)
            isolist[5] += mod_chemical_formula.get("Fe", 0)
            isolist[6] += mod_chemical_formula.get("K", 0)
            isolist[7] += mod_chemical_formula.get("Ca", 0)
            isolist[8] += mod_chemical_formula.get("Ni", 0)
            isolist[9] += mod_chemical_formula.get("Zn", 0)
            isolist[10] += mod_chemical_formula.get("Mg", 0)

            peptide = re.sub(fullseq_pattern, '', peptide)

        for aa in peptide:
            if aa not in aa_elementalcomp_dict:
                print(aa)
                raise ValueError("Unknown amino acid found in peptide")
            isolist += aa_elementalcomp_dict[aa]

        dist = isojim(isolist, isolen)
    except Exception as e:
        dist = None
    return dist


def peptide_to_mass(peptide):
    try:
        mod_matches = re.findall(fullseq_pattern, peptide)
        mod_mass = 0
        for mod in mod_matches:
            mod_mass = ms.Composition(mod.strip('[]')).mass() + mod_mass
        peptide = re.sub(fullseq_pattern, '', peptide)
        formula = ms.Composition(peptide)
        mass = formula.mass()
    except Exception:
        return None
    return mass + mod_mass

def add_nucleotide_to_dict(count, nuc_dict, seq_dict):
    for k,v in nuc_dict.items():
        seq_dict[k] += count * v

    return seq_dict

def rnaseq_to_mass(rna_seq):
    formula_list = np.array([0,2,0,2])
    n = len(rna_seq)

    formula_list += (n-1) * phospho

    for nuc in rna_seq:
        if nuc in rna_dict_sep:
            formula_list += rna_dict_sep[nuc]
        else:
            print("Unexpected nucleotide found: ", nuc)

    formula_string = ""
    for i,atom in enumerate(formula_list):
        formula_string = formula_string + str(rnaformula_dict_keys[i]) + str(formula_list[i]) + " "

    formula_string = formula_string + "P" + str(n-1)
    formula = molmass.Formula(formula_string)

    return formula.monoisotopic_mass

def rnaseq_to_formula(rna_seq):
    formula_list = np.array([0, 2, 0, 2])
    n = len(rna_seq)

    formula_list += (n - 1) * phospho

    for nuc in rna_seq:
        if nuc in rna_dict_sep:
            formula_list += rna_dict_sep[nuc]
        else:
            print("Unexpected nucleotide found: ", nuc)

    formula_string = ""
    for i, atom in enumerate(formula_list):
        formula_string = formula_string + str(rnaformula_dict_keys[i]) + str(formula_list[i]) + " "

    formula_string = formula_string + "P" + str(n - 1)
    return formula_string

# Convert peptide sequence into vector
def peptide_to_vector(peptide):
    try:
        vec = np.zeros(20)
        peptide = re.sub(fullseq_pattern, '', peptide)
        # Get counts of each amino acid
        counts = Counter(peptide)

        for i, aa in enumerate(aas):
            vec[i] = counts[aa]
    except Exception as e:
        vec = None
    return vec


def peptide_to_aacount(peptide):
    peptide = re.sub(fullseq_pattern, '', peptide)
    return len(peptide)


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


def rnaseq_to_isolist(rna_seq):
    isolist = np.zeros(11)
    for nt in rna_seq:
        try:
            isolist += rna_dict_sep[nt]
        except:
            print("Unexpected nucleotide", nt)
    return isolist


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
def rnaseq_to_dist(rnaseq, isolen=128):
    isolist = rnaseq_to_isolist(rnaseq)
    dist = isojim(isolist, isolen)
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
def rnamass_to_dist(mass, isolen=128):
    isolist = rna_makemass(mass)
    dist = isojim(isolist, isolen)
    return dist


def rnamass_to_isolen(mass):
    if mass < 18000:
        return 32
    elif mass > 18000 and mass < 75000:
        return 64
    else:
        return 128


def pepmass_to_isolen(mass):
    if mass < 11000:
        return 32
    elif mass > 11001 and mass < 55000:
        return 64
    elif mass < 120000:
        return 128
    else:
        return -1

#Take the total mass, then divide by the total number of atoms in the averagine
#Then multiply by
if __name__ == "__main__":
    seq = "ACGU"
    mass = rnaseq_to_mass(seq)
    exit()
