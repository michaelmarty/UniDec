import numpy as np
import math
import re
import os

aa_masses = {'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766,
             'G': 57.0519, 'H': 137.1411, 'I': 113.1594, 'K': 128.1741, 'L': 113.1594,
             'M': 131.1926, 'N': 114.1038, 'P': 97.1167, 'Q': 128.1307, 'R': 156.1875,
             'S': 87.0782, 'T': 101.1051, 'V': 99.1326, 'W': 186.2132, 'Y': 163.1760}

aa_masses_monoisotopic = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,
                          'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,
                          'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
                          'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333}

rna_masses = {'A': 329.2, 'U': 306.2, 'C': 305.2, 'G': 345.2, 'T': 306.2}

dna_masses = {'A': 313.2, 'T': 304.2, 'C': 289.2, 'G': 329.2, 'U': 304.2, }

s1 = "testtest"
s2 = "MKTVVLAVAVLFLTGSQARHFWQRDDPQTPWDRVKDFATVYVDAVKDSGREYVSQFETSALGKQLNLNLLENWDTLGSTVGRLQEQLGPVTQEFWDNLEKETEW" \
     "LRREMNKDLEEVKAKVQPYLDQFQTKWQEEVALYRQKMEPLGAELRDGARQKLQELQEKLTPLGEDLRDRMRHHVDALRTKMTPYSDQMRDRLAERLAQLKDSPTL" \
     "AEYHTKAADHLKAFGEKAKPALEDLRQGLMPVFESFKTRIMSMVEEASKKLNAQ"

mass_water = 18.0153
mass_OH = 17.008
mass_O = 15.9994
mass_HPO4 = 95.9793
mass_H = 1.00794
mass_proton = 1.00727647


def get_aa_mass(letter):
    if letter == " " or letter == "\t" or letter == "\n":
        return 0

    try:
        return aa_masses[letter]
    except Exception as exception:
        print("Bad Amino Acid Code:", letter)
        return 0


def get_rna_mass(letter):
    if letter == "T":
        print("Assuming T means U")

    try:
        return rna_masses[letter]
    except Exception as exception:
        print("Bad RNA Code:", letter)
        return 0


def get_dna_mass(letter):
    try:
        return dna_masses[letter]
    except Exception as exception:
        print("Bad DNA Code:", letter)
        return 0


def calc_pep_mass(sequence, allow_float=True, remove_nan=True, all_cyst_ox=False, pyroglu=False, round_to=2):
    if all_cyst_ox:
        # Count number of c in sequence
        c = sequence.lower().count("c")
        # Multiply by -1 * mass of H
        modmass = c * (-1 * mass_H)
    else:
        modmass = 0

    if remove_nan:
        if sequence.lower() == "nan":
            return 0.0

    if allow_float:
        try:
            mass = float(sequence)
        except Exception as exception:
            seq = sequence.upper()
            mass = np.sum([get_aa_mass(s) for s in seq]) + mass_water
    else:
        seq = sequence.upper()
        mass = np.sum([get_aa_mass(s) for s in seq]) + mass_water
    # print(sequence, mass)
    # Look for pyroglutamate mod if set
    if pyroglu:
        if sequence[0] == "E":
            modmass -= mass_water
        if sequence[0] == "Q":
            modmass -= mass_OH

    massoutput = np.round(mass + modmass, round_to)
    return massoutput


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

    return round(float(mass), 2)


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

    return round(float(mass), 2)


def read_fasta(path):
    f = open(path, 'r')
    lines = f.readlines()

    hre = re.compile('>(\S+)')
    lre = re.compile('^(\S+)$')

    gene = {}

    for line in lines:
        outh = hre.search(line)
        if outh:
            id = outh.group(1)
        else:
            outl = lre.search(line)
            if id in gene.keys():
                gene[id] += outl.group(1)
            else:
                gene[id] = outl.group(1)
    return gene


if __name__ == "__main__":
    #print(mass_HPO4 + mass_HPO4 - mass_O - mass_O + mass_H + mass_H)
    #os.chdir("..\\..\\Scripts\\old\\Jessica")
    #file = "ecoli.fasta"
    #genes = read_fasta(file)
    #print(genes)
    oligo = "aacauucaACgcugucggugAgu"
    mass = calc_rna_mass(oligo)
    print(mass)
