import numpy as np
import math

aa_masses = {'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766,
             'G': 57.0519, 'H': 137.1411, 'I': 113.1594, 'K': 128.1741, 'L': 113.1594,
             'M': 131.1926, 'N': 114.1038, 'P': 97.1167, 'Q': 128.1307, 'R': 156.1875,
             'S': 87.0782, 'T': 101.1051, 'V': 99.1326, 'W': 186.2132, 'Y': 163.1760}

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


def calc_pep_mass(sequence, allow_float=True, remove_nan=True):
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


if __name__ == "__main__":
    print(mass_HPO4 + mass_HPO4 - mass_O - mass_O + mass_H + mass_H)
