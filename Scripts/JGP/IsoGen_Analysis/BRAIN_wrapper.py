from brainpy import isotopic_variants
import matplotlib.pyplot as plt
from unidec.IsoDec.IsoGen.isogen_tools import *
from unidec.IsoDec.IsoGen.isogenc import *

#Ordered C,H,N,O,P
rnaveragine_formulacomp_numerical = np.array([9.50, 10.75, 3.75, 7.0, 1])
rnaveragine_formulamass = 320.28 #Da

def pepmass_to_brain_input(mass):
    num = mass / massavgine * avgine
    intnum = np.array([int(round(n)) for n in num])
    pep = {"C":intnum[0], "H":intnum[1], "N":intnum[2], "O":intnum[3], "S":intnum[4]}
    return pep

def rnamass_to_brain_input(mass):
    num = ((mass + 95.9534) / rnaveragine_formulamass) * rnaveragine_formulacomp_numerical
    intnum = np.array([int(round(n)) for n in num])
    intnum[3] -= 4
    intnum[4] -= 1
    rna = {"C": intnum[0], "H": intnum[1], "N": intnum[2], "O": intnum[3], "P": intnum[4]}
    return rna

def brain_rnamass_to_dist(mass, isolen=128):
    rna = rnamass_to_brain_input(mass)
    isotope_peaks = isotopic_variants(rna, npeaks=isolen, charge=1)
    intensities = np.zeros(isolen)
    max_intensity = 0
    for i in range(len(isotope_peaks)):
        if isotope_peaks[i].intensity > max_intensity:
            max_intensity = isotope_peaks[i].intensity
        intensities[i] = isotope_peaks[i].intensity
    intensities /= max_intensity
    return intensities

def brain_pepmass_to_dist(mass, isolen=128):
    pep = pepmass_to_brain_input(mass)
    isotope_peaks = isotopic_variants(pep, npeaks=isolen, charge=1)
    intensities = np.zeros(isolen)
    max_intensity = 0
    for i in range(len(isotope_peaks)):
        if isotope_peaks[i].intensity > max_intensity:
            max_intensity = isotope_peaks[i].intensity
        intensities[i] = isotope_peaks[i].intensity
    intensities /= max_intensity
    return intensities

def pepseq_to_brain_input(sequence):
    pep = {"C":0, "H":2, "N":0, "O":1, "S":0}

    for aa in sequence:
        if aa not in aa_elementalcomp_dict:
            print(aa)
            raise ValueError("Unknown amino acid found in peptide")
        else:
            add_aa_to_braindict(pep, aa_elementalcomp_dict[aa])

    return pep

def add_aa_to_braindict(braindict, aa_comp):
    braindict["C"] += aa_comp[0]
    braindict["H"] += aa_comp[1]
    braindict["N"] += aa_comp[2]
    braindict["O"] += aa_comp[3]
    braindict["S"] += aa_comp[4]

def brain_pepseq_to_dist(sequence, isolen=128):
    formula_dict = pepseq_to_brain_input(sequence)
    isotope_peaks = isotopic_variants(formula_dict, npeaks=isolen, charge=1)
    intensities = np.zeros(isolen)
    max_intensity = 0
    for i in range(len(isotope_peaks)):
        if isotope_peaks[i].intensity > max_intensity:
            max_intensity = isotope_peaks[i].intensity
        intensities[i] = isotope_peaks[i].intensity
    intensities /= max_intensity
    return intensities

def rnaseq_to_brain_input(rna_seq):
    formula_list = np.array([0, 2, 0, 2])
    n = len(rna_seq)

    formula_list += (n - 1) * phospho

    for nuc in rna_seq:
        if nuc in rna_dict_sep:
            formula_list += rna_dict_sep[nuc]
        else:
            print("Unexpected nucleotide found: ", nuc)

    formula_list = np.append(formula_list, [n-1])
    formula_dict = {"C":formula_list[0],
                    "H":formula_list[1],
                    "N":formula_list[2],
                    "O":formula_list[3],
                    "P":formula_list[4]}

    return formula_dict

def brain_rnaseq_to_dist(sequence, isolen=128):
    rna = rnaseq_to_brain_input(sequence)
    isotope_peaks = isotopic_variants(rna, npeaks=isolen, charge=1)
    intensities = np.zeros(isolen)
    max_intensity = 0
    for i in range(len(isotope_peaks)):
        if isotope_peaks[i].intensity > max_intensity:
            max_intensity = isotope_peaks[i].intensity
        intensities[i] = isotope_peaks[i].intensity
    intensities /= max_intensity
    return intensities


if __name__ == "__main__":
    sequence = "AGCUAGCUAGCU"
    brain_intensities = brain_rnaseq_to_dist(sequence)
    fft_intensities = fft_rna_seq_to_dist(sequence)

    plt.plot(brain_intensities, label="BRAIN")
    plt.plot(fft_intensities, label="FFT")
    plt.legend()
    plt.show()