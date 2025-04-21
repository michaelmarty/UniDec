from typing import Match, List

import numpy as np
import torch.cuda
from psims.mzid.components import IonType
from pyteomics import mass as ms
from pyteomics.mass import std_aa_mass, fast_mass
from sympy import sequence, false

#from Scripts.MTM.HR.uvpd_test import matched_peaks
from unidec.IsoDec.match import MatchedPeak


def gen_fragments(sequence):
    return np.array([[sequence[:i], sequence[i:]] for i in range(1, len(sequence) + 1)])


# Need to mod to take an input of ion types, and generate the appropriate lists AND fragmentation for ALL the corresponding
# single inputs**
def generateList(input, ion_type):
    dict = buildSubscriptDict(input)
    subscript_mass = 0
    sequence = sanitizeInput(input)
    specialInput = False
    specs = ['x+1', 'y-1', 'a+1', 'z']
    if ion_type in specs:
        specialInput = True
    ion_type = ion_type[0].lower() + ion_type[1:]
    #just grab the first char
    rip = ion_type
    ion_type = ion_type[0]
    firstHalf = ['a','b','c']
    #if special input and ion type == a or ion type == x or iontype == z, add 1.007276467
    #else we know we hit y-1 and subtract 1.007276467 from each iteration
    valid_types = {'a','b','c','x','y','z'}
    data = []
    if ion_type not in valid_types:
        return None
    count = len(sequence)
    for i in range(1, len(sequence)):
        key = sequence[i - 1] + str(i)
        addedMass = dict.get(key, 0)
        subscript_mass = subscript_mass + addedMass
        calc_sequence = ""
        #check first to see if iontype is in the first half of the second half
        if ion_type in firstHalf:
            total_addedMass = sum({k: v for k, v in dict.items() if int(k[1:]) <= i}.values())
            calc_sequence = sequence[:i]
        elif ion_type in firstHalf and specialInput:
            total_addedMass = sum({k: v for k, v in dict.items() if int(k[1:]) > i}.values()) + 1.007276467
            calc_sequence = sequence[:i]
        else:
            total_addedMass = sum({k: v for k, v in dict.items() if int(k[1:]) > i}.values())
            calc_sequence = sequence[i:]


        if specialInput:
            if ion_type != 'y':
                total_addedMass = ms.fast_mass(calc_sequence, ion_type[0])+total_addedMass + 1.007276467
            else:
                total_addedMass = ms.fast_mass(calc_sequence, ion_type[0]) + total_addedMass - 1.007276467
        else:
            total_addedMass = ms.fast_mass(calc_sequence, ion_type[0]) + total_addedMass
        if i < len(sequence):
            if ion_type in firstHalf:
                if specialInput:
                    data.append(f"{ion_type[0].upper()}{i}{rip[1:3]}:{total_addedMass:.8f}")
                else:
                    data.append(f"{ion_type[0].upper()}{i}:{total_addedMass:.8f}")

            else:
                if specialInput and ion_type[0] != 'z':
                    data.append(f"{ion_type[0].upper()}{len(sequence) - i}{rip[1:3]}:{total_addedMass:.8f}")
                else:
                    data.append(f"{ion_type[0].upper()}{len(sequence) - i}:{total_addedMass:.8f}")
    return data


def buildSubscriptDict(sequence):
    subscript_dict = {}
    i = 0
    ptr = 0
    while i < len(sequence):
        currSub = ""
        if i + 1 < len(sequence) and sequence[i + 1] == '[':
            currChar = sequence[i]
            i += 2
            while i < len(sequence) and sequence[i] != ']':
                if sequence[i].isdigit() or sequence[i] == '.':
                    currSub += sequence[i]
                i += 1
            if i < len(sequence) and sequence[i] == ']':
                key = currChar + str(ptr + 1)
                if currSub:
                    subscript_dict[key] = float(currSub)
                    ptr += 1
            i += 1
        else:
            if sequence[i].isalnum():
                ptr += 1
            i += 1
    return subscript_dict


def buildMatches(mass, ion_list, width):
    possible_matches = []
    for i in range(len(ion_list)):
        currList = ion_list[i].split(':')
        idx_mass = float(currList[1])
        if abs(idx_mass - mass) <= width:
            possible_matches.append(currList[0])

    return possible_matches


def getLists(sequence, ion_types):
    all_results = []
    for ion_type in ion_types:
        ls = generateList(sequence, ion_type)
        all_results.append(ls)
    return all_results

def evalMatchedPeaks(peaks: List[MatchedPeak], sequence, ion_type, tol_type="ppm", mass_tol=20.0) -> None:
    total_matches = 0
    combinations = generateList(sequence, ion_type)
    clean_sequence = sanitizeInput(sequence)
    firstHalf = ['a','b','c']
    fragPositions = []
    for peak in peaks:
        masses = peak.monoisos
        for mass in masses:
            width = mass_tol
            if tol_type == "ppm":
                width = mass * mass_tol / 1e6
            match = buildMatches(mass, combinations, width)
            for matches in match:
                total_matches += 1
                fragPositions.append(matches)
                peak.matchedions.append(matches)
    unique_breaks = []
    for frag in fragPositions:
        front = 0
        val = 0
        if frag[0].upper() in firstHalf or frag[0].lower() in firstHalf and '+' not in frag:
            val = int(frag[1:])
        elif '+' in frag or '-' in frag:
            frag = frag.split('+') if '+' in frag else frag.split('-')
            check = frag[0][0]
            check = check.lower()

            #account for the x and y
            if check != 'a':
                front = int(frag[0][1:])
                val = len(clean_sequence) - front
            else:
                val = int(frag[0][1:])
        #xyz
        else:
            front = int(frag[1:])
            val = len(clean_sequence) - front
        if val not in unique_breaks:
            unique_breaks.append(val)

    return sorted(unique_breaks), total_matches

def searchTD(matchedpeaks, seq, ion_types, tol_type="ppm", mass_tol=20.0, verbose=False):
    unique_frag_positions = []
    for ion_type in ion_types:
        positions, total_fragments = evalMatchedPeaks(matchedpeaks, seq, ion_type, tol_type, mass_tol)
        if verbose:
            print("Total", ion_type, "Fragments:", total_fragments)
        for pos in positions:
            if pos not in unique_frag_positions:
                unique_frag_positions.append(pos)
    clean_seq = sanitizeInput(seq)
    sequence_coverage = len(unique_frag_positions) / len(clean_seq) * 100

    peaks_with_matchedion = 0
    for peak in matchedpeaks:
        if len(peak.matchedions) > 0:
            peaks_with_matchedion += 1
    percentage_peaks_with_matchedion = peaks_with_matchedion / len(matchedpeaks) * 100
    if verbose:
        print("Total sequence coverage:", round(len(unique_frag_positions) / len(clean_seq) * 100, 2), "%")
        print("Percentage of peaks with a matched ion:", round(peaks_with_matchedion / len(matchedpeaks) * 100, 2), "%")
    return sequence_coverage, percentage_peaks_with_matchedion

def sanitizeInput(sequence):
    sequence = sequence.upper()
    seq = ""
    lets = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i in range(len(sequence)):
        if sequence[i] in lets:
            seq += sequence[i]
    return seq


def getAllMatchedPeaks(peaks: List[MatchedPeak], ion_types, sequence):
    masterList = []
    totaluniq = 0
    for ion in ion_types:
        ls,rip = evalMatchedPeaks(peaks, sequence, ion)
        totaluniq += rip
        for idx in ls:
            if idx not in masterList:
                masterList.append(idx)
    clean_seq = sanitizeInput(sequence)
    print(f"Total Unique Peaks: {round((totaluniq/len(peaks) * 100),2)}")
    print(f"Number of fragmentations accounted for: {len(masterList)}")
    print(f"Total Sequence coverage = " + str(len(masterList) / (len(clean_seq)) * 100) + "%")



def calc_truncation_masses(sequence):
    frags = gen_fragments(sequence)
    masses = np.array([[ms.fast_mass(f[0]), ms.fast_mass(f[1])] for f in frags])
    return masses


def calc_topdown_masses(sequence, charge):
    frags = gen_fragments(sequence)
    masses = np.array(
        [[ms.fast_mass(f[0], ion_type="b", charge=charge), ms.fast_mass(f[1], ion_type="y", charge=charge)] for f in
         frags])
    return masses


def calc_trunc_diff(sequence):
    topmass = ms.fast_mass(sequence)
    masses = np.array(calc_truncation_masses(sequence))
    diff = topmass - masses
    return diff


def calc_mass(sequence):
    return ms.fast_mass(sequence)


'''

rho="MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA"


print(calc_trunc_diff(rho))'''

if __name__ == '__main__':
    # rho = 'MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMV'

    seq = "AS[10.5]DG[23.3]"
    ls = generateList(seq,'c')
    print(ls)

    # test all list
    ion_types = ['a', 'b', 'c', 'x', 'y', 'z', 'x+1', 'a+1', 'y-1']
    matched_peaks = [
        MatchedPeak(z=0, mz=0),
        MatchedPeak(z=0, mz=0),
        MatchedPeak(z=0, mz=0),
    ]
    monoiso_values = [88.06366290,300.62263432, 145.0293]
    for peak, monoiso in zip(matched_peaks, monoiso_values):
        peak.monoiso = monoiso
        peak.monoisos.append(monoiso)

    print(evalMatchedPeaks(matched_peaks, seq, 'c'))





