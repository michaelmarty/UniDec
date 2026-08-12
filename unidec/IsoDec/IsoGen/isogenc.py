import ctypes
import numpy as np
import time
from scipy.spatial.distance import cosine

from Scripts.JGP.IsoGen_Analysis.distribution_assessment import calculate_cosinesimilarity
from unidec.IsoDec.IsoGen import isogen_tools as igs
import matplotlib as mpl
import matplotlib.pyplot as plt
from unidec.IsoDec.IsoGen.isogen_tools import get_dist_from_formula
from unidec.modules.isotopetools import fast_calc_averagine_isotope_dist

mpl.use('WxAgg')
np.set_printoptions(precision=6, suppress=True)
isogen_dll_path = r"C:\Python\UniDecDev\unidec\IsoDec\isogen.dll"

isogen_lib = ctypes.CDLL(isogen_dll_path)

isogen_lib.fft_pep_mass_to_dist.restype = ctypes.c_float
isogen_lib.fft_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]

isogen_lib.fft_rna_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_lib.fft_rna_mass_to_dist.restype = ctypes.c_float

isogen_lib.nn_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_lib.nn_pep_mass_to_dist.restype = ctypes.c_float

isogen_lib.fft_pep_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                           ctypes.c_int]
isogen_lib.fft_pep_seq_to_dist.restype = ctypes.c_float

isogen_lib.nn_pep_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                          ctypes.c_int]
isogen_lib.nn_pep_seq_to_dist.restype = ctypes.c_float

isogen_lib.fft_rna_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                           ctypes.c_int]
isogen_lib.fft_rna_seq_to_dist.restype = ctypes.c_float

isogen_lib.nn_rna_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_int]
isogen_lib.nn_rna_seq_to_dist.restype = ctypes.c_float

isogen_lib.nn_rna_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_lib.nn_rna_mass_to_dist.restype =  ctypes.c_float

isogen_lib.brain_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_lib.brain_pep_mass_to_dist.restype = ctypes.c_float

isogen_lib.brain_pep_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                           ctypes.c_int]
isogen_lib.brain_pep_seq_to_dist.restype = ctypes.c_float

isogen_lib.brain_rna_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                           ctypes.c_int]
isogen_lib.brain_rna_seq_to_dist.restype = ctypes.c_float

isogen_lib.brain_rna_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                           ctypes.c_int]
isogen_lib.brain_rna_mass_to_dist.restype = ctypes.c_float



# nn = neural net
# fft = fast fourier transform
# Everything should be normalized by the max intensity

# def nn_dist_from_mass(mass: float, isolen: int, type: str):
#     """Mass ->  to Neural Net dist"""
#     isodist = (ctypes.c_float * 128)()
#     mass_c = ctypes.c_float(mass)
#     isolen = ctypes.c_int(isolen)
#     isogen_lib.fft_pep_mass_to_dist(mass_c, isodist, isolen, 0)
#     return np.array(isodist, dtype=np.float32)


def fft_rna_seq_to_dist(seq, isolen=128, offset=0):
    """RNA Seq -> FFT IsoJim dist"""
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * isolen)()
    isolen_c = ctypes.c_int(isolen)
    offset_c = ctypes.c_int(offset)
    result = isogen_lib.fft_rna_seq_to_dist(seq_bytes, isodist, isolen_c, offset_c)
    if result == -1:
        return None
    return np.array(isodist).copy()


def fft_rna_mass_to_dist(mass, isolen=128):
    """RNA Mass -> IsoDist"""
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    isolen_c = ctypes.c_int(isolen)
    offset = ctypes.c_int(0)
    result = isogen_lib.fft_rna_mass_to_dist(fmass, isodist, isolen_c, offset)
    if result == -1:
        return None
    return np.array(isodist[:])


def brain_rna_seq_to_dist(seq, isolen=128, offset=0):
    """RNA Seq -> FFT IsoJim dist"""
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * isolen)()
    isolen_c = ctypes.c_int(isolen)
    offset_c = ctypes.c_int(offset)
    isogen_lib.brain_rna_seq_to_dist(seq_bytes, isodist, isolen_c, offset_c)
    return np.array(isodist).copy()


def brain_rna_mass_to_dist(mass, isolen=128):
    """RNA Mass -> IsoDist"""
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    isolen_c = ctypes.c_int(isolen)
    offset = ctypes.c_int(0)
    isogen_lib.brain_rna_mass_to_dist(fmass, isodist, isolen_c, offset)

    return np.array(isodist).copy()


def nn_rna_seq_to_dist(seq, isolen=128, offset=0):
    """RNA sequence -> Isotope Distribution"""
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * isolen)()
    offset = ctypes.c_int(offset)
    isolen = ctypes.c_int(isolen)
    result = isogen_lib.nn_rna_seq_to_dist(seq_bytes, isodist, isolen, offset)
    if result == -1:
        return None
    return np.array(isodist).copy()


def nn_rna_mass_to_dist(mass, isolen=128, offset=0):
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    offset = ctypes.c_int(offset)
    isolen = ctypes.c_int(isolen)
    result = isogen_lib.nn_rna_mass_to_dist(fmass, isodist, isolen, offset)
    if result == -1:
        return None
    return np.array(isodist).copy()


def nn_pep_seq_to_dist(seq, isolen=128, offset=0):
    """Proteoform Sequence -> Isotope Distribution"""
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * isolen)()
    isolen = ctypes.c_int(isolen)
    offset = ctypes.c_int(offset)
    result = isogen_lib.nn_pep_seq_to_dist(seq_bytes, isodist, isolen, offset)
    if result == -1:
        return None
    return np.array(isodist).copy()


def fft_pep_seq_to_dist(seq, isolen=128, offset=0):
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * isolen)()
    isolen_c = ctypes.c_int(isolen)
    offset_c = ctypes.c_int(offset)
    result = isogen_lib.fft_pep_seq_to_dist(seq_bytes, isodist, isolen_c, offset_c)
    if result == -1:
        return None
    return np.array(isodist).copy()

def brain_pep_seq_to_dist(seq, isolen=128, offset=0):
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * isolen)()
    isolen_c = ctypes.c_int(isolen)
    offset_c = ctypes.c_int(offset)
    result = isogen_lib.brain_pep_seq_to_dist(seq_bytes, isodist, isolen_c, offset_c)
    return np.array(isodist).copy()


def fft_pep_mass_to_dist(mass, isolen=128, offset=0):
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    offset = ctypes.c_int(0)
    isolen = ctypes.c_int(isolen)
    result = isogen_lib.fft_pep_mass_to_dist(fmass, isodist, isolen, offset)
    if result == -1:
        return None
    return np.array(isodist).copy()

def brain_pep_mass_to_dist(mass, isolen=128, offset=0):
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    offset = ctypes.c_int(0)
    isolen = ctypes.c_int(isolen)
    result = isogen_lib.brain_pep_mass_to_dist(fmass, isodist, isolen, offset)
    return np.array(isodist).copy()

def nn_pep_mass_to_dist(mass, isolen=128, offset=0):
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    offset = ctypes.c_int(0)
    isolen = ctypes.c_int(isolen)
    result = isogen_lib.nn_pep_mass_to_dist(fmass, isodist, isolen, offset)
    if result == -1:
        return None
    return np.array(isodist).copy()


def plot_protein_dist_speed_differences(massList, isolen):
    # Only plotting the last entry of massList for each
    FFTTotal = 0
    for mass in massList:
        FFTStart = time.perf_counter()
        #protein_dist = fast_calc_averagine_isotope_dist(mass, isolen=isolen)[:,1]
        protein_dist = fft_pep_mass_to_dist(mass, isolen)
        FFTEnd = time.perf_counter()
        FFTTotal += FFTEnd - FFTStart
        NNTotal = 0

        NNStart = time.perf_counter()
        nn_dist = nn_pep_mass_to_dist(mass, isolen)
        NNEnd = time.perf_counter()
        NNTotal += NNEnd - NNStart
        mse = np.mean((protein_dist - nn_dist) ** 2)
        print(f"MSE: {mse}")
        print(f"Cosine Similarity: {1 - cosine(protein_dist, nn_dist)}")

        plt.plot(protein_dist, label=f"FFT toDist (Time: {FFTTotal:.4f}s)", color="r")
        plt.plot(nn_dist, label=f"Neural Net toDist (Time: {NNTotal:.4f}s)", color="b")
        plt.legend()
        plt.show()

def plot_rna_dist_differences(mass, isolen):
    mass_dist = fft_rna_mass_to_dist(mass, isolen)
    print("Calculated dist!")
    # nn_dist = nn_rna_dist(mass, isolen)
    # mse = np.mean((mass_dist - nn_dist) ** 2)
    # print(f"MSE: {mse}")

    # print(f"Cosine Similarity: {1 - cosine(mass_dist, nn_dist)}")
    plt.plot(mass_dist, label="FFT RNA Dist from Mass", color="red")
    #plt.plot(nn_dist, label="NN RNA Dist from Mass", color="b")
    plt.title(mass)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    import random
    #Test peptide methods
    if False:
        test_sequences = ["CCCCCCC"]
        test_masses = [1200, 11205, 38000]

        print("Testing Peptide Sequence to Dist Methods")
        for i in range(len(test_sequences)):
            nn_dist = nn_pep_seq_to_dist(test_sequences[i])
            fft_dist = fft_pep_seq_to_dist(test_sequences[i])
            brain_dist = brain_pep_seq_to_dist(test_sequences[i])

            plt.plot(nn_dist, color='r', label='NN')
            plt.plot(fft_dist, color='b', label='FFT')
            plt.plot(brain_dist, color='k', label='BRAIN')
            plt.legend()
            plt.show()

        print("Testing Peptide Mass to Dist Methods")
        for i in range(len(test_masses)):
            nn_dist = nn_pep_mass_to_dist(test_masses[i])
            fft_dist = fft_pep_mass_to_dist(test_masses[i])
            brain_dist = brain_pep_mass_to_dist(test_masses[i])

            plt.plot(nn_dist, color='r', label='NN')
            plt.plot(fft_dist, color='b', label='FFT')
            plt.plot(brain_dist, color='k', label='BRAIN')
            plt.legend()
            plt.show()


    #Test RNA methods
    if False:
        # test_sequences = ["GUAC", "GUACGUACGUACGUAC", "GUACGUACGUACGUACGUACGUACGUACGUACGUAC",
        #             "GUACGUACGUACGCGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUAC",
        #             "AAAAAAAAAAAAA", "UUUUUUUUUUUUU"]
        test_sequences = ["A"*20, "A"*25, "A"*40, "A"*50, "A"*100, "A"*150]
        test_masses = [1200, 11205, 98765]

        print("Testing RNA Sequence to Dist Methods")
        for i in range(len(test_sequences)):
            nn_seq_dist = nn_rna_seq_to_dist(test_sequences[i], isolen=64)
            fft_seq_dist = fft_rna_seq_to_dist(test_sequences[i], isolen=64)

            mass = igs.rnaseq_to_mass(test_sequences[i])
            nn_mass_dist = nn_rna_mass_to_dist(mass)
            fft_mass_dist = fft_rna_mass_to_dist(mass)

            nn_seq_css = calculate_cosinesimilarity(nn_seq_dist, fft_seq_dist)
            nn_mass_css = calculate_cosinesimilarity(nn_mass_dist, fft_seq_dist)

            print("NN-Seq CSS: ", nn_seq_css)
            print("NN-Mass CSS: ", nn_mass_css)


            plt.title(test_sequences[i])
            plt.plot(nn_seq_dist, color='r', label='NN-Seq')
            plt.plot(fft_seq_dist, color='b', label='FFT-Seq')
            plt.plot(nn_mass_dist, color='g', label='NN-Mass')
            plt.plot(fft_mass_dist, color='k', label='FFT-Mass')
            plt.legend()
            plt.show()

        print("Testing RNA Mass to Dist Methods")
        for i in range(len(test_masses)):
            nn_dist = nn_rna_mass_to_dist(test_masses[i])
            fft_dist = fft_rna_mass_to_dist(test_masses[i])
            brain_dist = brain_rna_mass_to_dist(test_masses[i])

            plt.plot(nn_dist, color='r', label='NN')
            plt.plot(fft_dist, color='b', label='FFT')
            plt.plot(brain_dist, color='k', label='BRAIN')
            plt.legend()
            plt.show()


    if False:
        from isogen_tools import *

        seq = "GUAC"
        rnamass = rnaseq_to_mass(seq)

        print(rnamass)

        fft_mass_dist = fft_rna_mass_to_dist(rnamass)
        fft_seq_dist = fft_rna_seq_to_dist(seq)


        from unidec.modules.isotopetools import isojim

        #Do the calculation as done by valkenborg to see if I can match their model
        #Order is C,H,N,O,S,P
        A = np.array([5,5,5,0,0,0])
        C = np.array([4,5,3,1,0,0])
        G = np.array([5,5,5,1,0,0])
        U = np.array([4,4,2,2,0,0])
        Ribose = np.array([5,10,0,5,0,0])
        Phosphate = np.array([0,3,0,4,0,1])
        water = np.array([0,2,0,1,0,0])

        valkenborg_isolist = np.array([0,0,0,0,0,0])
        #Add G
        valkenborg_isolist += G + Ribose - water + Phosphate - water

        #Add U
        valkenborg_isolist += U + Ribose - water + Phosphate - water

        #Add A
        valkenborg_isolist += A + Ribose - water + Phosphate - water

        #Add C
        valkenborg_isolist += C + Ribose - water + Phosphate - water

        valkenborg_formula = ""
        valkenborg_formula += "C" + str(valkenborg_isolist[0]) + "H" + str(valkenborg_isolist[1])
        valkenborg_formula += "N" + str(valkenborg_isolist[2]) + "O" + str(valkenborg_isolist[3]) + "P" + str(valkenborg_isolist[5])

        me_formula = rnaseq_to_formula("GUAC")
        me_dist = get_dist_from_formula(me_formula)
        me_dist = me_dist / max(me_dist)

        print("Valkenborg isolist:")
        for i,v in enumerate(valkenborg_isolist):
            print("isolist " + str(i) + ":" + str(v))

        valkenborg_dist = get_dist_from_formula(valkenborg_formula)
        valkenborg_dist = valkenborg_dist / max(valkenborg_dist)

        print("Valkenborg:", valkenborg_formula)
        print("Me:", me_formula)


        # plt.plot(fft_mass_dist, label="FFT-Mass")
        # plt.plot(fft_seq_dist, label="FFT-Seq")
        plt.plot(me_dist, label="Me")
        plt.plot(valkenborg_dist, label="Valkenborg")
        plt.legend()
        plt.show()


    if True:
        from isogen_tools import *
        from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

        lengths = np.arange(500, 1000, 100)

        a_seqs = np.array([])
        u_seqs = np.array([])

        a_fft_mass_css_vals = np.array([])
        a_nn_seq_css_vals = np.array([])

        for l in lengths:
            a_seqs = np.append(a_seqs, "A"*l)
            u_seqs = np.append(u_seqs, "U"*l)

        for i in range(len(lengths)):
            a_seq_mass = rnaseq_to_mass(a_seqs[i])
            a_seq_formula = rnaseq_to_formula(a_seqs[i])
            a_seq_formula_dist = get_dist_from_formula(a_seq_formula)
            a_seq_formula_dist = a_seq_formula_dist / max(a_seq_formula_dist)

            u_seq_mass = rnaseq_to_mass(u_seqs[i])
            u_seq_formula = rnaseq_to_formula(u_seqs[i])
            u_seq_formula_dist = get_dist_from_formula(u_seq_formula)
            u_seq_formula_dist = u_seq_formula_dist / max(u_seq_formula_dist)

            a_fft_mass_dist = fft_rna_mass_to_dist(a_seq_mass)
            a_fft_seq_dist = fft_rna_seq_to_dist(a_seqs[i])
            a_nn_seq_dist = nn_rna_seq_to_dist(a_seqs[i])

            u_fft_mass_dist = fft_rna_mass_to_dist(u_seq_mass)
            u_fft_seq_dist = fft_rna_seq_to_dist(u_seqs[i])
            u_nn_seq_dist = nn_rna_seq_to_dist(u_seqs[i])

            a_fft_mass_css = calculate_cosinesimilarity(a_fft_mass_dist, a_fft_seq_dist)
            a_fft_mass_css_vals = np.append(a_fft_mass_css_vals, a_fft_mass_css)

            a_nn_seq_css = calculate_cosinesimilarity(a_nn_seq_dist, a_fft_seq_dist)
            a_nn_seq_css_vals = np.append(a_nn_seq_css_vals, a_nn_seq_css)

            plt.plot(a_fft_seq_dist, label="A-FFT-Seq")
            plt.plot(a_seq_formula_dist, label="A-Formula")
            plt.plot(u_fft_seq_dist, label="U-FFT-Seq")
            plt.plot(u_seq_formula_dist, label="U-Formula")
            plt.xlim(0)
            plt.ylim(0)
            plt.title("length:" + str(len(u_seqs[i])))
            plt.legend()
            plt.show()

        plt.plot(a_fft_mass_css_vals, label="FFT-Mass")
        plt.plot(a_nn_seq_css_vals, label="NN-Seq")
        plt.legend()
        plt.show()

    #Doing speed tests for large molecules
    if True:
        mass = 100000

        timestart = time.perf_counter()
        fft_mass_dist = fft_pep_mass_to_dist(mass, isolen=1024)
        fft_mass_time = round(time.perf_counter() - timestart, 5)


        timestart = time.perf_counter()
        brain_mass_dist = brain_pep_mass_to_dist(mass, isolen=1024)
        brain_mass_time = round(time.perf_counter() - timestart, 5)

        print("FFT Time: ", fft_mass_time)
        print("BRAIN Time: ", brain_mass_time)

        plt.plot(fft_mass_dist, label='FFT-Mass')
        plt.plot(brain_mass_dist, label='BRAIN-Mass')
        plt.legend()
        plt.show()



