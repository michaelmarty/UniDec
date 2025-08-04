import ctypes
import numpy as np
import time
from scipy.spatial.distance import cosine

from unidec.IsoDec.IsoGen import isogen_tools as igs
import matplotlib as mpl
import matplotlib.pyplot as plt
from unidec.IsoDec.IsoGen.isogen_tools import get_dist_from_formula
from unidec.modules.isotopetools import fast_calc_averagine_isotope_dist

mpl.use('WxAgg')
np.set_printoptions(precision=6, suppress=True)
isogen_dll_path = "C:/Python/UniDecDev/unidec/IsoDec/isogen.dll"

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

isogen_lib.nn_pep_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int]
isogen_lib.nn_pep_seq_to_dist.restype = ctypes.c_float

isogen_lib.fft_rna_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                           ctypes.c_int]
isogen_lib.fft_rna_seq_to_dist.restype = ctypes.c_float

isogen_lib.nn_rna_seq_to_dist.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_float), ctypes.c_int]
isogen_lib.nn_rna_seq_to_dist.restype = ctypes.c_float

isogen_lib.nn_rna_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_lib.nn_rna_mass_to_dist.restype =  ctypes.c_float



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
    isogen_lib.fft_rna_seq_to_dist(seq_bytes, isodist, isolen_c, offset_c)
    return np.array(isodist).copy()


def fft_rna_mass_to_dist(mass, isolen):
    """RNA Mass -> IsoDist"""
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    isolen_c = ctypes.c_int(isolen)
    offset = ctypes.c_int(0)

    isogen_lib.fft_rna_mass_to_dist(fmass, isodist, isolen_c, offset)

    return np.array(isodist[:])


def nn_rna_seq_to_dist(seq, offset=0):
    """RNA sequence -> Isotope Distribution"""
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * 64)()
    offset = ctypes.c_int(offset)
    isogen_lib.nn_rna_seq_to_dist(seq_bytes, isodist, offset)
    return np.array(isodist).copy()


def nn_rna_mass_to_dist(mass, offset=0):
    isolen = igs.rnamass_to_isolen(mass)

    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    offset = ctypes.c_int(0)
    isolen = ctypes.c_int(isolen)
    isogen_lib.nn_rna_mass_to_dist(fmass, isodist, isolen, offset)
    return np.array(isodist).copy()


def nn_pep_seq_to_dist(seq, offset=0):
    """Proteoform Sequence -> Isotope Distribution"""
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * 64)()
    offset = ctypes.c_int(offset)
    isogen_lib.nn_pep_seq_to_dist(seq_bytes, isodist, offset)
    return np.array(isodist).copy()


def fft_pep_seq_to_dist(seq, isolen=64, offset=0):
    seq_bytes = seq.encode("utf-8")
    isodist = (ctypes.c_float * isolen)()
    isolen_c = ctypes.c_int(isolen)
    offset_c = ctypes.c_int(offset)
    result = isogen_lib.fft_pep_seq_to_dist(seq_bytes, isodist, isolen_c, offset_c)
    return np.array(isodist).copy()


def fft_pep_mass_to_dist(mass, isolen, offset=0):
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    offset = ctypes.c_int(0)
    isolen = ctypes.c_int(isolen)
    isogen_lib.fft_pep_mass_to_dist(fmass, isodist, isolen, offset)
    # move to the c code
    return np.array(isodist).copy()


def nn_pep_mass_to_dist(mass, offset=0):
    isolen = igs.pepmass_to_isolen(mass)
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * isolen)()
    offset = ctypes.c_int(0)
    isolen = ctypes.c_int(isolen)
    isogen_lib.nn_pep_mass_to_dist(fmass, isodist, isolen, offset)
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
        test_sequences = ["AGCKLIHRK"]
        test_masses = [1200, 11205, 98765]

        print("Testing Peptide Sequence to Dist Methods")
        for i in range(len(test_sequences)):
            nn_dist = nn_pep_seq_to_dist(test_sequences[i])
            fft_dist = fft_pep_seq_to_dist(test_sequences[i], 128)

            plt.plot(nn_dist, color='r', label='NN')
            plt.plot(fft_dist, color='b', linestyle='--', label='FFT')
            plt.legend()
            plt.show()

        for i in range(len(test_masses)):
            nn_dist = nn_pep_mass_to_dist(test_masses[i])
            fft_dist = fft_pep_mass_to_dist(test_masses[i], 128)

            plt.plot(nn_dist, color='r', label='NN')
            plt.plot(fft_dist, color='b', linestyle='--', label='FFT')
            plt.legend()
            plt.show()


    #Test RNA methods
    if True:
        test_sequences = ["GUAC", "GUACGUACGUACGUAC", "GUACGUACGUACGUACGUACGUACGUACGUACGUAC",
                    "GUACGUACGUACGCGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUAC",
                    "AAAAAAAAAAAAA", "UUUUUUUUUUUUU"]
        test_masses = [1200, 11205, 98765]

        print("Testing RNA Sequence to Dist Methods")
        for i in range(len(test_sequences)):
            nn_dist = nn_rna_seq_to_dist(test_sequences[i])
            fft_dist = fft_rna_seq_to_dist(test_sequences[i], 128)

            plt.plot(nn_dist, color='r', label='NN')
            plt.plot(fft_dist, color='b', linestyle='--', label='FFT')
            plt.legend()
            plt.show()

        print("Testing RNA Mass to Dist Methods")
        for i in range(len(test_masses)):
            nn_dist = nn_rna_mass_to_dist(test_masses[i])
            fft_dist = fft_rna_mass_to_dist(test_masses[i], 128)

            plt.plot(nn_dist, color='r', label='NN')
            plt.plot(fft_dist, color='b', linestyle='--', label='FFT')
            plt.legend()
            plt.show()
