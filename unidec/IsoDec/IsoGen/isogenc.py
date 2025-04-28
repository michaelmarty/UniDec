import ctypes
import numpy as np
import time
from scipy.spatial.distance import cosine

from unidec.IsoDec.IsoGen import isogen_tools as igs
import matplotlib as mpl
import matplotlib.pyplot as plt

from unidec.IsoDec.IsoGen.isogen_tools import get_dist_from_formula

mpl.use('WxAgg')
np.set_printoptions(precision=6, suppress=True)
isogen_dll_path = "C:/Python/UniDec3/unidec/IsoDec/isogen.dll"

isogen_lib = ctypes.CDLL(isogen_dll_path)
isogen_lib.fft_pep_mass_to_dist.restype = ctypes.c_float
isogen_lib.fft_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]

isogen_lib.mass_to_formula_averaging.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_int)]
isogen_lib.mass_to_formula_averaging.restype = ctypes.POINTER(ctypes.c_int)

isogen_lib.fft_rna_seq_to_dist.argtypes = [ctypes.c_char_p]
isogen_lib.fft_rna_seq_to_dist.restype = ctypes.POINTER(ctypes.c_float)

isogen_lib.rnaToVector.argtypes = [ctypes.c_char_p]
isogen_lib.rnaToVector.restype = ctypes.POINTER(ctypes.c_float)

isogen_lib.rnaVectorToMass.argtypes = [ctypes.POINTER(ctypes.c_float)]
isogen_lib.rnaVectorToMass.restype = ctypes.c_float

isogen_lib.fft_rna_mass_to_dist.argtypes = [ctypes.c_float]
isogen_lib.fft_rna_mass_to_dist.restype = ctypes.POINTER(ctypes.c_float)

isogen_lib.isogen_nn.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_char_p]
isogen_lib.isogen_nn.restype = None


# nn = neural net
# fft = fast fourier transform
# Everything should be normalized by the max intensity

def nn_dist_from_mass(mass: float, isolen: int, type: str):
    """Mass ->  to Neural Net dist"""
    isodist = (ctypes.c_float * 128)()
    mass_c = ctypes.c_float(mass)
    isolen = ctypes.c_int(isolen)
    isogen_lib.isogen_nn(mass_c, isodist, isolen, b"PROTEIN")
    return np.array(isodist, dtype=np.float32)


def fft_rna_seq_to_dist(seq):
    """RNA Seq -> FFT IsoJim dist"""
    seq_bytes = seq.encode("utf-8")
    dist_ptr = isogen_lib.fft_rna_seq_to_dist(ctypes.c_char_p(seq_bytes))
    dist = np.ctypeslib.as_array(dist_ptr, shape=(128,))
    return np.array(dist.copy())


def fft_rna_mass_to_dist(mass):
    """RNA Mass -> IsoDist"""
    mass_c = ctypes.c_float(mass)
    isolist_ptr = isogen_lib.fft_rna_mass_to_dist(mass_c)
    isolist = np.ctypeslib.as_array(isolist_ptr, shape=(128,))
    return isolist


def rna_to_vector(seq):
    """RNA Seq -> IsoList"""
    seq_bytes = seq.encode("utf-8")
    isolist_ptr = isogen_lib.rnaToVector(ctypes.c_char_p(seq_bytes))
    isolist = np.ctypeslib.as_array(isolist_ptr, shape=(4,))
    return isolist


def fft_get_protein_vector(mass):
    """Protein Mass -> IsoList"""
    sequence = (ctypes.c_int * 5)()
    mass_c = ctypes.c_float(mass)
    isogen_lib.mass_to_formula_averaging(mass_c, sequence)
    res = np.ctypeslib.as_array(sequence)
    return res


def rna_vec_to_mass(rna_vec):
    """RNA vector-> Mass"""
    rna_vec_np = np.array(rna_vec, dtype=np.float32)
    ptr = rna_vec_np.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    result = isogen_lib.rnaVectorToMass(ptr)
    return result


def fft_pep_mass_to_dist(mass):
    fmass = ctypes.c_float(mass)
    isodist = (ctypes.c_float * 128)()
    offset = ctypes.c_int(0)
    isolen = ctypes.c_int(128)
    isogen_lib.fft_pep_mass_to_dist(fmass, isodist, isolen, offset)
    # move to the c code
    return np.array(isodist).copy()


def nn_rna_dist(mass, isolen):
    isodist = (ctypes.c_float * 128)()
    mass_c = ctypes.c_float(mass)
    isogen_lib.isogen_nn(mass_c, isodist, isolen, b"RNA")
    return np.array(isodist, dtype=np.float32)


def fft_grab_atom_isodist(sequence):
    seq = ctypes.c_char_p(sequence.encode("utf-8"))
    isodist = (ctypes.c_float * 128)()
    isogen_lib.fft_atom_sequence_to_dist(seq, isodist, 128)
    return np.array(isodist).copy()


def nn_atom_dist(sequence, isolen):
    isodist = (ctypes.c_float * 128)()
    sequence = sequence.encode("utf-8")
    isogen_lib.isogen_atom(ctypes.c_char_p(sequence), isodist, isolen)
    return np.array(isodist, dtype=np.float32)


def plot_protein_dist_speed_differences(massList, isolen):
    # Only plotting the last entry of massList for each
    FFTTotal = 0
    for mass in massList:
        FFTStart = time.perf_counter()
        protein_dist = fft_pep_mass_to_dist(mass)
        FFTEnd = time.perf_counter()
        FFTTotal += FFTEnd - FFTStart
        NNTotal = 0

        NNStart = time.perf_counter()
        nn_dist = nn_dist_from_mass(mass, isolen, b"PEPTIDE")
        NNEnd = time.perf_counter()
        NNTotal += NNEnd - NNStart
        mse = np.mean((protein_dist - nn_dist) ** 2)
        print(f"MSE: {mse}")
        print(f"Cosine Similarity: {1 - cosine(protein_dist, nn_dist)}")

        plt.plot(protein_dist, label=f"FFT toDist (Time: {FFTTotal:.4f}s)", color="r")
        plt.plot(nn_dist, label=f"Neural Net toDist (Time: {NNTotal:.4f}s)", color="b")
        plt.legend()
        plt.show()


def plot_rna_dist_differences(seq, isolen):
    mass = rna_vec_to_mass(rna_to_vector(seq))
    mass_dist = fft_rna_mass_to_dist(mass)

    nn_dist = nn_rna_dist(mass, isolen)
    mse = np.mean((mass_dist - nn_dist) ** 2)
    print(f"MSE: {mse}")
    print(f"Cosine Similarity: {1 - cosine(mass_dist, nn_dist)}")
    plt.plot(mass_dist, label="FFT RNA Dist from Mass", color="red")
    plt.plot(nn_dist, label="NN RNA Dist from Mass", color="b")
    plt.title(f"RNA Sequence: {seq} with isolen {isolen}")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # Protein comparison
    if True:
        import random
        isolen = 8
        if isolen == 8:
            massList = [random.uniform(10, 2000) for i in range(100)]
            plot_protein_dist_speed_differences(massList, 8)
        elif isolen == 32:
            massList = [random.uniform(1000, 10000) for i in range(100)]
            plot_protein_dist_speed_differences(massList, 32)
        elif isolen == 64:
            massList = [random.uniform(1000, 10000) for i in range(100)]
            plot_protein_dist_speed_differences(massList, 64)


    # RNA comparison
    if False:
        rna_sequences = [
            "GAUC", "UAGC", "CGAU", "GAUCA", "UGACU", "CGAUGC", "GAUCUG", "UAUGCAU",
            "GUACUGA", "GAUCGUA", "CGAUGCAU", "GAUCGACG", "UACGUGAU", "GAUCUGACU", "UGACGUAGC"]
        for i in rna_sequences:
            plot_rna_dist_differences(i, 32)

    # Atom comparison
    if True:
        isolen = 32
        testformulas = [
            "C6H12O6", "C2H6O", "C3H6O", "CH3COOH", "C7H6O2", "C3H7NO2", "C5H9NO4", "C6H14N4O2", "C5H11NO2", "CHCl3",
            "CHBr3", "CF4", "CCl4", "C2F6", "Pb(NO3)2", "UO2(NO3)2", "Bi2S3", "WCl6", "Fe2O3", "NaCl",
            "H2SO4", "HNO3", "C60", "C70", "C8H10N4O2", "C9H8O4", "C27H46O", "C12H22O11",
            "C10H16N5O13P3",
            "C55H72MgN4O5", "O2", "N2", "CO2", "H2O"
        ]

        for formula in testformulas:

            nn_dist = nn_atom_dist(formula, isolen)
            fft_dist = fft_grab_atom_isodist(formula)

            # print(fft_dist)

            plt.plot(fft_dist, label="FFT Atom Dist", color="red")
            plt.plot(nn_dist, label="NN Atom Dist", color="blue")
            plt.title(f"Sequence: {formula} with isolen {isolen}")
            plt.legend()
            # plt.plot(truedist, label="True Dist", color="green")
            plt.show()
