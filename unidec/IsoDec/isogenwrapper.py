import ctypes
import os
import numpy as np
import time
import platform
from unidec.tools import start_at_iso

current_path = os.path.dirname(os.path.realpath(__file__))

if platform.system() == "Windows":
    dllname = "isogen.dll"
elif platform.system() == "Linux":
    dllname = "isogen.so"

else:
    print("Not yet implemented for MacOS")
    dllname = "isogen.dylib"


dllpath = start_at_iso(dllname, guess=current_path)

if not dllpath:
    print("DLL not found anywhere:", dllname)
# else:
#     print("DLL Found at:", dllpath)

isodist = ctypes.c_float * 64

isogen_c_lib = ctypes.CDLL(dllpath)

isogen_c_lib.fft_rna_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_c_lib.fft_rna_mass_to_dist.restype = ctypes.c_float

isogen_c_lib.nn_rna_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_int]
isogen_c_lib.nn_rna_mass_to_dist.restype = ctypes.c_float

isogen_c_lib.fft_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int,
                                            ctypes.c_int]
isogen_c_lib.fft_pep_mass_to_dist.restype = ctypes.c_float

isogen_c_lib.nn_pep_mass_to_dist.argtypes = [ctypes.c_float, ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.c_int]
isogen_c_lib.nn_pep_mass_to_dist.restype = ctypes.c_float


def nn_gen_isodist(mass, type="PEPTIDE", isolen=64):

    # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    res = b""
    if type is None:
        return
    if type == "RNA":
        res = b"Rna"
        # Call the C function
        isogen_c_lib.nn_rna_mass_to_dist(
            ctypes.c_float(mass),
            ptr,
            ctypes.c_int(isolen),
            ctypes.c_int(0)  # offset, not used in this case
        )
    elif type == "PEPTIDE":
        res = b"Peptide"
        # Call the C function
        isogen_c_lib.nn_pep_mass_to_dist(
            ctypes.c_float(mass),
            ptr,
            ctypes.c_int(isolen),
            ctypes.c_int(0)  # offset, not used in this case
        )

    # Convert isodist to numpy
    isodist = np.ctypeslib.as_array(isodist)
    return isodist

def fft_gen_isodist(mass, type="PEPTIDE", isolen = 128):
    # Create empty array
    isodist = np.zeros(isolen).astype(np.float32)
    ptr = isodist.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    if type == "RNA":
        result = isogen_c_lib.fft_rna_mass_to_dist(ctypes.c_float(mass), ptr, ctypes.c_int(isolen), ctypes.c_int(0))
        isodist = np.ctypeslib.as_array(isodist)

    elif type == "PEPTIDE":
        result = isogen_c_lib.fft_pep_mass_to_dist(ctypes.c_float(mass), ptr, ctypes.c_int(isolen), ctypes.c_int(0))
        isodist = np.ctypeslib.as_array(isodist)

    else:
        print("Unknown type for FFT generation:", type)
        return None

    return np.array(isodist)


if __name__ == "__main__":

    m =20000
    d1 = fft_gen_isodist(m, type="PEPTIDE")
    print(d1)

    # import matplotlib.pyplot as plt
    # plt.plot(d1/np.amax(d1), label="Isogen")
    # plt.legend()
    # plt.show()
    exit()
    #
    #
    # eng = IsoGenWrapper(dllpath=dllpath)
    # n = 10000
    # random_masses = np.random.uniform(1000, 60000, n)
    # starttime = time.perf_counter()
    # for mass in random_masses:
    #     dist = eng.gen_isodist(mass)
    #
    # print("Isogen Time:", time.perf_counter()-starttime)
    # print("Microseconds Per:", (time.perf_counter()-starttime)/n* 1e6)
    #
    # starttime = time.perf_counter()
    # for mass in random_masses:
    #     dist = eng.gen_isomike(mass)
    #
    # print("Isomike Time:", time.perf_counter()-starttime)
    # print("Microseconds Per:", (time.perf_counter()-starttime)/n* 1e6)
    #
    # starttime = time.perf_counter()
    # for mass in random_masses:
    #     dist = eng.gen_isofft(mass)
    #
    # print("FFT Time:", time.perf_counter()-starttime)
    # print("Microseconds Per:", (time.perf_counter()-starttime)/n * 1e6)
