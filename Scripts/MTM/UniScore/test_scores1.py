import os
from unidec import engine
import matplotlib.pyplot as plt
from unidec.modules.isolated_packages import score_window as sw

directory = "C:\\Data\\AlgorithmTesting"
#directory = "C:\\Python\\UniDec3\\bin\\Example Data"
os.chdir(directory)

#file = "ADH.txt"
#file = "GroEL HCD unidec.txt"

save = False
runall = True


if True:
    file = "250313_AQPZ_POPC_100.txt"
    eng = engine.UniDec()
    eng.open_file(file, directory)
    eng.unidec_imports(everything=True)
    eng.pick_peaks()
    eng.dscore()
    sw.score_plots2(eng, save=save)

if runall:
    file = "BSA.txt"
    eng = engine.UniDec()
    eng.open_file(file, directory)
    eng.unidec_imports(everything=True)
    eng.pick_peaks()
    eng.dscore()
    sw.score_plots2(eng, save=save)

if runall:
    file = "0.txt"
    eng = engine.UniDec()
    eng.open_file(file, directory)
    eng.unidec_imports(everything=True)
    eng.pick_peaks()
    eng.dscore()
    sw.score_plots2(eng, save=save)


if runall:
    file = "POPC_Nanodiscs.txt"
    eng = engine.UniDec()
    eng.open_file(file, directory)
    eng.unidec_imports(everything=True)
    eng.pick_peaks()
    eng.dscore()
    sw.score_plots2(eng, save=save)

if runall:
    file = "myo.txt"
    eng = engine.UniDec()
    eng.open_file(file, directory)
    eng.unidec_imports(everything=True)
    eng.pick_peaks()
    eng.dscore()
    sw.score_plots2(eng, save=save)

if False:
    file = "ADH.txt"
    eng = engine.UniDec()
    eng.open_file(file, directory)
    eng.unidec_imports(everything=True)
    eng.pick_peaks()
    eng.dscore()
    sw.score_plots2(eng, save=save)

    eng.estimate_FDR()

plt.show()