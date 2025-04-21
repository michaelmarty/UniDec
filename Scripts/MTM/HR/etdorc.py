from sympy.integrals.intpoly import point_sort

from Scripts.MTM.proteomicstools import *
from unidec.IsoDec.engine import IsoDecEngine
import os
import time
import matplotlib.pyplot as plt
import unidec.tools as ud
from unidec.IsoDec.match import *
from unidec.UniDecImporter.ImporterFactory import ImporterFactory

file2 = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
os.chdir(os.path.dirname(file2))

mzrange = [700, 720]
mzrange = None

eng = IsoDecEngine(phaseres=4, use_wrapper=False, verbose=False)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
eng.predmode = 2

reader = ImporterFactory.create_importer(file2)
spectrum = reader.grab_centroid_data(1)
if mzrange is not None:
    spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])

timestart = time.time()


def main():
    actual_pks = read_manual_annotations(delimiter=' ')
    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    ppmtol = 50

    recalls = []
    precisions = []
    for threshold in thresholds:
        eng.config.css_thresh = threshold
        eng.batch_process_spectrum(spectrum, centroided=True)
        # Get the number of matches
        truepos, d, z = compare_annotated(actual_pks, eng.pks, 50, 10)
        print(len(truepos))
        truepos = len(truepos)
        falseneg = max(0, len(actual_pks) - truepos)

        falsepos = len(eng.pks) - truepos
        print(f"truepos : {truepos}")
        print(f"falseneg : {falseneg}")
        print(f"falsepos : {falsepos}")
        # Calculate recall and precision
        recall = truepos / (truepos + falseneg) if (truepos + falseneg) > 0 else 0

        precision = truepos / (truepos + falsepos) if (truepos + falsepos) > 0 else 0

        recalls.append(recall)
        precisions.append(precision)
        eng.pks = MatchedCollection()
        print(
            f"Threshold: {threshold}, Matches: {truepos}, True Positives: {truepos}, False Negatives: {falseneg}, False Positives: {falsepos}, Recall: {recall}, Precision: {precision}")

    print("Recalls:", recalls)
    print("Precisions:", precisions)

    plt.figure()
    plt.plot(recalls, precisions, "bo-", markersize=10, label="IsoDec")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Recall vs. Precision")
    plt.grid()
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend()
    plt.show()


main()
