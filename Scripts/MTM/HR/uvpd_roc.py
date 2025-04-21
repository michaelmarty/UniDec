from unidec.IsoDec.engine import IsoDecEngine
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *
import pandas as pd
from unidec.IsoDec.match import *
from unidec.UniDecImporter.ImporterFactory import ImporterFactory

# Set backend to Agg
mpl.use('WxAgg')

#file2 = "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"
file2 = "D:\\Johnny\\IsoDec\\UVPD_Test\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"
os.chdir(os.path.dirname(file2))


eng8 = IsoDecEngine(phaseres=8, use_wrapper=False)
eng8.config.activescan = 1
eng8.config.activescanrt = -1
eng8.config.activescanorder = 2
eng8.predmode = 0

eng4 = IsoDecEngine(phaseres=4, use_wrapper=False)
eng4.config.activescan = 1
eng4.config.activescanrt = -1
eng4.config.activescanorder = 2
eng4.predmode = 2

reader = ImporterFactory.create_importer(file2)
spectrum = reader.grab_centroid_data(1)
empty = MatchedCollection()

annotated_peaks = read_manual_annotations("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt", delimiter='\t')

uvpd_eTHRASH = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\MASH_Workflows\\10uMCA25+UVPDLPT1,0mJ1pFullP_3_20240813_102854_eTHRASH_1_1\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.msalign", spectrum)

shared, u1, u2 = compare_annotated(annotated_peaks, uvpd_eTHRASH, ppmtol=50, maxshift=3)

truepos = len(shared)
falseneg = len(u1)
falsepos = len(u2)

recall = truepos / (truepos + falseneg)
print("eTHRASH Recall: " + str(recall))
precision = truepos / (truepos + falsepos)
print("eTHRASH Precision: " + str(precision))

recalls8 = []
precisions8 = []

recalls4 = []
precisions4 = []

css_areacovered_kdrounds = \
    [(0.1, 0.01, 15),
     (0.5, 0.2, 15),
     (0.6, 0.2, 13),
     (0.7, 0.2, 10),
     (0.7, 0.2, 9),
     (0.7, 0.2, 8),
     (0.7, 0.2, 7),
     (0.7, 0.2, 6),
     (0.7, 0.2, 5),
     (0.75, 0.2, 5),
     (0.8, 0.2, 5),
     (0.85, 0.2, 5),
     (0.9, 0.2, 5),
     (0.95, 0.2, 5),
     (0.98, 0.2, 5),
     (0.99, 0.2, 5)]



timestart = time.time()
for set in css_areacovered_kdrounds:
    eng8.config.css_thresh = set[0]
    eng8.config.minareacovered = set[1]
    eng8.config.knockdown_rounds = set[2]
    eng8.batch_process_spectrum(spectrum, centroided=True)

    shared8, unique_annotated8, unique_experimental8 = compare_annotated(annotated_peaks, eng8.pks, ppmtol=50, maxshift=3)

    eng8.pks = MatchedCollection()
    truepos = len(shared8)
    falseneg = len(unique_annotated8)
    falsepos = len(unique_experimental8)
    recall = truepos / (truepos + falseneg)
    print("8-bit Recall: " + str(recall))
    precision = truepos / (truepos + falsepos)
    print("8-bit Precision: " + str(precision))
    F1_score = 2 * (precision * recall) / (precision + recall)
    print("8-bit F1 Score: " + str(F1_score))

    recalls8.append(recall)
    precisions8.append(precision)

    eng4.config.css_thresh = set[0]
    eng4.config.minareacovered = set[1]
    eng4.config.knockdown_rounds = set[2]
    eng4.batch_process_spectrum(spectrum, centroided=True)

    shared4, unique_annotated4, unique_experimental4 = compare_annotated(annotated_peaks, eng4.pks, ppmtol=50,
                                                                         maxshift=3)

    eng4.pks = MatchedCollection()
    truepos = len(shared4)
    falseneg = len(unique_annotated4)
    falsepos = len(unique_experimental4)
    recall = truepos / (truepos + falseneg)
    print("4-bit Recall: " + str(recall))
    precision = truepos / (truepos + falsepos)
    print("4-bit Precision: " + str(precision))
    F1_score = 2 * (precision * recall) / (precision + recall)
    print("4-bit F1 Score: " + str(F1_score))

    recalls4.append(recall)
    precisions4.append(precision)

print("Time:", time.time() - timestart)

plt.figure()
# Create ROC plot
plt.plot(recalls8, precisions8, "bo", markersize=10, label="IsoDec-8bit")
plt.plot(recalls4, precisions4, "ro", markersize=10, label="IsoDec-4bit")
plt.legend()
plt.xlabel("Recall", fontsize=24, fontweight='bold')
plt.ylabel("Precision", fontsize=24, fontweight='bold')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xticks(np.arange(0,1.1,0.5), fontsize=16)
plt.yticks(np.arange(0,1.1,0.5), fontsize=16)
plt.tight_layout()
plt.savefig("Z:\\\Group Share\\\JGP\\\IsoDec_figs\\uvpd_roc.png", dpi=300, transparent=True)
plt.show()
plt.figure()

exit()
