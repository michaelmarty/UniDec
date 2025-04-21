from unidec.IsoDec.engine import IsoDecEngine
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import pandas as pd
from unidec.IsoDec.match import *

# Set backend to Agg
mpl.use('WxAgg')

etd_file = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
uvpd_file = "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"

etd_reader = ImporterFactory.create_importer(etd_file)
etd_spectrum = etd_reader.grab_centroid_data(1)

uvpd_reader = ImporterFactory.create_importer(uvpd_file)
uvpd_spectrum = uvpd_reader.grab_centroid_data(1)

etd_annotated = read_manual_annotations("Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt", delimiter=' ')
uvpd_annotated = read_manual_annotations("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt", delimiter='\t')

eng8 = IsoDecEngine(phaseres=8, verbose=False, use_wrapper=False)
eng8.config.activescan = 1
eng8.config.activescanrt = -1
eng8.config.activescanorder = 2

recalls8 = []
precisions8 = []

eng4 = IsoDecEngine(phaseres=4, verbose=False, use_wrapper=False)
eng4.config.activescan = 1
eng4.config.activescanrt = -1
eng4.config.activescanorder = 2

recalls4 = []
precisions4 = []

css_areacovered_kdrounds = \
    [(0.5, 0.2, 10),
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
for css, ac, kd in css_areacovered_kdrounds:
    print("CSS: " + str(round(css, 2)) + ", Minareacovered: " + str(round(ac, 2)) + ", " + str(kd) + " KD Rounds")
    truepos = 0
    falseneg = 0
    falsepos = 0

    eng8.config.css_thresh = css
    eng8.config.minareacovered = ac
    eng8.config.knockdown_rounds = kd

    eng8.batch_process_spectrum(etd_spectrum, centroided=True)

    shared, u1, u2 = compare_annotated(etd_annotated, eng8.pks, ppmtol=50, maxshift=3)

    eng8.pks = MatchedCollection()

    truepos += len(shared)
    falseneg += len(u1)
    falsepos += len(u2)

    eng8.batch_process_spectrum(uvpd_spectrum, centroided=True)
    shared, u1, u2 = compare_annotated(uvpd_annotated, eng8.pks, ppmtol=50, maxshift=3)

    truepos += len(shared)
    falseneg += len(u1)
    falsepos += len(u2)

    eng8.pks = MatchedCollection()

    recall = truepos / (truepos + falseneg)
    print("8bit Recall: " + str(recall))
    precision = truepos / (truepos + falsepos)
    print("8bit Precision: " + str(precision))
    F1_score = 2 * (precision * recall) / (precision + recall)
    print("8-bit F1 Score: " + str(F1_score))

    recalls8.append(recall)
    precisions8.append(precision)

    truepos = 0
    falseneg = 0
    falsepos = 0

    eng4.config.css_thresh = css
    eng4.config.minareacovered = ac
    eng4.config.knockdown_rounds = kd

    eng4.batch_process_spectrum(etd_spectrum, centroided=True)

    shared, u1, u2 = compare_annotated(etd_annotated, eng4.pks, ppmtol=50, maxshift=3)

    eng4.pks = MatchedCollection()

    truepos += len(shared)
    falseneg += len(u1)
    falsepos += len(u2)

    eng4.batch_process_spectrum(uvpd_spectrum, centroided=True)
    shared, u1, u2 = compare_annotated(uvpd_annotated, eng4.pks, ppmtol=50, maxshift=3)

    truepos += len(shared)
    falseneg += len(u1)
    falsepos += len(u2)

    eng4.pks = MatchedCollection()

    recall = truepos / (truepos + falseneg)
    print("4bit Recall: " + str(recall))
    precision = truepos / (truepos + falsepos)
    print("4bit Precision: " + str(precision))
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
plt.savefig("Z:\\\Group Share\\\JGP\\\IsoDec_figs\\uvpd_etd_roc.png", dpi=300, transparent=True)
plt.show()
plt.figure()

exit()

plot_pks(u1, forcecolor="r", show=False, plotmass=False, nocentroids=True)
plot_pks(u2, forcecolor="b", show=False, plotmass=False, nocentroids=True)
plot_pks(shared, centroids=spectrum, ccolor="k", forcecolor="g", show=False, plotmass=False, nocentroids=True)
plt.show()
