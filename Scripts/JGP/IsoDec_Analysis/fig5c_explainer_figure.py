import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from unidec.IsoDec.match import read_manual_annotations
import unidec.tools as ud
import os
from unidec.IsoDec.engine import *
import matplotlib.patches as patches
from matplotlib import rcParams

# Set backend to Agg
mpl.use('WxAgg')
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

def plot_predvec(predvec, highlighted_zs, round):
    plt.figure(figsize=(3, 4))
    rects = []
    for z in highlighted_zs:
        # Create a box around the max index
        rect = patches.Rectangle((-0.5, z - 0.5), 1, 1, linewidth=1.5, edgecolor='r', facecolor='none')
        rects.append(rect)

    plt.imshow(np.transpose([predvec]), aspect=0.1, cmap="viridis", interpolation="none")
    for rectangle in rects:
        plt.gca().add_patch(rectangle)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks([])
    plt.tight_layout()
    if False:
        plt.savefig("C:\\Users\\Admin\\University of Arizona\\martylab - "
                    "Documents\\Papers\\IsoDec\\Figures\\predvec_z"+str(highlighted_zs[0])+ "_round" + str(round) +".png",
                    dpi=300, transparent=True)
        plt.savefig("C:\\Users\\Admin\\University of Arizona\\martylab - "
                    "Documents\\Papers\\IsoDec\\Figures\\predvec_z"+str(highlighted_zs[0])+ "_round" + str(round) +".svg",
                    dpi=300, transparent=True, format="svg")
    plt.show()

file2 = "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"
os.chdir(os.path.dirname(file2))

spectrum = np.loadtxt("uvpd_spectrum.txt")

uvpd_manual_annotations = read_manual_annotations(
    "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt",
    delimiter='\t')

mz_range = [910,917]

if mz_range is not None:
    spectrum = ud.datachop(spectrum, mz_range[0], mz_range[1])
    uvpd_manual_annotations.peaks = [peak for peak in uvpd_manual_annotations.peaks if peak.mz > mz_range[0] and peak.mz < mz_range[1]]



eng = IsoDecEngine(phaseres=8, use_wrapper=False)
eng.config.verbose = True
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
eng.config.knockdown_rounds  = 1
eng.config.matchtol = 5


phasemodel = PhaseModel()
phasemodel.modelid = 1

manual_rounds = 20
for round in range(manual_rounds):
    if round> 0:
        eng.config.peakwindow = eng.config.peakwindow * 0.5
    eng.pks = MatchedCollection()

    max_intensity = np.max(spectrum[:,1])
    max_plotted = max_intensity * 1.2
    # Pick peaks
    peaks = fastpeakdetect(spectrum, window=int(eng.config.peakwindow), threshold=eng.config.peakthresh)

    encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(spectrum, peaks,
                                                                            lowmz=eng.config.mzwindow[
                                                                                0],
                                                                            highmz=eng.config.mzwindow[
                                                                                1],
                                                                            minpeaks=2,
                                                                            datathresh=eng.config.datathreshold)

    indices_to_remove = np.array([])

    rectangles = []
    for i in range(len(goodpeaks)):
        if goodpeaks[i][0] >913 and goodpeaks[i][0] <915:
            rect = patches.Rectangle((goodpeaks[i][0] - 1.05, -max_plotted), 3.10, 2 * max_plotted, color='b',
                                     alpha=0.25)
            rectangles.append(rect)


        predvec = phasemodel.predict_returnvec(encodingcentroids[i])[0]
        top_two_indices = predvec.argsort()[-2:][::-1]
        print("Top score: ", predvec[top_two_indices[0]], " Top 2: ", top_two_indices, predvec[top_two_indices[1]])
        report_two = predvec[top_two_indices[1]] / predvec[top_two_indices[0]] > eng.config.zscore_threshold
        if report_two:
            matchedindexes = eng.get_matches_multiple_z(outcentroids[i], top_two_indices, goodpeaks[i][0], pks=eng.pks)
        else:
            matchedindexes = eng.get_matches(outcentroids[i], int(top_two_indices[0]), goodpeaks[i][0], pks=eng.pks)

        if len(matchedindexes) > 0:
            indval = indexes[i]
            matchindvals = indval[matchedindexes]
            eng.pks.peaks[-1].matchedindexes = np.array(matchindvals)

            if report_two:
                plot_predvec(predvec, top_two_indices, round)
            else:
                plot_predvec(predvec, [int(top_two_indices[0])], round)
        else:
            kindex = fastnearest(spectrum[:, 0], goodpeaks[i][0])
            indices_to_remove = np.append(indices_to_remove, kindex)



    #Now remove the matchedcentroids from the original spectrum.
    indices_to_remove = np.array([])
    for i in range(len(eng.pks.peaks)):
        for index in eng.pks.peaks[i].matchedindexes:
            if index == 0:
                continue
            else:
                indices_to_remove = np.append(indices_to_remove, index)

    indices_to_remove = indices_to_remove.astype(int)
    removed_spectrum = spectrum[indices_to_remove]
    spectrum = np.delete(spectrum, indices_to_remove, axis=0)


    for i in range(len(spectrum)):
        plt.vlines(spectrum[i,0], 0, spectrum[i,1], color="k")

    for i in range(len(removed_spectrum)):
        plt.vlines(removed_spectrum[i,0], 0, removed_spectrum[i,1], color="r")

    for i in range(len(eng.pks.peaks)):
        for j in range(len(eng.pks.peaks[i].isodist)):
            plt.vlines(eng.pks.peaks[i].isodist[j,0], -eng.pks.peaks[i].isodist[j,1], 0, color="b")

    for rectangle in rectangles:
        plt.gca().add_patch(rectangle)

    plt.hlines(0, mz_range[0], mz_range[1], color="k")
    plt.xticks(np.arange(910,921,1), fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(910,917)
    plt.ylim(-max_plotted,max_plotted)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlabel("m/z", fontstyle='italic', fontsize=24)
    plt.ylabel("Intensity", fontsize=24)
    plt.tight_layout()
    plt.savefig("C:\\Users\\Admin\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\5c_explainer_iteration" +
                str(round) + ".png", dpi=300, transparent=True)
    plt.savefig("C:\\Users\\Admin\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\5c_explainer_iteration" +
                str(round) + ".pdf", dpi=300, transparent=True)
    plt.show()


