import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib import rcParams
from matplotlib import patches
from unidec.IsoDec.datatools import *
from unidec.IsoDec.encoding import *
from unidec.modules.unidecstructure import IsoDecConfig
from unidec.IsoDec.models import PhaseModel
from unidec.IsoDec.altdecon import *


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

    charges = range(len(predvec))
    plt.plot(charges, predvec, color="k")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlim(0, 30)
    plt.xlabel("Charge", fontsize=24)
    plt.ylabel("Score", fontsize=24)
    plt.xticks(np.arange(0, 31, 5), fontsize=18)
    plt.yticks(np.arange(-6,10,3), fontsize=18)
    plt.ylim(-6,9)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\Figures\\nn_thrash_m_um_nn_zmap.png",
                dpi=300, transparent=True)
    plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - "
                "Documents\\Papers\\IsoDec\\Figures\\nn_thrash_m_um_nn_zmap.pdf",
                dpi=300, transparent=True)
    plt.show()

# Set backend to Agg
mpl.use('WxAgg')
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

config = IsoDecConfig()
config.verbose = True
config.activescan = 1
config.activescanrt = -1
config.activescanorder = 2

file2 = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
os.chdir(os.path.dirname(file2))

spectrum = np.loadtxt("etd_spectrum.txt")
spectrum = remove_noise_cdata(spectrum, 100, factor=1.5, mode="median")

mz_range = [1167, 1172]

spectrum = ud.datachop(spectrum, mz_range[0], mz_range[1])

phasemodel = PhaseModel()
phasemodel.modelid = 1

peaks = fastpeakdetect(spectrum, window=int(config.peakwindow), threshold=config.peakthresh)

encodingcentroids, goodpeaks, outcentroids, indexes = extract_centroids(spectrum, peaks,
                                                                            lowmz=config.mzwindow[
                                                                                0],
                                                                            highmz=config.mzwindow[
                                                                                1],
                                                                            minpeaks=2,
                                                                            datathresh=config.datathreshold)

print("Total centroids:", len(encodingcentroids))
for i in range(len(encodingcentroids)):
    print("Peak:", goodpeaks[i][0])
    predvec = phasemodel.predict_returnvec(encodingcentroids[i])[0]
    predz = int(predvec.argmax())
    print("NN z:",predz)
    for j in range(len(outcentroids[i])):
        plt.vlines(outcentroids[i][j, 0], 0, outcentroids[i][j, 1], color="k")
    plt.show()
    plot_predvec(predvec, [predz], i)

    pk = optimize_shift2(config, outcentroids[i], predz, goodpeaks[i][0])
    if pk is not None:
        for j in range(len(outcentroids[i])):
            plt.vlines(outcentroids[i][j, 0], 0, outcentroids[i][j, 1], color="k")
        for j in range(len(pk.isodist)):
            plt.vlines(pk.isodist[j,0], -pk.isodist[j,1], 0, color="b")

        plt.hlines(0, mz_range[0], mz_range[1], color="k")
        plt.xticks(np.arange(mz_range[0],mz_range[1]+1,1), fontsize=18)
        plt.xlim(mz_range[0], mz_range[1])
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xlabel("m/z", fontstyle='italic', fontsize=24)
        plt.ylabel("Intensity", fontsize=24)
        plt.ylim(-10000,10000)
        plt.yticks(np.arange(-10000,10001, 5000), fontsize=18)
        plt.tight_layout()
        # plt.savefig("C:\\Users\\Admin\\University of Arizona\\martylab - "
        #             "Documents\\Papers\\IsoDec\\Figures\\nn_thrash_m_um_spectrum.png",
        #             dpi=300, transparent=True)
        # plt.savefig("C:\\Users\\Admin\\University of Arizona\\martylab - "
        #             "Documents\\Papers\\IsoDec\\Figures\\nn_thrash_m_um_spectrum.pdf",
        #             dpi=300, transparent=True)
        plt.show()



    linear_data, ac, fft, mul, fftx, acx = gen_thrash_arrays(encodingcentroids[i])
    thrash_predz = np.argmax(mul)
    charge = round(acx[thrash_predz])
    print("THRASH z:", charge)
    plt.plot(acx, mul, color="k")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.xlim(0,30)
    plt.xlabel("Charge", fontsize=24)
    plt.ylabel("Score", fontsize=24)
    plt.xticks(np.arange(0,31,5), fontsize=18)
    plt.yticks(np.arange(0,0.21,0.05), fontsize=18)
    plt.ylim(0,0.2)
    plt.tight_layout()
    plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - "
                "Documents\\Papers\\IsoDec\\Figures\\nn_thrash_m_um_fp_zmap.png",
                dpi=300, transparent=True)
    plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - "
                "Documents\\Papers\\IsoDec\\Figures\\nn_thrash_m_um_fp_zmap.pdf",
                dpi=300, transparent=True)
    plt.show()

