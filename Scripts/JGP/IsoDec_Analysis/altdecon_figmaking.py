import os
from UniDecImporter.ImporterFactory import ImporterFactory
import numpy as np
import pyopenms
import matplotlib.pyplot as plt
import matplotlib as mpl
from unidec.IsoDec.datatools import *
import unidec.tools as ud
from unidec.IsoDec.altdecon import *
from unidec.IsoDec.engine import IsoDecEngine
import matplotlib.patches as patches

# Set backend to Agg
mpl.use('WxAgg')

if False:
    file2 = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
    os.chdir(os.path.dirname(file2))
    reader = ImporterFactory.create_importer(file2)
    spectrum = reader.grab_centroid_data(1)


    mz_range = [701,705]
    spectrum = ud.datachop(spectrum, mz_range[0], mz_range[1])
    np.savetxt("Z:\\Group Share\\JGP\\js8b05641_si_001\\etd_chopped.txt", spectrum)
    exit()

spectrum = np.loadtxt("Z:\\Group Share\\JGP\\js8b05641_si_001\\etd_chopped.txt")
#spectrum = remove_noise_cdata(spectrum)
#normalize spectrum intensities to 100
spectrum[:,1] = spectrum[:,1]/np.max(spectrum[:,1])*100

if False:
    plt.figure(figsize=(3,4))
    for i in range(len(spectrum)):
        plt.plot([spectrum[i, 0], spectrum[i, 0]], [0, spectrum[i, 1]], color="black")

    plt.xticks(np.arange(702,705,1), fontsize=18)
    plt.yticks(np.arange(0,101,50), labels=["", "%", "100"], fontsize=18)
    plt.xlabel("m/z", fontsize=24, fontstyle='italic')
    plt.ylim(0,100)
    plt.xlim(702,704)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\etd_subspec.png", dpi=300)
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\etd_subspec.pdf", dpi=300)
    plt.show()
    exit()

if False:
    eng = IsoDecEngine(phaseres=8, use_wrapper=True, verbose=False)
    eng.config.activescan = 1
    eng.config.activescanrt = -1
    eng.config.activescanorder = 2

    eng.batch_process_spectrum(spectrum)
    for p in eng.pks.peaks:
        print(p.z, p.monoiso)


linear_data, ac, fft, mul, fftx, acx = gen_thrash_arrays(spectrum)

#Plot the fourier-patterson charge map
if False:
    plt.figure(figsize=(3, 4))
    plt.plot(acx, fft, color='k')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks(np.arange(0,51,10), fontsize=14)
    plt.xlabel("Charge", fontsize=18)
    plt.title("Fourier Charge Map", fontsize=24)
    plt.tight_layout()
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\fourier_zmap.png", dpi=300)
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\fourier_zmap.pdf", dpi=300)
    plt.show()

if False:
    plt.figure(figsize=(3, 4))
    plt.plot(acx, ac, color='k')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks(np.arange(0, 51, 10), fontsize=14)
    plt.xlabel("Charge", fontsize=18)
    plt.title("Patterson Charge Map", fontsize=24)
    plt.tight_layout()
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\patterson_zmap.png", dpi=300)
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\patterson_zmap.pdf", dpi=300)
    plt.show()

if False:
    plt.figure(figsize=(3, 4))

    #Create a rectangle between 5 and 12 on the x axis and from 0 to the maximum on the y
    rect = patches.Rectangle((5, 0), 5, np.max(mul)+0.05*np.max(mul), color='b', alpha=0.25)

    plt.plot(acx, mul, color='k')
    plt.gca().add_patch(rect)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks(np.arange(0, 51, 10), fontsize=14)
    plt.xlabel("Charge", fontsize=18)
    plt.tight_layout()
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\combination_zmap.png", dpi=300)
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\combination_zmap.pdf", dpi=300)
    plt.show()

if True:
    plt.figure(figsize=(3, 4))
    plt.plot(acx, mul, color='k')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.yticks([])
    plt.xticks(np.arange(5, 11, 1), fontsize=14)
    plt.xlim(5,10)
    plt.tight_layout()
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\zoomed_combination_zmap.png", dpi=300)
    plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\zoomed_combination_zmap.pdf", dpi=300)
    plt.show()

exit()