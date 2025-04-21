import pyopenms

from unidec.IsoDec.engine import IsoDecEngine
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud
from unidec.IsoDec.plots import plot_pks, gen_plot_fig
from unidec.IsoDec.datatools import *
from unidec.IsoDec.match import *
import pyperclip
from Scripts.MTM.proteomicstools import *
from unidec.UniDecImporter.ImporterFactory import ImporterFactory


def get_precision_recall(annotated_mc, exp_mc):
    shared, u_annotated, u_exp = compare_annotated(annotated_mc, exp_mc, ppmtol=50, maxshift=3)

    precision = len(shared) / (len(shared) + len(u_exp))
    recall = len(shared) / (len(shared) + len(u_annotated))
    f1 = 2 * precision * recall / (precision + recall)

    true_positives = len(shared)
    false_positives = len(u_exp)
    false_negatives = len(u_annotated)
    print("TP:", true_positives)
    print("FP:", false_positives)
    print("FN:", false_negatives)

    print("Precision:", precision)
    print("Recall:", recall)
    print("F1:", f1)

# Set backend to Agg
mpl.use('WxAgg')

file2 = "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"
os.chdir(os.path.dirname(file2))


eng = IsoDecEngine(phaseres=8, use_wrapper=False)

eng.config.verbose = False
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
# eng.config.knockdown_rounds = 5
# eng.config.css_thresh = 0.5
# eng.config.minareacovered = 0.1

eng.predmode = 0

uvpd_manual_annotations = read_manual_annotations(
    "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt",
    delimiter='\t')



if False:
    reader = ImporterFactory.create_importer(file2)
    spectrum = reader.grab_centroid_data(1)
    np.savetxt("uvpd_spectrum.txt", spectrum)

else:
    spectrum = np.loadtxt("uvpd_spectrum.txt")

    # mz = spectrum[:,0]
    # intensities = spectrum[:,1]
    #
    # exp = pyopenms.MSExperiment()
    #
    # spec = pyopenms.MSSpectrum()
    # spec.set_peaks((mz, intensities))
    # spec.setType(pyopenms.MSSpectrum.MSCentroid)
    # exp.addSpectrum(spec)
    #
    # pyopenms.MzMLFile().store("uvpd_spectrum_centroided.mzML", exp)
    # exit()


mzrange = [910,920]
mzrange = None

if mzrange is not None:
    spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])
    uvpd_manual_annotations.peaks = [peak for peak in uvpd_manual_annotations.peaks if mzrange[0] <= peak.mz <= mzrange[1]]

eng.batch_process_spectrum(spectrum, centroided=True, window=80)



if True:
    ca_seq = 'S[42.010565]HHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEATSRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLVHWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNVLDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRGFPK'
    ion_types = ['a','b','c','x','y','z','a+1','x+1','y-1']
    searchTD(eng.pks.peaks, ca_seq, ion_types, "ppm", 20.0, verbose=True)

#
# print("Written", len(eng.pks.peaks), "peaks to pstest.txt")
# print("Accepted clusters: " + str(len(eng.pks.peaks)))


if False:
    uvpd_mc = read_manual_annotations("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt", '\t', data=spectrum)
    uvpd_topfd = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\MASH_Workflows\\10uMCA25+UVPDLPT1,0mJ1pFullP_3_20240813_142813_TopFD_1_1\\input_topfd_ms2.msalign", spectrum)
    uvpd_flash = read_fd_tsv_to_matchedcollection("Z:\\Group Share\\JGP\\js8b05641_si_001\\FLASHDeconv_results\\uvpd_spectrum_centroided_ms1.tsv", spectrum)


    if mzrange is not None:
        eng.pks.peaks = [peak for peak in eng.pks.peaks if mzrange[0] <= peak.mz <= mzrange[1]]
        uvpd_mc.peaks = [peak for peak in uvpd_mc.peaks if mzrange[0] <= peak.mz <= mzrange[1]]

    print("IsoDec...")
    eng_shared, unique_annotated, unique_eng = compare_annotated(uvpd_mc, eng.pks, ppmtol=50, maxshift=3)
    precision = len(eng_shared) / (len(eng_shared) + len(unique_eng))
    recall = len(eng_shared) / (len(eng_shared) + len(unique_annotated))
    print("Precision:", precision)
    print("Recall:", recall)

    print("TopFD...")
    tfd_shared, unique_annotated, unique_tfd = compare_annotated(uvpd_mc, uvpd_topfd, ppmtol=50, maxshift=3)
    precision = len(tfd_shared) / (len(tfd_shared) + len(unique_tfd))
    recall = len(tfd_shared) / (len(tfd_shared) + len(unique_annotated))
    print("Precision:", precision)
    print("Recall:", recall)

    print("FLASHDeconv...")
    fld_shared, unique_annotated, unique_fld = compare_annotated(uvpd_mc, uvpd_flash, ppmtol=50, maxshift=3)
    precision = len(fld_shared) / (len(fld_shared) + len(unique_fld))
    recall = len(fld_shared) / (len(fld_shared) + len(unique_annotated))
    print("Precision:", precision)
    print("Recall:", recall)




    # plot_pks(unique_eng, ccolor='k', forcecolor="b", show=False, plotmass=False, nocentroids=True, labelz=True)
    plot_pks(unique_tfd, forcecolor="r", show=False, plotmass=False, nocentroids=True, labelz=True, labelpeaks=True)
    # plot_pks(unique_fld, forcecolor="r", show=False, plotmass=False, nocentroids=True, labelz=True, labelpeaks=True)
    plot_pks(uvpd_mc, centroids=spectrum, ccolor='k', forcecolor="g", show=False, plotmass=False, nocentroids=True, labelz=True, labelpeaks=True)
    plt.show()
    exit()

if False:
    spectrum = ud.datachop(spectrum, 913.25, 915.5)
    gen_plot_fig(unique_annotated, centroids=spectrum)
    plt.show()
    exit()


if True:
    print("TopFD Results:")
    uvpd_topfd = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\MASH_Workflows\\10uMCA25+UVPDLPT1,0mJ1pFullP_3_20240813_142813_TopFD_1_1\\input_topfd_ms2.msalign", spectrum)
    searchTD(uvpd_topfd, ca_seq, ion_types, "ppm", 20.0, verbose=True)
    get_precision_recall(uvpd_manual_annotations, uvpd_topfd)

    print("eTHRASH Results:")
    uvpd_eTHRASH = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\MASH_Workflows\\10uMCA25+UVPDLPT1,0mJ1pFullP_3_20240813_102854_eTHRASH_1_1\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.msalign", spectrum)
    searchTD(uvpd_eTHRASH, ca_seq, ion_types, "ppm", 20.0, verbose=True)
    get_precision_recall(uvpd_manual_annotations, uvpd_eTHRASH)

    print("FLASHDeconv Results:")
    uvpd_flash = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\js8b05641_si_001\\FLASHDeconv_results\\uvpd_spectrum_ms1.msalign", spectrum)
    searchTD(uvpd_flash, ca_seq, ion_types, "ppm", 20.0, verbose=True)
    get_precision_recall(uvpd_manual_annotations, uvpd_flash)



