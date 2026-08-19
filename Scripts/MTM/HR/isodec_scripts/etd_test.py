from Scripts.MTM.proteomicstools import *
from unidec.IsoDec.compare_tools import read_manual_annotations, compare_annotated, read_msalign_to_matchedcollection, \
    compare_matchedcollections
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *
import pyperclip
from unidec.IsoDec.match import *
from UniDecImporter.ImporterFactory import ImporterFactory

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

file2 = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
os.chdir(os.path.dirname(file2))

mzrange = [1100, 1200]
mzrange = None

eng = IsoDecEngine(phaseres=8, use_wrapper=True, verbose=False)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2
eng.predmode = 0

eng_thrash = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False)
eng_thrash.config.activescan = 1
eng_thrash.config.activescanrt = -1
eng_thrash.config.activescanorder = 2
eng_thrash.config.knockdown_rounds = 5
eng_thrash.predmode = 0

if False:
    reader = ImporterFactory.create_importer(file2)
    spectrum = reader.grab_centroid_data(1)
    np.savetxt("etd_spectrum.txt", spectrum)

else:
    spectrum = np.loadtxt("etd_spectrum.txt")

# if mzrange is not None:
#     spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])
#     np.savetxt("chopped.txt", spectrum)

eng.batch_process_spectrum(spectrum, centroided=True)
eng_thrash.batch_process_spectrum(spectrum, centroided=True)
print("Peaks:", len(eng.pks.peaks))

exit()
if False:
    ca_seq = 'S[42.010565]HHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEATSRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLVHWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNVLDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRGFPK'
    ion_types = ['c', 'z']
    searchTD(eng.pks.peaks, ca_seq, ion_types, "ppm", 20.0, verbose=True)

if True:
    actual = read_manual_annotations("Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt", delimiter=' ', data=spectrum)

    shared, u1, u2 = compare_annotated(actual, eng.pks,  ppmtol=50, maxshift=3)
    print("IsoDec 8 Default Parameters Results")
    get_precision_recall(actual, eng.pks)

    shared_thrash, u1_thrash, u2_thrash = compare_annotated(actual, eng_thrash.pks,  ppmtol=50, maxshift=3)

    shared_eng, u_nn, u_thrash = compare_matchedcollections(shared, shared_thrash, ppmtol=50, maxshift=3)


    plot_pks(shared_eng, centroids=spectrum, ccolor="k", forcecolor="g", show=False, plotmass=False, nocentroids=True, labelz=True)
    plot_pks(u_nn, forcecolor="r", centroids=spectrum, ccolor="k", show=False, plotmass=False, nocentroids=True, labelz=True)
    plot_pks(u_thrash, forcecolor="b", show=False, plotmass=False, nocentroids=True, labelz=True)
    plt.show()

if True:
    print("TopFD Results:")
    etd_topfd = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\js8b05641_si_001\\MASH_Workflows\\1500_scans_200K_16 fills-qb1_20240815_120332_TopFD_1_1\\1500_scans_200K_16 fills-qb1_ms2.msalign", spectrum)
    searchTD(etd_topfd, ca_seq, ion_types, "ppm", 20.0, verbose=True)
    get_precision_recall(actual, etd_topfd)

    print("eTHRASH Results:")
    etd_eTHRASH = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\js8b05641_si_001\\MASH_Workflows\\1500_scans_200K_16 fills-qb1_20240816_104358_eTHRASH_1_1\\1500_scans_200K_16 fills-qb1.msalign", spectrum)
    searchTD(etd_eTHRASH, ca_seq, ion_types, "ppm", 20.0, verbose=True)
    get_precision_recall(actual, etd_eTHRASH)

    print("FLASHDeconv Results:")
    etd_flash = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\js8b05641_si_001\\FLASHDeconv_results\\1500_scans_200K_16 fills-qb1_ms2.msalign", spectrum)
    searchTD(etd_flash, ca_seq, ion_types, "ppm", 20.0, verbose=True)
    get_precision_recall(actual, etd_flash)

    # print("TopFD Comparison:")
    # id_shared, id_unique, tfd_shared, tfd_unique = compare_matched_ions(eng.pks, etd_topfd)
    # print("eTHRASH Comparison:")
    # compare_matched_ions(eng.pks, etd_eTHRASH)
    #
    # plot_pks(id_shared, forcecolor="g", show=False, plotmass=False, nocentroids=True, labelz=True)
    # plot_pks(id_unique, centroids=spectrum, ccolor="k", forcecolor="b", show=False, plotmass=False, nocentroids=True,
    #          labelz=True)
    # plot_pks(tfd_unique, forcecolor="r", show=False, plotmass=False, nocentroids=True, labelz=True)
    # plt.show()

if False:
    for peak in etd_topfd:
        eng.pks.add_peak(peak)
    for peak in etd_eTHRASH:
        eng.pks.add_peak(peak)

    for peak in eng.pks.peaks:
        peak.matchedions = []
    print("Combined coverage:")
    searchTD(eng.pks.peaks, ca_seq, ion_types, "ppm", 20.0, verbose=True)



