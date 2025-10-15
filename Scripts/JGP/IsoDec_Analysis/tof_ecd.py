from Scripts.MTM.proteomicstools import *
from unidec.IsoDec.compare_tools import read_manual_annotations, compare_annotated, read_msalign_to_matchedcollection
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.datatools import *
import pyperclip
from unidec.IsoDec.match import *
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
from unidec.UniDecImporter.MZML.mzML import *
from pyteomics import mgf




# def check_cal_drift(scans):
#     reference = scans[0]
#     ref_peaks = reference.peaks[:10]
#     for i in range(1, len(scans)):
#         #Compare the top 10 most intense peaks in each scan to the reference scan
#         scan_pks = scans[i].peaks[:10]
#         #Now see if any of the scan peaks are within 0.1 Da of the reference peaks
#         for ref_pk in ref_peaks:
#             for scan_pk in scan_pks:
#                 if abs(ref_pk.mz - scan_pk.mz) < 0.1:
#                     print("Matched Peak found, mz drift: ", ref_pk.mz - scan_pk.mz)

# Set backend to Agg
mpl.use('WxAgg')

file = "Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002.d"
os.chdir(os.path.dirname(file))


eng = IsoDecEngine(phaseres=4, use_wrapper=False, verbose=False)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

eng.config.knockdown_rounds = 5
eng.config.matchtol = 10
eng.config.minareacovered = 0.25

eng.predmode = 2

spectrum = None

if True:
    merge_bins  = 0.01 #Th
    merge_type = "Integrate"

    fwhm = 0.05 #Th
    centroiding_window = int((fwhm/merge_bins) / 2)

    centroiding_window = 5
    if True:
        spectrum = np.loadtxt("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum_full.txt")
    else:
        reader = ImporterFactory.create_importer(file)
        info = reader.msrun.scan_info()
        ms2_scan_indices = []
        ms2_scan_indices = [item[2] for i, item in enumerate(info) if item[3] == "MS2" and item[2] not in ms2_scan_indices]
        print("Loading data...")
        data = reader.get_all_scans()
        ms2_data = [data[i] for i in ms2_scan_indices]
        print("Merging data...")
        spectrum = merge_spectra(ms2_data, mzbins=merge_bins, type=merge_type)

    raw_spectrum = np.copy(spectrum)
    # np.savetxt("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum_full2.txt", raw_spectrum)
    spectrum = get_all_centroids(raw_spectrum, window=centroiding_window, threshold=eng.config.peakthresh, background=100, moving_average_smoothing=4)
    # plt.plot(raw_spectrum[:, 0], raw_spectrum[:, 1])
    # for c in spectrum:
    #     plt.plot([c[0], c[0]], [0, c[1]], color="k")
    # plt.show()
    # np.savetxt("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum.txt", spectrum)

    if False:
        spectra = [{'id': 'spectrum1',
                    'msLevel': 2,
                    'centroided': 1,
                    'm/z array': spectrum[:, 0],
                    'intensity array': spectrum[:, 1],
                    'params': {'id': 'spectrum1',
                               'centroided': 1,
                               'msLevel': 2,}}]

        mgf.write(spectra, "ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002.mgf")
        np.savetxt("ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum.txt", spectrum)
        exit()



if spectrum is None:
    spectrum = np.loadtxt("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum.txt")


# mzrange = [845, 855]
mzrange = None

if mzrange is not None:
    spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])

timestart = time.time()
eng.batch_process_spectrum(spectrum, centroided=True)
if True:
    apo_seq = "GLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELGFQG"
    searchTD(eng.pks.peaks, apo_seq, ["b", "y", "c", "z"], "ppm", 20.0, verbose=True)
    exit()


    #get peaks where there length of the peaks matchedions list is not 0
    matchedpeaks = [peak for peak in eng.pks.peaks if len(peak.matchedions) != 0]
    matchedpeaks_collection = MatchedCollection()
    matchedpeaks_collection.add_peaks(matchedpeaks)
    unmatchedpeaks = [peak for peak in eng.pks.peaks if len(peak.matchedions) == 0]
    unmatchedpeaks_collection = MatchedCollection()
    unmatchedpeaks_collection.add_peaks(unmatchedpeaks)

    topFD_tof_ecd = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\MASH_Workflows\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_20250102_120613_TopFD_1_1\\input_topfd_ms2.msalign")
    print("Top FD Results:")
    searchTD(topFD_tof_ecd, apo_seq, ["b", "y", "c", "z"], "ppm", 20.0, verbose=True)


    print("FLASHDeconv Results:")
    fd_tof_ecd = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\js8b05641_si_001\\FLASHDeconv_results\\tof_ecd_spectrum_ms1.msalign")
    searchTD(fd_tof_ecd, apo_seq, ["b", "y", "c", "z"], "ppm", 50.0, verbose=True)

# top_fd = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\MSV000094542\\raw\\RAW\\MASH_Workflows\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_20241113_112215_TopFD_1_1\\input_topfd_ms2.msalign",
#                                            data=spectrum)
#
# shared,u_id, u_tfd  = compare_matchedcollections(eng.pks, top_fd, ppmtol=50)

manual_annotations = read_manual_annotations("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\apomyoglobin_tof_manual_annotations.txt",
                                             delimiter='\t')


shared, unique_annotated, unique_experimental = compare_annotated(manual_annotations, eng.pks, ppmtol=50, maxshift=3)

print("Recall: ", len(shared.peaks) / len(manual_annotations))
print("Precision: ", len(shared.peaks) / len(eng.pks.peaks))

print("FLASHDeconv Results:")
shared, unique_annotated, unique_experimental = compare_annotated(manual_annotations, fd_tof_ecd, ppmtol=50, maxshift=3)
print("Recall: ", len(shared.peaks) / len(manual_annotations))
print("Precision: ", len(shared.peaks) / len(fd_tof_ecd.peaks))

print("TopFD Results:")
shared, unique_annotated, unique_experimental = compare_annotated(manual_annotations, topFD_tof_ecd, ppmtol=50, maxshift=3)
print("Recall: ", len(shared.peaks) / len(manual_annotations))
print("Precision: ", len(shared.peaks) / len(topFD_tof_ecd.peaks))
