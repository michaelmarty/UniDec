from unidec.IsoDec.compare_tools import read_msalign_to_matchedcollection, peak_mz_z_df_to_matchedcollection, \
    compare_matchedcollections
from unidec.IsoDec.engine import IsoDecEngine
import pyperclip
import pandas as pd

from UniDecImporter.Importer import ImporterFactory


def process_file(file, comp_MC, remove_noise=False):
    eng = IsoDecEngine(phaseres=8)
    eng.use_wrapper = False

    eng.config.activescan = 1
    eng.config.activescanrt = -1
    eng.config.activescanorder = 2

    reader = ImporterFactory.create_importer(file)
    if file.endswith(".raw"):
        spectrum = reader.grab_centroid_data(1)
    elif file.endswith(".mzML"):
        spectrum = reader.grab_scan_data(0)

    if remove_noise:
        spectrum = remove_noise_cdata(spectrum, 100, factor=1.5, mode="median")

    eng.batch_process_spectrum(spectrum, centroided=True, window=80)

    u1, u2, shared = compare_matchedcollections(eng.pks, comp_MC, ppmtol=50, objecttocompare="monoiso", maxshift=3)

    return u1, u2, shared



write_individual_results = True
write_combined_results = True

etd_file = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
etd_df = pd.read_csv("Z:\\Group Share\\JGP\\js8b05641_si_001\\CAII_All_Ions.txt", sep="\t")
etd_mc = peak_mz_z_df_to_matchedcollection(etd_df, None)
etd_shared, etd_u1, etd_u2  = process_file(etd_file, etd_mc)

uvpd_file = "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"
uvpd_mc = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\MASH_Workflows\\10uMCA25+UVPDLPT1,0mJ1pFullP_3_20240813_102854_eTHRASH_1_1\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.msalign", None)
uvpd_shared, uvpd_u1, uvpd_u2 = process_file(uvpd_file, uvpd_mc)

aTpm_CID = "Z:\\Group Share\\JGP\\PXD018043\\aTpm_CID_20180612_Rh2426_F2-6_CID_21V_1259.5mz_3width_2M_30%_10.0s_62scans_000001.mzML"
aTpm_CID_mc = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\PXD018043\\expert_annotations_2020_04_29\\expert_annotations_2020_04_29\\aTpm\\CID\\Expert Annotations\\aTpm_CID_Ground_Truth.msalign", None)
aTpm_CID_shared, aTpm_CID_u1, aTpm_CID_u2 = process_file(aTpm_CID, aTpm_CID_mc)

aTpm_ECD = "Z:\\Group Share\\JGP\\PXD018043\\aTpm_ECD_20170616_F1-8_042917_ECD_1364.4mz_3width_0.7V_20ms_1.0s_1000scans_000001.mzML"
aTpm_ECD_mc = read_msalign_to_matchedcollection("Z:\\Group Share\\JGP\\PXD018043\\expert_annotations_2020_04_29\\expert_annotations_2020_04_29\\aTpm\\ECD\\Expert Annotations\\aTpm_ECD_Ground_Truth.msalign", None)
aTpm_ECD_shared, aTpm_ECD_u1, aTpm_ECD_u2 = process_file(aTpm_ECD, aTpm_ECD_mc)


if write_individual_results:
    print("UVPD Peaks:", str(len(uvpd_u1) + len(uvpd_shared)))
    uvpd_perc_annotated_peaks_matched = (len(uvpd_mc.peaks) - len(uvpd_u2)) / len(uvpd_mc.peaks) * 100
    print("UVPD Percent Annotated Peaks Matched:", uvpd_perc_annotated_peaks_matched)
    uvpd_perc_clusters_matched = len(uvpd_shared) / (len(uvpd_u1) + len(uvpd_shared)) * 100
    print("UVPD Percent Clusters Matched:", uvpd_perc_clusters_matched)

    print("ETD Peaks:", str(len(etd_u1) + len(etd_shared)))
    etd_perc_annotated_peaks_matched = (len(etd_mc.peaks) - len(etd_u2))  / len(etd_mc.peaks) * 100
    print("ETD Percent Annotated Peaks Matched:", etd_perc_annotated_peaks_matched)
    etd_perc_clusters_matched = len(etd_shared) / (len(etd_u1) + len(etd_shared)) * 100
    print("ETD Percent Clusters Matched:", etd_perc_clusters_matched)

    print("aTpm CID Peaks:", str(len(aTpm_CID_u1) + len(aTpm_CID_shared)))
    aTpm_CID_perc_annotated_peaks_matched = (len(aTpm_CID_mc.peaks) - len(aTpm_CID_u2)) / len(aTpm_CID_mc.peaks) * 100
    print("aTpm CID Percent Annotated Peaks Matched:", aTpm_CID_perc_annotated_peaks_matched)
    aTpm_CID_perc_clusters_matched = len(aTpm_CID_shared) / (len(aTpm_CID_u1) + len(aTpm_CID_shared)) * 100
    print("aTpm CID Percent Clusters Matched:", aTpm_CID_perc_clusters_matched)

    print("aTpm ECD Peaks:", str(len(aTpm_ECD_u1) + len(aTpm_ECD_shared)))
    aTpm_ECD_perc_annotated_peaks_matched = (len(aTpm_ECD_mc.peaks) - len(aTpm_ECD_u2)) / len(aTpm_ECD_mc.peaks) * 100
    print("aTpm ECD Percent Annotated Peaks Matched:", aTpm_ECD_perc_annotated_peaks_matched)
    aTpm_ECD_perc_clusters_matched = len(aTpm_ECD_shared) / (len(aTpm_ECD_u1) + len(aTpm_ECD_shared)) * 100
    print("aTpm ECD Percent Clusters Matched:", aTpm_ECD_perc_clusters_matched)


if write_combined_results:
    combined_shared = len(etd_shared) + len(uvpd_shared) + len(aTpm_CID_shared) + len(aTpm_ECD_shared)
    combined_u1 = len(etd_u1) + len(uvpd_u1) + len(aTpm_CID_u1) + len(aTpm_ECD_u1)
    combined_u2 = len(etd_u2) + len(uvpd_u2) + len(aTpm_CID_u2) + len(aTpm_ECD_u2)

    print("\n\nCombined Peaks:", str(combined_u1 + combined_shared))
    combined_perc_annotated_peaks_matched = (((len(uvpd_mc.peaks) - len(uvpd_u2)) + (len(etd_mc.peaks) - len(etd_u2)) +
                                             (len(aTpm_CID_mc.peaks) - len(aTpm_CID_u2)) + (len(aTpm_ECD_mc.peaks) - len(aTpm_ECD_u2)))
                                             / (len(etd_mc.peaks) + len(uvpd_mc.peaks) + len(aTpm_CID_mc.peaks) + len(aTpm_ECD_mc.peaks)) * 100)

    print("Combined Percent Annotated Peaks Matched:", combined_perc_annotated_peaks_matched)
    combined_perc_clusters_matched = combined_shared / (combined_u1 + combined_shared) * 100
    print("Combined Percent Clusters Matched:", combined_perc_clusters_matched)
