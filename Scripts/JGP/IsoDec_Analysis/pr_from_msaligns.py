from unidec.IsoDec.match import *

tof_msalign = "Z:\\Group Share\\JGP\\js8b05641_si_001\\FLASHDeconv_results\\tof_ecd_spectrum_ms1.msalign"
uvpd_msalign = "Z:\\Group Share\\JGP\\js8b05641_si_001\\FLASHDeconv_results\\uvpd_spectrum_ms1.msalign"
etd_msalign = "Z:\\Group Share\\JGP\\js8b05641_si_001\\FLASHDeconv_results\\1500_scans_200K_16 fills-qb1_ms2.msalign"

tof_annotated = "Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\apomyoglobin_tof_manual_annotations.txt"
uvpd_annotated = "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt"
etd_annotated = "Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt"

tof_mc = read_msalign_to_matchedcollection(tof_msalign)
uvpd_mc = read_msalign_to_matchedcollection(uvpd_msalign)
etd_mc = read_msalign_to_matchedcollection(etd_msalign)

tof_annotated = read_manual_annotations(tof_annotated, "\t")
uvpd_annotated = read_manual_annotations(uvpd_annotated, "\t")
etd_annotated = read_manual_annotations(etd_annotated, " ")

total_shared = 0
total_u1 = 0
total_u2 = 0

shared, u1, u2 = compare_annotated(tof_mc, tof_annotated,ppmtol=20, maxshift=3)
total_shared += len(shared)
total_u1 += len(u1)
total_u2 += len(u2)

shared, u1, u2 = compare_annotated(uvpd_mc, uvpd_annotated,ppmtol=20, maxshift=3)
total_shared += len(shared)
total_u1 += len(u1)
total_u2 += len(u2)

shared, u1, u2 = compare_annotated(etd_mc, etd_annotated,ppmtol=20, maxshift=3)
total_shared += len(shared)
total_u1 += len(u1)
total_u2 += len(u2)

precision = total_shared / (total_shared + total_u1)
recall = total_shared / (total_shared + total_u2)
f1 = 2 * (precision * recall) / (precision + recall)

print("Precision:", precision)
print("Recall:", recall)
print("F1 Score:", f1)