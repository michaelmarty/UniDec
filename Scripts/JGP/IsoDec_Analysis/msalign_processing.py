from unidec.IsoDec.match import *

ms1_msalign = read_msalign_to_matchedcollection("D:\\Johnny\\IsoDec\\PXD019247\\FLASHDeconv_msalign\\20170105_L_MaD_ColC4_Ecoli20161116-12_01_ms1.msalign")
ms2_msalign = read_msalign_to_matchedcollection("D:\\Johnny\\IsoDec\\PXD019247\\FLASHDeconv_msalign\\20170105_L_MaD_ColC4_Ecoli20161116-12_01_ms2.msalign")

print("MS1 Peaks: ", len(ms1_msalign.peaks))
print("MS2 Peaks: ", len(ms2_msalign.peaks))
