from unidec.IsoDec.compare_tools import read_msalign_to_matchedcollection

ms1_msalign_file = "D:\\Johnny\\IsoDec\\PXD019247\\topFD_msalign\\TopFD_20170105_L_MaD_ColC4_Ecoli20161108-10_01_ms1.msalign"
ms2_msalign_file = "D:\\Johnny\\IsoDec\\PXD019247\\topFD_msalign\\TopFD_20170105_L_MaD_ColC4_Ecoli20161108-10_01_ms2.msalign"

ms1_MC = read_msalign_to_matchedcollection(ms1_msalign_file, None)
ms2_MC = read_msalign_to_matchedcollection(ms2_msalign_file, None)

print("MS1 Peaks:", len(ms1_MC.peaks))
print("MS2 Peaks:", len(ms2_MC.peaks))