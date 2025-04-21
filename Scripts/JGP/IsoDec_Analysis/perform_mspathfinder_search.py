import os
import subprocess

os.chdir("C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247")

mspathfinder_exe_path = "C:\\Program Files\\Informed-Proteomics\\MSPathFinderT.exe"

datafile_path = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\20170105_L_MaD_ColC4_Ecoli20161108-10_01.raw"
database_path = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\uniprot-taxonomy3A83333.fasta"

ms1_msalign = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\topFD_msalign\\20170105_L_MaD_ColC4_Ecoli20161108-10_01_ms2.msalign"
ms2_msalign = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\topFD_msalign\\20170105_L_MaD_ColC4_Ecoli20161108-10_01_ms2.msalign"

command = [mspathfinder_exe_path, '-s',datafile_path, '-d',database_path]

result = subprocess.run(command, capture_output=True, text=True)