import os
import numpy as np
import matplotlib.pyplot as plt
from unidec.IsoDec.trainingdata import *


os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

file = "C:\\Data\\TabbData\\20170105_L_MaD_ColC4_Ecoli20161108-10_01.raw"
file = "C:\\Data\\TabbData\\LPT_CYTO_GFPb.RAW"
#file = "C:\\Data\\TabbData\\20170105_L_MaD_ColC4_Ecoli20161116-12_01.raw"
#file = "D:\\Data\\TopDown\\denatured_discovery_input_topfd.mzxml"

zdat = process_file(file)
'''
# list files in directory
dir = "D:\\Data\\TopDown\\"
dir = "Z:\\Group Share\\JGP\\MSV000091923\\raw"
dir = "Z:\\Group Share\\JGP\\MSV000090488"
dir = "Z:\\Group Share\\JGP\\PXD019247"
# dir = "Z:\\Group Share\\JGP\\RibosomalPfms_Td_Control\\"

zdat = process_dir(dir)'''

plt.hist(zdat, bins=np.arange(0.5, 50.5, 1))
plt.show()
