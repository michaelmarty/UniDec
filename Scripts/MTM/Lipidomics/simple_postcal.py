from mzml_postcal import *

file = "C:\Data\Lipidomics\\220808_SONAR_IMS_data\SONAR\\SWATH.mzML"
mzml_out = file[:-5] + "_fixed.mzML"
postcal(file, mzml_out, 0.03348)
