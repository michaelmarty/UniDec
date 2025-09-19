from unidec.modules.thermo_reader.RawFileReader import RawFileReader as rr
from unidec.modules.thermo_reader.ThermoImporter import ThermoDataImporter as TI
#matplotlib.use("TkAgg")

path = "Z:\\Group Share\\Marius Kostelic\\Tiny tips and teflon\\09272021 BSA test\\202100924_MMK_BSA_normal_needle_new_3_3min_aqu.RAW"
#path = "C:\\Data\\223321.RAW"

ti = TI(path)

#TI.get_tic(ti)

rr.GetChromatogram(ti,scanrange=None, trace_number=0)
print("chromatogram:",ti)