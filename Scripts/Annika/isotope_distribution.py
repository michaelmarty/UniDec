from unidec.IsoDec.IsoGen.isogen_tools import *

popc = "C42H82NO8P"

isodist = get_dist_from_formula(popc, 5)
print(isodist)

dpc = "C41H73 D7 NO8P"

d7pcdist = get_dist_from_formula(dpc, 5)
print(d7pcdist)