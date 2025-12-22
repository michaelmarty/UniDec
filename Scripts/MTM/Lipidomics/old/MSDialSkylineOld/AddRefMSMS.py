import os
from matchms.importing import load_from_msp
# from matchms.exporting import save_as_msp

from LipidFunctions import *


topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
os.chdir(topdir)

posfile = "PosIDsC1.csv"
posdf = pd.read_csv(posfile)
refposdf = parse_xml("Pos.xml")
posdf = add_ref_from_xml(posdf, refposdf)

negfile = "NegIDsC1.csv"
negdf = pd.read_csv(negfile)
refnegdf = parse_xml("Neg.xml")
negdf = add_ref_from_xml(negdf, refnegdf)

rttol = 0.1
posdf, negdf = merge_pos_neg(posdf, negdf, rttol)

# posdf.to_csv("PosIDsC1_Ref.csv", index=False)
# negdf.to_csv("NegIDsC1_Ref.csv", index=False)

poslibfile = "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp"
neglibfile = "MSDIAL-TandemMassSpectralAtlas-VS69-Neg-C1.msp"
# libfile = neglibfile
# outlibfile = libfile.replace(".msp", "_Ref.msp")
posdf, outdata = add_ref_spectra(posdf, poslibfile, rttol=3)
#
# # Save outdata to new msp file
# save_as_msp(outdata, outlibfile, mode="w")
# # Save dfpos to new csv file
posdf.to_csv(posfile.replace(".csv", "_Ref2.csv"), index=False)




