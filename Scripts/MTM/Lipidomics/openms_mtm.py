import sys
pathroot ="C:\\Python\\lipyd\\src"
sys.path.append(pathroot)

from pyopenms import *
import os
import unidec.tools as ud

def filter(mzarray, mztarget, mztolerance, rtarray, rttarget, rttolerance):
    b1 = np.abs(mzarray-mztarget)<mztolerance
    b2 = np.abs(rtarray-rttarget)<rttolerance
    b3 = b1 * b2
    return np.any(b3)

toppath = "C:\\Data\\Lipidomics\\OpenMS\\"
os.chdir(toppath)

mzML_files_top = ud.match_files(toppath, "*.mzml", exclude="aligned")

feature_maps=[]
mzML_files = [file[:-5]+"_aligned.mzML" for file in mzML_files_top]
mgf_files = [file[:-5]+"_aligned.mgf" for file in mzML_files_top]


consensusXML_file = "FeatureMatrix.consensusXML"

for f in mzML_files_top:
    feature_map = FeatureMap()
    FeatureXMLFile().load(f[:-4]+"featureXML", feature_map)
    feature_maps.append(feature_map)

consensus_map = ConsensusMap()
ConsensusXMLFile().load(consensusXML_file, consensus_map)
cdf = consensus_map.get_df()
mzs = cdf["mz"].to_numpy()
rts = cdf["RT"].to_numpy()

for file in mzML_files:
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    filtered = MSExperiment()
    for s in exp:
        if s.getMSLevel()==2:
            p = s.getPrecursors()[0]
            pmz = p.getMZ()
            p.setMZ(701)
            p.setCharge(-1)
            s.setPrecursors([p])
            prt = s.getRT()
            polarity = s.getInstrumentSettings().getPolarity()
            s.setName(str(pmz) + " " + str(prt))
            if filter(mzs, pmz, 0.01, rts, prt, 10) and pmz > 700 and pmz < 750:
                print(pmz, prt)
                filtered.addSpectrum(s)

    MzMLFile().store(file[:-5] + "_conM2.mzML", filtered)
    #MSPGenericFile().store(file[:-5] + "_conM2.msp", filtered)
    # MascotGenericFile().store(file[:-5] + "_conM2.mgf", filtered)
exit()
for f in consensus_map:
    mz = f.getMZ()
    rt = f.getRT()
    print(mz, rt)

