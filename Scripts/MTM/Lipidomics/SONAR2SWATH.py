from pyopenms import *
import os
import pandas as pd
import math

def read_sonar_cal(file):
    cals = []
    with open(file) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if "Calibration=" in line:
                l = line.strip()
                l = l.split("=")[1]
                l = l.split(",")
                cals.append(l)
                print(line.strip())
    return np.array(cals).astype(float)




sonarfile = "C:\Data\Lipidomics\\220808_SONAR_IMS_data\\220808_Ecoli_Lightsplash_SONAR_Res_Neg.raw\_sonar.inf"
cals = read_sonar_cal(sonarfile)
print(cals)

highcal = np.polyfit(cals[:,1], cals[:,0], 1)
midcal = np.polyfit(cals[:,2], cals[:,0], 1)
lowcal = np.polyfit(cals[:,3], cals[:,0], 1)

print(lowcal, midcal, highcal)

df = pd.DataFrame()
cycleid = np.arange(0, 201)
df["CycleID"]=cycleid
type = ["SWATH" for c in cycleid]
type[0]="SCAN"
df["Type"]= type
mzbegin = cycleid * lowcal[0] + lowcal[1]
mzbegin[0]=400
mzend = cycleid * highcal[0] + highcal[1]
mzend[0]=1600
df["m/z begin"] = mzbegin.astype(int)
df["m/z end"] = mzend.astype(int)
print(df)


file = "C:\Data\Lipidomics\\220808_SONAR_IMS_data\SONAR\\220808_Ecoli_Lightsplash_SONAR_Res_Neg_sep.mzML"
os.chdir(os.path.dirname(file))
df.to_csv("SWATH_exp.txt", sep="\t", index=False)
#file = "test.mzML"
exit()

exp = MSExperiment()
MzMLFile().load(file, exp)
filtered = MSExperiment()

merger = SpectraMerger()
params = merger.getDefaults()
# params.setValue("mz_binning_width", float(mz_binning_width))
# params.setValue("mz_binning_width_unit", b"ppm" if is_ppm else b"Da")
params.setValue("block_method:rt_block_size", 201)
# params.setValue("block_method:ms_levels", [ms_level])
merger.setParameters(params)

fswitch = False

for s in exp:
    # rt=s.getRT()
    # if rt < 30:
    #    filtered.addSpectrum(s)
    id = s.getNativeID()
    ids = id.split()
    d = dict()
    for thing in ids:
        thing = thing.split("=")
        d[thing[0]] = int(thing[1])
    func = d["function"]
    scan = d["scan"]
    s1 = scan - 1
    e = s1 % 200 + 1
    c = math.floor(scan/200)
    #print(d)
    if not fswitch and func == 1:
        fswitch = True
        ms1exp = MSExperiment()
    elif fswitch and func == 2:
        fswitch = False
        ms1exp.updateRanges()
        merger.mergeSpectraBlockWise(ms1exp)
        news = ms1exp[0]
        newid = "sample=1 period=1" + " cycle=" + str(c) + " experiment=1"
        news.setNativeID(newid)
        filtered.addSpectrum(news)
        print(d, newid)

    if func == 1:
        ms1exp.addSpectrum(s)
    elif func == 2:
        newid = "sample=1 period=1" + " cycle=" + str(c) + " experiment=" + str(e+1)
        s.setNativeID(newid)
        s.setType(1)
        filtered.addSpectrum(s)
        #print(d, newid)
        #exit()

MzMLFile().store("SWATH.mzML", filtered)
