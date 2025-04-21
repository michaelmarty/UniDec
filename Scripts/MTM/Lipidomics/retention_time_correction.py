import matplotlib.pyplot as plt
import pandas as pd
import unidec.tools as ud
from mzml_postcal import *

toppath = "C:\\Data\\Lipidomics\\OpenMS\\"
os.chdir(toppath)

df = pd.read_csv("Filtered.csv")
print(df.keys())

b1 = df.Ontology == "PG"
b2 = df.Ontology == "PE"
b3 = df.Ontology == "CL"
#df = df[b1+b2+b3]

refmz = df["Reference m/z"].to_numpy()
mz = df["Average Mz"].to_numpy()
rt = df["Average Rt(min)"].to_numpy()
refrt = df["Reference RT"].to_numpy()
ontcodes = df["Ontology"].astype("category").cat.codes


plt.scatter(rt, refrt, marker="o", c=ontcodes)
# plt.scatter(kmass[b2], kmd[b2], marker="o", c="red")

fit = np.polyfit(rt*60, refrt*60, 1)
fitdat = fit[0]*rt + fit[1]/60
print(fit)

plt.plot(rt, fitdat)


mzML_files = ud.match_files(toppath, "*.mzml", exclude=None)
print(mzML_files)

for file in mzML_files:
    mzml_in = file
    mzml_out = file[:-5] + "_rtfixed.mzML"
    postcal_RT(mzml_in, mzml_out, fit[0], fit[1])

plt.show()