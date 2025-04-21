import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import unidec.tools as ud
from unidec.IsoDec.match import *
from unidec.IsoDec.plots import *
from unidec.UniDecImporter.Importer import ImporterFactory

file2 = "Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\10uMCA25+UVPDLPT1,0mJ1pFullP_3.raw"
os.chdir(os.path.dirname(file2))

df1 = pd.read_csv("IsoDec.tsv", sep="\t")
df2 = pd.read_csv("ethrash60_v2.tsv", sep="\t")


mc1 = df_to_matchedcollection(df1)

mc2 = df_to_matchedcollection(df2)

sharedmf, u1mf, u2mf = get_unique_matchedions(mc1, mc2)
print("Shared:" + str(len(sharedmf)))
for p in sharedmf:
    print(p.mz, p.monoiso, p.matchedion)

print("Unique 1:" + str(len(u1mf)))
for p in u1mf:
    print(p.mz, p.monoiso, p.matchedion)

print("Unique 2:" + str(len(u2mf)))
for p in u2mf:
    print(p.mz, p.monoiso, p.matchedion)


shared, u1, u2 = compare_matchedcollections(mc1, mc2)

print(len(shared), len(u1), len(u2))

reader = ImporterFactory.create_importer(file2)
spectrum = reader.grab_centroid_data(1)


plot_pks(u1, forcecolor="r", show=False, plotmass=False, nocentroids=True)
plot_pks(u2, forcecolor="b", show=False, plotmass=False, nocentroids=True)
plot_pks(shared, centroids=spectrum, ccolor="k", forcecolor="g", show=False, plotmass=False, nocentroids=True)
plt.show()