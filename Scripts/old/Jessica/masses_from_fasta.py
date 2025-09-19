import os
import numpy as np
from pyteomics import mass
from pyteomics import fasta
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/Scripts/old/Jessica")

file = "ecoli.fasta"
outfile = "ecoli_masses.csv"

mods=[-131.040485]
array = []
mvals=[]
with fasta.read(file,parser=fasta.std_parsers['uniprotkb']) as db:
    for descr, seq in db:
        #print(descr)
        try:
            mval = mass.calculate_mass(sequence=seq)
        except:
            mval = 0
        #print(mval)
        subarray=[]
        subarray.append(descr["entry"])
        subarray.append(mval)
        for m in mods:
            subarray.append(mval+m)
        #print(subarray)
        array.append(subarray)
        mvals.append(mval)

df = pd.DataFrame(array)
df.to_csv(outfile)

mvals=np.array(mvals)
mvals=np.sort(mvals)
diffs=np.diff(mvals)

print(np.average(diffs))
print(np.std(diffs))



plt.plot(diffs)
plt.show()

#np.savetxt("ecoli_masses.csv", array, delimiter=",")

