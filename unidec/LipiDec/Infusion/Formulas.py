import numpy as np
import pandas as pd
from molmass import Formula

def get_count(comp, element):
    try:
        return comp[element].count
    except:
        return 0

# Read the lib file
libfile = "C:\\Data\\Lipidomics\\Libraries\\nonder1.npz"
npz = np.load(libfile, allow_pickle=True)
#libdf = pd.DataFrame(libdf)
libdf= pd.DataFrame.from_dict({item: npz[item] for item in npz.files}, orient='columns')


formulas = libdf["Formula"].to_numpy()
print(formulas)

elist = ["C", "N", "O", "H", "P"]

countlist = []
for form in formulas:
    f = Formula(form)
    comp = f.composition()
    counts = [get_count(comp, c) for c in elist]
    countlist.append(counts)
countlist = np.array(countlist)
print(elist)
print(np.amin(countlist, axis=0))
print(np.amax(countlist, axis=0))


