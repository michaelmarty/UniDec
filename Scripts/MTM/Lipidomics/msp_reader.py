import os
from matchms import importing as im
import pickle
import numpy as np
import pandas as pd
import time

negfile = "Z:\Group Share\MTM\MSDIAL-TandemMassSpectralAtlas-VS69-Neg.msp"
posfile = "Z:\Group Share\MTM\MSDIAL-TandemMassSpectralAtlas-VS69-Pos.msp"
negpkl = "neglib.pkl"
pospkl = "poslib.pkl"
os.chdir(os.path.dirname(negfile))

if False:
    if os.path.isfile(negpkl):
        ndat = pickle.load(open(negpkl, "rb"))
        print("Loaded Neg Pkl:", len(ndat))
    else:
        ndat = im.load_from_msp(negfile)
        pickle.dump(list(ndat), open(negpkl, "wb"))
        print("Saved Neg pkl:", len(list(ndat)))
    df = None
    for i, n in enumerate(ndat):
        dict = n.metadata
        if df is None:
            df = pd.DataFrame(columns=dict.keys())
            columns = dict.keys()

        if dict.keys() == columns:
            vals = dict.values()
        else:
            vals = []
            for c in columns:
                try:
                    val = dict[c]
                except:
                    val = ""
                vals.append(val)
        if i % 1000 == 0:
            print(i)

        df.loc[len(df.index)] = vals

    df.to_csv("neglib.csv")
else:
    if os.path.isfile(pospkl):
        ndat = pickle.load(open(pospkl, "rb"))
        print("Loaded Pos Pkl:", len(ndat))
    else:
        ndat = im.load_from_msp(posfile)
        pickle.dump(list(ndat), open(pospkl, "wb"))
        print("Saved Pos pkl:", len(list(ndat)))

    dlen = len(list(ndat))
    tstart = time.perf_counter()
    df = None
    lines = []
    for i, n in enumerate(ndat):
        dict = n.metadata
        if i == 0:
            columns = dict.keys()
        if dict.keys() == columns:
            vals = dict.values()
        else:
            vals = []
            for c in columns:
                try:
                    val = dict[c]
                except:
                    val = ""
                vals.append(val)
        if i % 1000 == 0 and i != 0:
            newt = time.perf_counter()
            time_passed = newt - tstart
            time_total = time_passed * dlen / i
            time_remaining = time_total - time_passed
            print(i, time_passed, time_total, time_remaining)

        lines.append(vals)

    df = pd.DataFrame(lines, columns=dict.keys())

    df.to_csv("poslib.csv")
