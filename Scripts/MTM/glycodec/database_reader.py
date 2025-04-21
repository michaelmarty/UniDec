import numpy as np
import os
import pandas as pd
import glycodec_eng as gd
import time


def remove_redundant(df):
    names = df["IUPAC Condensed"]
    numbers = df["Accession Number"].to_numpy()
    b1 = np.ones(len(names)).astype(bool)
    '''
    alts = []
    for i, s1 in enumerate(names):
        for j, s2 in enumerate(names):
            if i < j and b1[j] and b1[i]:
                if s1 == s2:
                    b1[j] = False
    df["Unique"] = b1
    '''
    df.reset_index()

    structs = []
    masses = []
    for i, row in df.iterrows():
        eng = gd.GlycoDecEngine()
        name = row["IUPAC Condensed"]
        try:
            s1 = eng.load_iupac(name, attrs=row)
            m = s1.mass()
        except:
            s1 = False
            m = 0
        structs.append(s1)
        masses.append(m)

    alts = []
    for i, s1 in enumerate(structs):
        alt = ""
        for j, s2 in enumerate(structs):
            if i < j and b1[j] and b1[i] and masses[i] == masses[j]:
                same = s1 == s2
                if same or not s2:
                    b1[j] = False
                    alt = alt + " " + numbers[j]
        alts.append(alt)
    df["Alt Accession"] = alts
    return df[b1.astype(bool)]

def find_n_core(df):
    b1 = np.zeros(len(df)).astype(bool)
    for i, row in df.iterrows():
        eng = gd.GlycoDecEngine()
        name = row["IUPAC Condensed"]
        try:
            eng.load_iupac(name, attrs=row)
            r = eng.test_n_core()
            b1[i] = r
        except:
            pass
    df["NCore"] = b1
    return df

starttime = time.perf_counter()
enzyme_sequence = ["ABS", "BTG", "GUH"]

dbpath = "glycosmos_glycans_human_nonredundant.csv"

outpath = "filtered_human.csv"

df = pd.read_csv(dbpath)
#find_n_core(df)
#df.to_csv(dbpath)


df = df[df['IUPAC Condensed'].notna()]
'''
df = df[df["Taxonomy"] == "Homo sapiens:9606"]
print("Filter Taxonomy:", len(df))

df = remove_redundant(df)
print("Removed Redundants:", len(df))
df.to_csv("glycosmos_glyans_human_nonredundant.csv")
print("Time:", time.perf_counter() - starttime)
exit()
'''

df["AllCleaved"] = False
df["AnyCleaved"] = False
df["Mass0"] = 0
for i, e in enumerate(enzyme_sequence):
    cname = "Mass" + str(i + 1)
    df[cname] = 0
    cname2 = "DM" + str(i + 1)
    df[cname2] = 0

for i, row in df.iterrows():
    eng = gd.GlycoDecEngine()
    name = row["IUPAC Condensed"]
    try:
        eng.load_iupac(name, attrs=row)
        allc, anyc, masses = eng.test_enzyme_sequence(enzyme_sequence)
        for j, m in enumerate(masses):
            df.loc[i, "Mass" + str(j)] = m
            if i > 0:
                df.loc[i, "DM" + str(j)] = m - masses[j - 1]

        df.loc[i, "AllCleaved"] = allc
        df.loc[i, "AnyCleaved"] = anyc
        # print("Result", row)
    except Exception as e:
        print("Error with", name)
        print(e)

df2 = df[df["AnyCleaved"]]
df2.to_csv(outpath)
print("Removed uncleaved and bad imports:", len(df2))
print("Time:", time.perf_counter() - starttime)
