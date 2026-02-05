import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def find_match(libfile, compoundname):
    with open(libfile) as f:
        while True:
            line = f.readline()
            if not line:
                break
            #line = next((l for l in f if compoundname in l), None)
            if compoundname in line:
                print(line.strip())
    print("Done")

'''
libcsvpath = "C:\\Data\\Lipidomics\\neglib.csv"
os.chdir(os.path.dirname(libcsvpath))

df = pd.read_csv(libcsvpath)
print(df.keys())

df["compoundclass"].astype("category")
df["compoundcodes"] = df["compoundclass"].astype("category").cat.codes

mz = df["precursor_mz"].to_numpy()
rt = df["retention_time"].to_numpy()
ccs = df["ccs"].to_numpy()'''

libfile = "C:\\Data\\Lipidomics\\neglib.csv"
libfile = "C:\\Data\\Lipidomics\\MSDIAL-TandemMassSpectralAtlas-VS69-Neg.msp"
libfile = "C:\\Users\\marty\\Downloads\MSDIAL ver.4.90 Windows\MSDIAL ver.4.90 Windows\Msp20220322084048_converted.lbm2"

toppath = "C:\\Data\\Lipidomics\\OpenMS\\"
os.chdir(toppath)

df = pd.read_csv("Filtered.csv")
print(df.keys())

n=25
name1 = df.iloc[n]["Full Name"]
print(name1)
print(df.iloc[n])

find_match(libfile, name1)
