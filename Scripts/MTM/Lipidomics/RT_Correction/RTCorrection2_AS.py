import numpy as np
import pickle
import os

os.chdir(r"Z:\Group Share\Annika\Stellar\HEK")

poslibfile = r"MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp"
neglibfile = r"MSDIAL-TandemMassSpectralAtlas-VS69-Neg-C1.msp"

interpolationfunction = r"RTInterpolation.pkl"
with open(interpolationfunction, 'rb') as f:
    p = pickle.load(f)

# Rename outfile to replace C1 with Stellar 1
posoutfile = poslibfile.replace("-C1", "-Stellar1")
negoutfile = neglibfile.replace("-C1", "-Stellar1")

# Read library file
with open(poslibfile, 'r') as f:
    lines = f.readlines()

# Write new library file correcting RTs
with open(posoutfile, 'w') as f:
    for i, line in enumerate(lines):
        if line.startswith("RETENTIONTIME:"):
            rt = float(line.split(":")[1].strip())
            newrt = p(rt)
            f.write(f"RETENTIONTIME: {newrt:.4f}\n")
            print(f"Corrected RT: {rt} -> {newrt:.4f}")
        else:
            f.write(line)

with open(neglibfile, 'r') as f:
    lines = f.readlines()

# Write new library file correcting RTs
with open(negoutfile, 'w') as f:
    for i, line in enumerate(lines):
        if line.startswith("RETENTIONTIME:"):
            rt = float(line.split(":")[1].strip())
            newrt = p(rt)
            f.write(f"RETENTIONTIME: {newrt:.4f}\n")
            print(f"Corrected RT: {rt} -> {newrt:.4f}")
        else:
            f.write(line)
