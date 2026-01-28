import numpy as np
import pickle
import os

os.chdir(r"Z:\Group Share\Annika\Stellar\Untargeted DDA\HEK")

poslibfile = r"MSDIAL-TandemMassSpectralAtlas-VS69-Pos-StellarRTP2.msp"
neglibfile = r"MSDIAL-TandemMassSpectralAtlas-VS69-Neg-StellarRTP2.msp"

interpolationfunction = r"RTInterpolation.pkl"
with open(interpolationfunction, 'rb') as f:
    p = pickle.load(f)

posoutfile = poslibfile.replace("-RTP2", "-RTP3")
negoutfile = neglibfile.replace("-RTP2", "-RTP3")

# -------- POS --------
with open(poslibfile, 'r') as f:
    poslines = f.readlines()

with open(posoutfile, 'w') as f:
    for line in poslines:
        if line.startswith("RETENTIONTIME:"):
            rt = float(line.split(":")[1].strip())
            newrt = p(rt)
            f.write(f"RETENTIONTIME: {newrt:.4f}\n")
        else:
            f.write(line)

# -------- NEG --------
with open(neglibfile, 'r') as f:
    neglines = f.readlines()

with open(negoutfile, 'w') as f:
    for line in neglines:
        if line.startswith("RETENTIONTIME:"):
            rt = float(line.split(":")[1].strip())
            newrt = p(rt)
            f.write(f"RETENTIONTIME: {newrt:.4f}\n")
        else:
            f.write(line)
