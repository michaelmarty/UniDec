import numpy as np
import pickle

# fit = [ 2.42393107e-05, -1.59066957e-03,  3.61629338e-02, -3.54583941e-01,
#   2.46833668e+00, -1.96586339e+00]
# p = np.poly1d(fit)

interpolationfunction = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\RTInterpolation.pkl"
with open(interpolationfunction, 'rb') as f:
    p = pickle.load(f)

libfile = r"C:\Users\mm96978\Downloads\MSDIAL-TandemMassSpectralAtlas-VS69-Pos.msp"
outfile = r"C:\Users\mm96978\Downloads\MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp"

# Read library file
with open(libfile, 'r') as f:
    lines = f.readlines()

# Write new library file correcting RTs
with open(outfile, 'w') as f:
    for i, line in enumerate(lines):
        if line.startswith("RETENTIONTIME:"):
            rt = float(line.split(":")[1].strip())
            newrt = p(rt)
            f.write(f"RETENTIONTIME: {newrt:.4f}\n")
            print(f"Corrected RT: {rt} -> {newrt:.4f}")
        else:
            f.write(line)
