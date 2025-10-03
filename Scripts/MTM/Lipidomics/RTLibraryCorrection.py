import numpy as np

fit = [-1.20739753e-03,  7.55428678e-02, -1.72632321e+00,  1.83605245e+01,
 -9.03502485e+01,  1.67554242e+02]
p = np.poly1d(fit)

libfile = r"C:\Users\mm96978\Downloads\MSDIAL-TandemMassSpectralAtlas-VS69-Pos.msp"
outfile = r"C:\Users\mm96978\Downloads\MSDIAL-TandemMassSpectralAtlas-VS69-Pos-Corrected.msp"

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
