import numpy as np
import pickle
import os

os.chdir(r"Z:\Group Share\Annika\IQ-X\Feb 2026 UntDDA\MSDial\Marked Low Score")

poslibfile = r"IQX_pos_OTOT_RTP2.msp"
neglibfile = r"IIQX_neg_OTOT_RTP2.msp"

interpolationfunction = r"RTInterpolation.pkl"
with open(interpolationfunction, 'rb') as f:
    p = pickle.load(f)

# Rename outfile to replace C1 with Stellar 1
posoutfile = poslibfile.replace("-OTOT_RTP2", "-OTIT_RTP2")
negoutfile = neglibfile.replace("-OTOT_RTP2", "-OTIT_RTP2")

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
            # print(f"Corrected RT: {rt} -> {newrt:.4f}")
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
            # print(f"Corrected RT: {rt} -> {newrt:.4f}")
        else:
            f.write(line)
