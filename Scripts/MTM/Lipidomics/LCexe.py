import sys, os, subprocess
import pandas as pd

exe = r"C:\Users\mm96978\AppData\Local\Apps\2.0\ND2QYAB6.QT6\VET5K44V.LYE\skyl..tion_9286511f3362df93_0019.0001_31a2f7d2aa6cc5be\Tools\SkyE3AF\LipidCreator.exe"

def run_LC(infile, outfile):
    # helpcall = [exe, "transitionlist", "help"]
    tlcall = [exe, "transitionlist", infile, outfile, "-h", "0"]
    proc = subprocess.Popen(tlcall, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    # print(err)
    # print(out.decode())
    return out.decode()

def ParseErrors(errstr):
    errlines = errstr.split("\n")
    badlipids = []
    for line in errlines:
        if "could not parsed" in line:
            # print("LipidCreator Error:", line)
            # Select everything betwee ''
            lipid = line.split("'")[1]
            print("  Problematic Lipid:", lipid)
            badlipids.append(lipid)
    return badlipids


topdir = r"Y:\Data\Lipidomics\Libraries"
topdir = r"Z:\Group Share\Annika\Stellar"
topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
os.chdir(topdir)

infile1 = "Combined_Lipids.csv"

df = pd.read_csv(infile1)
# Drop some bad classes
badclasses = ["ASM", "GM3", "HexCer_NDS"]
for bc in badclasses:
    df = df[df["Molecule List Name"] != bc]


classes = df["Molecule List Name"].unique()
print(classes)

infile = "Combined_LipidList.csv"
outfile = "LipidCreatorTransitionList.csv"

inpath = os.path.join(topdir, infile)
outpath = os.path.join(topdir, outfile)
inpath2 = os.path.join(topdir, infile.replace(".csv", "_Edited.csv"))

# Drop class "ASM" from df and write others to inpath2

df[["Full Name Adduct"]].to_csv(inpath2, index=False, header=False)
errors = run_LC(inpath2, outpath)

badlipids = ParseErrors(errors)

df = pd.read_csv(infile1)
# Combine badclasses and badlipids from topdf into subdf
subdf = df[df["Molecule List Name"].isin(badclasses)]
# Subdf2 is rows where Full Name Adduct is in badlipids
subdf2 = df[df["Full Name Adduct"].isin(badlipids)]
subdf = pd.concat([subdf, subdf2])

print(subdf)
# Write subdf to new file of bad lipids
badfile = infile1.replace(".csv", "_BadLipids.csv")
subdf.to_csv(badfile, index=False)

# Read Transition List file and say how many unique Molecule Names are in it
tldf = pd.read_csv(outpath)
print(f"Successfully created {outpath} with {tldf['PrecursorName'].nunique()} unique Molecule Names")

# for cls in classes:
#     print(f"Processing class {cls}")
#     df2 = df[df["Molecule List Name"] == cls]
#     # Write Full Name Adduct to inpath2
#     df2[["Full Name Adduct"]].to_csv(inpath2, index=False, header=False)
#     run_LC(inpath2, outpath.replace(".csv", f"_{cls}_Transitions.csv"))



#
# #Loop through the lines and write each to inpath2
# with open(inpath, 'r') as f:
#     lines = f.readlines()
#
# # Remove \r\n from the end of each line and write to CSV
# with open(inpath2, 'w', newline="") as f:
#     for i, line in enumerate(lines):
#         if "PE" not in line:
#             continue
#         line = line.strip()
#         print(line)
#         f.write(line)
#         f.write("\n")



# tcall = [exe, "translate", inpath, inpath2]
# proc = subprocess.Popen(tcall, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# out, err = proc.communicate()
# print(err)#, out.decode())
# print(out.decode())
# print(f"Translated {inpath} to {inpath2}")


