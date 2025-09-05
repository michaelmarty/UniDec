import os
import pandas as pd
import pickle
import re

henc = {'ADGGA': 0, 'AHexBRS': 1, 'AHexCAS': 2, 'AHexCS': 3, 'AHexCer': 4, 'AHexSIS': 5, 'AHexSTS': 6, 'ASM': 7, 'Ac2PIM1': 8, 'Ac2PIM2': 9, 'Ac3PIM2': 10, 'Ac4PIM2': 11, 'BASulfate': 12, 'BMP': 13, 'BRSE': 14, 'BileAcid': 15, 'CAR': 16, 'CASE': 17, 'CE': 18, 'CL': 19, 'CerP': 20, 'Cer_ADS': 21, 'Cer_AP': 22, 'Cer_AS': 23, 'Cer_BDS': 24, 'Cer_BS': 25, 'Cer_EBDS': 26, 'Cer_EODS': 27, 'Cer_EOS': 28, 'Cer_NDS': 29, 'Cer_NP': 30, 'Cer_NS': 31, 'CoQ': 32, 'DCAE': 33, 'DG': 34, 'DGCC': 35, 'DGDG': 36, 'DGGA': 37, 'DGTS': 38, 'DHSph': 39, 'DLCL': 40, 'EtherDG': 41, 'EtherDGDG': 42, 'EtherLPC': 43, 'EtherLPE': 44, 'EtherLPG': 45, 'EtherMGDG': 46, 'EtherOxPC': 47, 'EtherOxPE': 48, 'EtherPC': 49, 'EtherPE': 50, 'EtherPG': 51, 'EtherPI': 52, 'EtherPS': 53, 'EtherSMGDG': 54, 'EtherTG': 55, 'FA': 56, 'FAHFA': 57, 'GDCAE': 58, 'GLCAE': 59, 'GM3': 60, 'HBMP': 61, 'Hex2Cer': 62, 'Hex3Cer': 63, 'HexCer_AP': 64, 'HexCer_EOS': 65, 'HexCer_HDS': 66, 'HexCer_HS': 67, 'HexCer_NDS': 68, 'HexCer_NS': 69, 'LDGCC': 70, 'LDGTS': 71, 'LNAPE': 72, 'LNAPS': 73, 'LPA': 74, 'LPC': 75, 'LPE': 76, 'LPG': 77, 'LPI': 78, 'LPS': 79, 'MG': 80, 'MGDG': 81, 'MLCL': 82, 'NAE': 83, 'NAGly': 84, 'NAGlySer': 85, 'NAOrn': 86, 'OxPC': 87, 'OxPE': 88, 'OxPG': 89, 'OxPI': 90, 'OxPS': 91, 'PA': 92, 'PC': 93, 'PE': 94, 'PE_Cer': 95, 'PEtOH': 96, 'PG': 97, 'PI': 98, 'PI_Cer': 99, 'PMeOH': 100, 'PS': 101, 'PhytoSph': 102, 'SHex': 103, 'SHexCer': 104, 'SISE': 105, 'SL': 106, 'SM': 107, 'SQDG': 108, 'SSulfate': 109, 'STSE': 110, 'Sph': 111, 'TDCAE': 112, 'TG': 113, 'TLCAE': 114, 'VAE': 115, 'Vitamin_D': 116, 'Vitamin_E': 117}

def lipid_name_to_vector(name):
    # Parse a lipid name into a one-hot vector representation
    parts = name.split(" ")
    headgroup = parts[0]
    tails = re.split(r"[\/_]", parts[1]) if len(parts) > 1 else []
    tail_info = []
    for tail in tails:
        if ":" in tail:
            nums = tail.split(":")
            try:
                tail_info.append((int(nums[0]), int(nums[1])))
            except Exception as e:
                print("Error parsing tail:", tail, e)
                continue
        else:
            tail_info.append((0, 0))
        enc = henc[headgroup]
        print(name, headgroup, enc, tail_info)



topdir = r"Z:\Group Share\Annika\Stellar"
file1 = "Thermo Merged LipidList.csv"

df = pd.read_csv(os.path.join(topdir, file1))

molecules = df["Molecule Name"].tolist()
print("Total molecules:", len(molecules))

for m in molecules:
    lipid_name_to_vector(m)

# headgroups = [m.split(" ")[0] for m in molecules]
# uniqueheadgroups = set(headgroups)
# print("Unique headgroups:", len(uniqueheadgroups))
# print(uniqueheadgroups)

#
# libdir = r"Y:\Data\Lipidomics\Libraries"
# file2 = "neglib.pkl"
# ndat = pickle.load(open(os.path.join(libdir, file2), "rb"))
# print("Loaded Pkl:", len(ndat))
# for i, n in enumerate(ndat):
#     dict = n.metadata
#
#     if i == 0:
#         columns = dict.keys()
#     if dict.keys() == columns:
#         vals = dict.values()
#     print(columns)
#     print(n)
#     exit()