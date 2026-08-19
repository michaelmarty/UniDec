from UniDecImporter import ImporterFactory
import pandas as pd
import os
import subprocess
import numpy as np

def get_lipidcreator_name(name, polarity):
    if "light" in name:
        split = name.split(" (light)")
    elif "heavy" in name:
        split = name.split(" (heavy)")

    if polarity.ToString() == "Negative":
        adduct = "[M-H]1-"
    if polarity.ToString() == "Positive":
        adduct = "[M+H]1+"

    fullname = str(split[0] + adduct)
    name = split[0]
    return fullname, name

directory = r"Z:\Group Share\Annika\2025-04-29_OE240_Stellar_UTAustin\Analytical Flow\Stellar\Data"
os.chdir(directory)
importer = ImporterFactory.create_importer("2025-03-13_Stellar_UTAustin_Sample_A_posneg_tMS2_01.raw")

lipid_prec_dict = {}

for i in importer.scans:
    scanfilter = importer.msrun.source.GetFilterForScanNumber(int(i))
    polarity = scanfilter.get_Polarity()
    scanEvent = importer.msrun.source.GetScanEventForScanNumber(int(i))
    reaction = scanEvent.GetReaction(0)
    prec_mass = reaction.PrecursorMass
    name = scanfilter.Name
    if name not in lipid_prec_dict:
        lipid_prec_dict[name] = {}
        lipid_prec_dict[name][polarity] = prec_mass
    else:
        if polarity in lipid_prec_dict[name]:
            continue
        else:
            lipid_prec_dict[name][polarity] = prec_mass

polarities = []
masses = []
lipidcreator_names = []
fullnames = []

for lipid in lipid_prec_dict:
    for polarity in lipid_prec_dict[lipid]:
        polarities.append(polarity)
        masses.append(lipid_prec_dict[lipid][polarity])
        lipidcreator_name, fullname = get_lipidcreator_name(lipid,polarity)
        lipidcreator_names.append(lipidcreator_name)
        fullnames.append(fullname)

unique_lipids = np.array(lipidcreator_names)
np.savetxt("unique_lipids.csv", unique_lipids, fmt="%s", delimiter=",")

lipidcreator_path = r"C:\Users\Admin\Downloads\LipidCreator\LipidCreator.exe"

infile = "unique_lipids.csv"
outfile = "transition_list.csv"

command = [lipidcreator_path, "transitionlist", infile, outfile]
print("Running command:", " ".join(command))
subprocess.run(command)

transitionlist_df = pd.read_csv(outfile)
prec_names = transitionlist_df["PrecursorName"]
product_names = transitionlist_df["ProductName"]
product_mzs = transitionlist_df["ProductMz"]

prec_name_transitions = {}
for i in range(len(prec_names)):
    if prec_names[i]in prec_name_transitions:
        prec_name_transitions[prec_names[i]].append(product_mzs[i])
    else:
        prec_name_transitions[prec_names[i]] = [product_mzs[i]]


print("Finished creating transition list!")