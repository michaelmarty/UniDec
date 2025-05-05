import pandas as pd
import numpy as np
import os

adduct_translator = np.array([["[M+H]+", "[M+H]", "+1", 1, "[M+H]1+"],
                              ["[M+NH4]+", "[M+NH4]", "+1", 1, "[M+NH4]1+"],
                              ["[M+2H]2+", "[M+2H]", "+2", 2, "[M+2H]2+"],
                              ["[M-H]-", "[M-H]", "-1", -1, "[M-H]1-"],
                              ["[M-HCOO]-", "[M-HCOO]", "-1", -1, "[M-HCOO]1-"],
                              ["[M-CH3OO]-", "[M-CH3OO]", "-1", -1, "[M-CH3OO]1-"],
                              ["[M-2H]2-", "[M-2H]", "-2", -2, "[M-2H]2-"]])

library_file = r"Z:\Group Share\JGP\Lipidomics_MassIVE_folder\Raw\DDA\calibrated_mzML\rt_calibrated_filtered_library.csv"
library_df = pd.read_csv(library_file)
os.chdir(os.path.dirname(library_file))

library_names = library_df["Lipid name"]
adducts = library_df["Adduct"]

full_name_adducts = []
for i,name in enumerate(library_names):
    if adducts[i] in adduct_translator[:,0]:
        adduct_index = np.where(adduct_translator[:,0] == adducts[i])[0][0]
        adduct_row = adduct_translator[adduct_index]
        aname = adduct_row[-1]
        acharge = adduct_row[-2]
    else:
        aname = ""
        acharge = 0

    if "|" in name:
        names = name.split("|")
        full_name = names[1]
    else:
        full_name = name

    full_name_adducts.append(str(full_name) + str(aname))

output_df = pd.DataFrame()
output_df["lipid names"] = full_name_adducts
output_df.to_csv("lipid_list.csv", index=False)

