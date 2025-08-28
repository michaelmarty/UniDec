import numpy as np

adduct_translator = np.array([["[M+H]+", "[M+H]", "+1", 1, "[M+H]1+"],
                              ["[M+NH4]+", "[M+NH4]", "+1", 1, "[M+NH4]1+"],
                              ["[M+2H]2+", "[M+2H]", "+2", 2, "[M+2H]2+"],
                              ["[M-H]-", "[M-H]", "-1", -1, "[M-H]1-"],
                              ["[M-HCOO]-", "[M-HCOO]", "-1", -1, "[M-HCOO]1-"],
                              ["[M-CH3OO]-", "[M-CH3OO]", "-1", -1, "[M-CH3OO]1-"],
                              ["[M-2H]2-", "[M-2H]", "-2", -2, "[M-2H]2-"],
                              ["[M+CH3COO]-", "[M+CH3COO]", "-1", -1, "[M+CH3COO]1-"]])

def lipidlist_to_LC(df):
    outnames = []
    fixednames = []
    for i, row in df.iterrows():
        name = row["Molecule Name"]
        name = fix_dash(name)
        fixednames.append(name)
        adduct = row["Precursor Adduct"]
        print(name, adduct)
        if adduct in adduct_translator[:, 0]:
            index = np.where(adduct_translator[:, 0] == adduct)[0][0]
            adduct_row = adduct_translator[index]
            adduct = adduct_row[1]

        outname = f"{name} {adduct}"
        outnames.append(outname)
    df["Molecule Name"] = fixednames
    df["Full Name Adduct"] = outnames
    return df

def fix_dash(name):
    if "-" in name:
        return name.replace("-", "_")
    return name