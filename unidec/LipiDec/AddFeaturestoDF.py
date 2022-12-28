import pandas as pd
import numpy as np
import re


def get_overall_tails(name):
    if "|" in name:
        try:
            tname = name.split("|")[0] # Remove before vertical line
            tname = tname.split()[1] # Remove after space to get tails
            tl = tname.split(":")[0] # Get length
            tu = tname.split(":")[1] # Get unsaturation
        except:
            return name, -1, -1
        try:
            tl = int(tl)
        except:
            tl = -1
        try:
            tu = int(tu)
        except:
            tu = -1
    else:
        tname = name.split()[1]
        fas = tname.split("_")
        tls=0
        tus=0
        for f in fas:
            tl = f.split(":")[0]  # Get length
            try:
                tu = f.split(":")[1]  # Get unsaturation
            except:
                print(f, name)
                exit()
            try:
                tl = int(tl)
            except:
                tls = -1
                break
            try:
                tu = int(tu)
            except:
                tus = -1
                break
            tls += tl
            tus += tu
        tl = tls
        tu = tus

    return tname, tl, tu


def set_tails(df, name="Metabolite name"):
    df = df.copy(deep=True)
    names = df[name].to_numpy()
    tus = []
    tls = []
    tnames = []
    for n in names:
        total, tl, tu = get_overall_tails(n)
        tus.append(tu)
        tls.append(tl)
        tnames.append(total)
    df["Tail Lengths"] = tls
    df["Tail Unsaturation"] = tus
    df["Tail Names"] = tnames
    return df


adduct_translator = np.array([["[M+H]+", "[M+H]", "+1", 1, "[M+H]1+"],
                              ["[M+NH4]+", "[M+NH4]", "+1", 1, "[M+NH4]1+"],
                              ["[M+2H]2+", "[M+2H]", "+2", 2, "[M+2H]2+"],
                              ["[M-H]-", "[M-H]", "-1", -1, "[M-H]1-"],
                              ["[M-HCOO]-", "[M-HCOO]", "-1", -1, "[M-HCOO]1-"],
                              ["[M-CH3OO]-", "[M-CH3OO]", "-1", -1, "[M-CH3OO]1-"],
                              ["[M-2H]2-", "[M-2H]", "-2", -2, "[M-2H]2-"]])


def parse_names(df, namecolumn="Metabolite name", parse_adduct=True):
    simp_names = []
    full_names = []
    full_names_adduct = []
    charges = []
    for i, row in df.iterrows():
        #Parse Name
        mol_name = row[namecolumn]
        if "|" in mol_name:
            names = mol_name.split("|")
            simp_name = names[0]
            full_name = names[1]
        else:
            simp_name, full_name = mol_name, mol_name
        simp_names.append(simp_name)
        full_names.append(full_name)

        #Parse Adduct
        if parse_adduct:
            prec_adduct = row["Adduct type"]
            if prec_adduct in adduct_translator[:, 0]:
                index = np.where(adduct_translator[:, 0] == prec_adduct)[0][0]
                adduct_row = adduct_translator[index]
                aname = adduct_row[-1]
                charge = adduct_row[-2]
            else:
                aname = ""
                charge = 0
            charges.append(charge)
            full_names_adduct.append(str(full_name) + str(aname))

    df["Simp Name"] = simp_names
    df["Full Name"] = full_names
    if parse_adduct:
        df["Full Name Adduct"] = full_names_adduct
        df["Charge"] = charges
    return df


kendrick_ref_mass = 14.01565
mass_h = 1.007276467


def kendrick_analysis(mz, charge):
    mass = np.abs(mz * charge) - charge * mass_h
    kmass = mass * 14 / kendrick_ref_mass
    kmd = kmass - np.round(kmass)
    return mass, kmass, kmd


def mass_defect_df(df, ref_col="Reference m/z", ex_col="Average Mz", z_col="Charge"):
    refmz = df[ref_col].to_numpy()
    mz = df[ex_col].to_numpy()
    charge = df[z_col].to_numpy().astype(float)

    mass, kmass, kmd = kendrick_analysis(mz, charge)
    refmass, refkmass, refkmd = kendrick_analysis(refmz, charge)
    deltamd = np.abs(kmd - refkmd)

    df["mass"] = mass
    df["kmass"] = kmass
    df["kmassd"] = kmd
    df["refmass"] = refmass
    df["refkmass"] = refkmass
    df["refkmassd"] = refkmd
    df["deltamd"] = deltamd
    return df


def isotope_ratios(df):
    ratios = []
    for i, row in df.iterrows():
        sraw = row["MS1 isotopic spectrum"]
        s = np.array(re.split(" |:", sraw))
        try:
            ints = s[1::2].astype(float)
            ratio = ints[1]/ints[0]
        except:
            ratio = -1
        ratios.append(ratio)
    df["ms1iratios"] = ratios
    return df


def fill_out_df(df):
    df = isotope_ratios(df)
    df = set_tails(df)
    df = parse_names(df)
    df = mass_defect_df(df)
    return df


if __name__ == "__main__":
    file = "D:\\Data\\Lipidomics\\OpenMS\\FilteredFixedRT_Tails.csv"

    df = pd.read_csv(file)

    df = fill_out_df(df)
    df.to_csv("D:\\Data\\Lipidomics\\OpenMS\\FilteredFixedAll.csv")
