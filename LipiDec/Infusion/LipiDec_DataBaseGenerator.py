import time
from copy import deepcopy
import pandas as pd
import numpy as np
import unidec_modules.unidectools as ud
import os
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem.Descriptors import *
from molmass import Formula
import re


def construct_mol(smiles, subs):
    subdict = {}
    for s in subs:
        key = "${" + s[0] + "}"
        sm = s[4]
        subdict[key] = sm
    m = Chem.MolFromSmiles(smiles, replacements=subdict)
    return m


def massfromstruct(struct, subs, basemods=[]):
    if struct.attrib["type"] == "SMILES":
        sm = struct.text
        sm = re.sub("D", '[2H]', sm)

        if len(basemods) > 0:
            for bm in basemods:
                sm = re.sub(bm[0], bm[1], sm)
        return massfromsmiles(sm, subs)
    else:
        print("ERROR:", struct.tag, struct.attrib, struct.text)
        return -1


def find_basemods(struct):
    basemods = l.findall("baseModification")
    basemodlist = []
    if len(basemods) > 0:
        for bm in basemods:
            mods = bm.findall("modification")
            for submod in mods:
                target = submod.find("targetStructure").text
                # target = re.sub(r"\\", '', target)
                target = re.sub("\$", '', target)
                replace = submod.find("replaceStructure").text
                replace = re.sub(r"D", '[2H]', replace)
                basemodlist.append([target, replace])
    return basemodlist


def massfromsmiles(smiles, subs):
    mol = construct_mol(smiles, subs)
    mass = ExactMolWt(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol, abbreviateHIsotopes=True, separateIsotopes=True)
    return mass, formula


def get_subs(att, subs, cond, mods):
    outputlist = []
    site = att["siteKey"]
    base = att["baseKey"]
    direct = att["direction"]
    conds = cond.attrib
    # print(conds)

    for subtype in subs:
        if subtype.attrib["sgKey"] == base:
            for s in subtype:
                key = s.attrib["subKey"]
                key = re.sub('P-', '', key)
                key = re.sub('O-', '', key)
                # print("Key:", key)

                if "d" in key or "t" in key:
                    key = re.sub('d', '', key)
                    key = re.sub('t', '', key)
                    numC = int(key.split(":")[0]) - 3  # No sure about this one
                    numDB = int(key.split(":")[1])
                else:
                    numC = int(key.split(":")[0])
                    numDB = int(key.split(":")[1])
                if direct == "l":
                    struct = s.find("leftStructure").text
                elif direct == "r":
                    struct = s.find("rightStructure").text
                else:
                    print("ERROR in direction", direct)
                    struct = None

                for mod in mods:
                    submods = mod.findall("modification")
                    for submod in submods:
                        target = submod.find("targetStructure").text
                        target = re.sub(r"\\", '', target)
                        target = re.sub(r"\$", '', target)
                        replace = submod.find("replaceStructure").text
                        replace = re.sub(r"D", '[2H]', replace)
                        if struct == target:
                            struct = replace
                        elif struct[-4:] == target:
                            struct = struct[:-4] + replace
                outputlist.append([site, key, numC, numDB, struct])
    outputlist = np.array(outputlist)
    return outputlist


def make_subs_list(lclass, subs):
    """
    Get the list of substituents for a given class
    :param lclass: Lipid Class XML object
    :param subs: Substituents XML object
    :return: DataFrame with substituent info, numpy array of similar info
    """
    subclass = lclass.attrib["subclassKey"]
    exempt = ["GM3(D5)", "GM1(D5)"]

    subcond = lclass.find("substituentCondition")
    # dbs = subcond.find("doubleBondCondition")[0].attrib
    # print(dbs)  # Note, I don't understand this
    try:
        crange = [int(subcond.attrib["minSumCNum"]), int(subcond.attrib["maxSumCNum"])]
    except:
        crange = [0, 1000]
    try:
        dbrange = [int(subcond.attrib["minSumUNum"]), int(subcond.attrib["maxSumUNum"])]
    except:
        dbrange = [0, 100]

    targets = subcond.findall("substituentTarget")
    sublists = []
    lengths = []
    for t in targets:
        mods = t.findall("subModification")
        sublist = get_subs(t.attrib, subs, subcond, mods)
        sublists.append(sublist)
        lengths.append(len(sublist))

    subindexes = list(np.ndindex(*lengths))

    combinedlist = []
    totalclist = []
    totalulist = []
    for s in subindexes:
        allsubs = []
        totalC = 0
        totaldb = 0
        for i, index in enumerate(s):
            line = sublists[i][index]
            allsubs.append(line)
            totalC += int(line[2])
            totaldb += int(line[3])

        if (crange[0] <= totalC <= crange[1] and dbrange[0] <= totaldb <= dbrange[1]) or len(
                subindexes) == 1 or subclass in exempt:
            combinedlist.append(allsubs)
            totalclist.append(totalC)
            totalulist.append(totaldb)

    combinedlist = np.array(combinedlist)

    cnames = ["Site", "Name", "C", "U", "SMILES"]
    # print(combinedlist)
    cdf = pd.DataFrame()
    for i in range(len(lengths)):
        for j, c in enumerate(cnames):
            cdf["Site" + str(i) + "_" + c] = combinedlist[:, i, j]

    cdf["TotalC"] = totalclist
    cdf["TotalU"] = totalulist

    sc1 = np.core.defchararray.add(cdf["TotalC"].to_numpy(dtype=str), ":")
    cdf["SumComp"] = np.core.defchararray.add(sc1, cdf["TotalU"].to_numpy(dtype=str))

    print("Number of Substituent Combinations", len(cdf))
    return cdf, combinedlist


def adduct_mass(name):
    """
    Calculate the mass shift of an adduct from the name
    :param name: The name from the adduct database
    :return: Monoisotopic Mass
    """
    atoms = name.split()
    totalmass = 0
    for atom in atoms:
        if "-" in atom:
            sign = -1
            formula = atom.replace("-", "")
        else:
            sign = 1
            formula = atom
        mass = Formula(formula).monoisotopic_mass * sign
        totalmass += mass
    return totalmass


def make_adduct_db(adducts):
    # Create a DataFrame from the attributes
    attributelist = [a.attrib for a in adducts]
    adf = pd.DataFrame(attributelist)
    # Calculate and add the mass shift for each adduct
    adf["Mass"] = np.array([adduct_mass(name) for name in adf["name"].to_numpy()])
    return adf


def sum_comp_only(df):
    columns = df.keys()
    droplist = ["site", "fa-chain"]
    for c in columns:
        if "Site" in c:
            droplist.append(c)
    df = df.drop(columns=droplist)
    # df = df.reset_index()

    df = df.drop_duplicates(subset=["charge", "SumComp", "Formula", "Adduct", "subclassKey"], ignore_index=True)
    df = df.sort_values("Mz")

    return df


libfile = "C:\\Data\\Lipidomics\\Libraries\\191219-190502-LSDB5035-YHR-Human-NH4HCOO-NonDer_200929.xml"
outfile = "C:\\Data\\Lipidomics\\Libraries\\nonder1.xlsx"

libfile = "C:\\Data\\Lipidomics\\Libraries\\191219-190502-LSDB5035-YHR-Human-NH4HCOO-Der_200929.xml"
outfile = "C:\\Data\\Lipidomics\\Libraries\\der1.xlsx"

tree = ET.parse(libfile)
root = tree.getroot()

st = time.perf_counter()

# Unpack top database
lclasses = []
for r in root:
    if r.tag == "adductIonDefinition":
        adducts = r
    elif r.tag == "substituentDefinition":
        subs = r
    elif r.tag == "lipidIonCondition":
        lclasses.append(r)
    else:
        print("Unable to parse:", r.tag)

adf = make_adduct_db(adducts)

classlist = ["PC", "PG", "CL", "PE", "PS", 'PC(D7)', "PC-P"]
classlist = ["PC-P", "PC", "PE", "PE-P", "PC(D7)"]
classlist = None
#classlist = ["PE-P(D9)"]

ignore_list = ["TG", "TG-O", "TG-P", "TG(D5)", "AcCer"]
ignore_list = []

classes = [l.attrib["subclassKey"] for l in lclasses]
classchar = [len(l.attrib["subclassKey"]) for l in lclasses]
sortindexes = np.argsort(classchar)
print(classes)

# Loop through lipid classes
topdf = pd.DataFrame()
for i in sortindexes:
    l = lclasses[i]
    att = l.attrib
    lclass = att["subclassKey"]
    if (classlist is None or lclass in classlist) and lclass not in ignore_list:
        print(lclass)
        # Base structure of head group
        struct = l.find("baseStructure")

        basemodlist = find_basemods(struct)

        # Substituents
        moldf, sublist = make_subs_list(l, subs)
        # Fill out DF
        for d in att:
            moldf[d] = att[d]

        # Get Masses for all lipids
        masses = []
        formulas = []
        for s in sublist:
            try:
                m, f = massfromstruct(struct, s, basemods=basemodlist)
                masses.append(m)
                formulas.append(f)
            except Exception as e:
                print("ERROR:", e)
                print(struct, s)

                masses.append(0)
                formulas.append("")
        moldf["Mass"] = masses
        moldf["Formula"] = formulas

        # Find the adduct ion conditions
        adducts = l.find("adductIonCondition").text
        adductlist = adducts.split(";")
        # Loop through the adducts, add their mass and information into the growing df
        for adduct in adductlist:
            row = adf[adf["ionKey"] == adduct]
            charge = float(row["charge"])
            adductmass = float(row["Mass"])

            newdf = deepcopy(moldf)
            newdf["charge"] = charge
            newdf["Mz"] = (newdf["Mass"] + adductmass) / np.abs(charge)
            newdf["Adduct"] = adduct
            newdf["AdductMass"] = adductmass

            topdf = pd.concat([topdf, newdf])

topdf = sum_comp_only(topdf)

outdict = {item: topdf[item] for item in topdf.keys()}
np.savez(outfile[:-5] + ".npz", **outdict)

# Write output file
topdf.to_excel(outfile)
print("Database Size:", len(topdf))
# print(topdf.keys())
print("Done:", time.perf_counter() - st)
