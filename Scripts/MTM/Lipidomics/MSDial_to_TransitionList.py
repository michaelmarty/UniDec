import os
import pandas as pd
import numpy as np

def get_charge(adduct):
    if adduct in ["[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+"]:
        return 1
    elif adduct in ["[M+2H]2+"]:
        return 2
    elif adduct in ["[M-H]-", "[M+HCOO]-", "[M+CH3COO]-"]:
        return -1
    elif adduct in ["[M-2H]2-"]:
        return -2
    else:
        print("Unknown adduct", adduct)
        return 0

def add_product_ions(row, cutoff=0.1):
    msmsspec = row["MS/MS spectrum"]
    # Parse spectrum
    msmsspec = [a.split(":") for a in msmsspec.split(" ")]
    msmsspec = np.array([(float(a[0]), float(a[1])) for a in msmsspec])

    # Filter by cutoff
    b1 = msmsspec[:, 1] >= cutoff * np.max(msmsspec[:, 1])
    msmsspec = msmsspec[b1]

    mz = msmsspec[:, 0]
    newrows = []
    for m in mz:
        newrow = row.copy()
        charge = row["Precursor Charge"]
        if charge > 1:
            # Assume charge 1 for product ions
            charge = 1
        elif charge < -1:
            charge = -1
        newrow["Product Mz"] = m
        newrow["Product Charge"] = charge
        newrows.append(newrow)

    return newrows

def msdial_to_transitionlist(df):
    # Drop all where "Ontology" is Others or Unknown
    # df = df[~df["Ontology"].isin(["Others", "Unknown"])]
    newdf = df.copy()
    # Export Df with "Molecule List Name", "Molecule Name", "Precursor Adduct"
    newdf = newdf[["Ontology", "Metabolite name", "Adduct type", "Average Mz", "Formula", "Average Rt(min)", "MS/MS spectrum"]]
    newdf = newdf.rename(columns={"Ontology": "Molecule List Name", "Metabolite name": "Molecule Name",
                            "Adduct type": "Precursor Adduct", "Average Mz": "Precursor Mz",
                                  "Average Rt(min)":"Explicit Retention Time", "Formula":"Molecule Formula"})

    charges = []
    for i, row in newdf.iterrows():
        prec_adduct = row["Precursor Adduct"]
        charge = get_charge(prec_adduct)
        charges.append(charge)
    newdf["Precursor Charge"] = charges

    # Set Explicit Retention Time Window to 4
    newdf["Explicit Retention Time Window"] = 4

    # Sort by Molecule List Name and Molecule Name
    newdf = newdf.sort_values(by=["Molecule List Name", "Molecule Name"])

    newrows = []
    # Get product ion list
    for i, row in newdf.iterrows():
        prodrows = add_product_ions(row)
        newrows.extend(prodrows)

    # New rows to df
    proddf = pd.DataFrame(newrows)

    # Drop MS/MS spectrum column
    proddf = proddf.drop(columns=["MS/MS spectrum"])


    print(proddf.columns)
    return proddf
    pass

topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
os.chdir(topdir)

infile = "NegIDsC1.csv"
outfile = "Neg_TransitionList.csv"

df = pd.read_csv(infile)

print(df.columns)

df = msdial_to_transitionlist(df)
df.to_csv(outfile, index=False)

# Repeate for Positive Mode
infile = "PosIDsC1.csv"
outfile = "Pos_TransitionList.csv"

posdf = pd.read_csv(infile)
posdf = msdial_to_transitionlist(posdf)
posdf.to_csv(outfile, index=False)

# Merge pos and neg
df = pd.concat([df, posdf])
# Write to csv
df.to_csv("Combined_TransitionList.csv", index=False)

# Needed: Precursor Ion Formula	Precursor Neutral Formula
# Needed: Product Mz	Product Charge	Product Ion Formula	Product Neutral Formula	Product Adduct