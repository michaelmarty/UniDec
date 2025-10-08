import os
import pandas as pd
from matchms.importing import load_from_msp
from matchms.exporting import save_as_msp
import numpy as np

def get_specstring(spectrum):
    spec = spectrum.peaks
    spec = np.array(spec)
    specstr = ""
    try:
        for mz, inten in zip(spec[:, 0], spec[:, 1]):
            specstr += f"{mz:.4f}:{inten:.1f} "
    except:
        print("Error in spectrum:", spec)
    specstr = specstr.strip()
    return specstr

def load_df(filename):
    df = pd.read_csv(posfile)
    names = df["Metabolite name"].values
    redone_names = []
    for n in names:
        if "|" in n:
            n = n.split("|")[1]
        redone_names.append(n)
    df["Metabolite name Corr"] = redone_names
    print(df.columns)
    return df

def add_ref_spectra(df, libfile):
    lib = load_from_msp(libfile)
    outdata = []
    for l in lib:

        name = l.metadata["compound_name"]
        print(name)
        if name in df["Metabolite name Corr"].values:
            adduct = l.metadata["adduct"]
            rt = l.metadata["retention_time"]

            # Find row in dfpos that matches name and adduct
            rows = df[df["Metabolite name Corr"] == name]
            rows = rows[rows["Adduct type"] == adduct]
            if len(rows) == 0:
                continue

            print(name, adduct, rt, l)
            outdata.append(l)
            if len(rows) > 1:
                print("Multiple matches for", name, adduct, rows)
            specstr = get_specstring(l)
            for index, row in rows.iterrows():
                df.at[index, "Ref Spec"] = specstr

    print(len(outdata), len(df))
    return df, outdata

def parse_xml(xmlfile):
    pass

topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
os.chdir(topdir)

posfile = "PosIDsC1.csv"
negfile = "NegIDsC1.csv"
dffile = negfile

poslibfile = "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp"
neglibfile = "MSDIAL-TandemMassSpectralAtlas-VS69-Neg-C1.msp"
libfile = neglibfile
outlibfile = libfile.replace(".msp", "_Ref.msp")
outcsvfile = posfile.replace(".csv", "_withRef.csv")

df = load_df(dffile)

exit()
df, outdata = add_ref_spectra(df, libfile)

# Save outdata to new msp file
save_as_msp(outdata, outlibfile, mode="w")
# Save dfpos to new csv file
df.to_csv(outcsvfile, index=False)




