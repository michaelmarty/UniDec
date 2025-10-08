import os
import pandas as pd
from matchms.importing import load_from_msp
from matchms.exporting import save_as_msp
import numpy as np
import lxml

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
    df = pd.read_csv(filename)
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

def add_ref_from_xml(df, refdf):
    # Match SpotID in refdf to Alignment ID in df
    for index, row in refdf.iterrows():
        spotid = row["SpotID"]
        spec = row["MS/MS spectrum ref"]
        # Find rows in df that match Alignment ID
        rows = df[df["Alignment ID"] == int(spotid)]
        if len(rows) == 0:
            continue
        for i, r in rows.iterrows():
            df.at[i, "Ref Spec"] = spec

    # Print how many rows have Ref Spec empty
    n_empty = len(df[df["Ref Spec"].isnull()])
    print("N rows without Ref Spec:", n_empty, "out of", len(df))
    return df

def parse_xml(xmlfile):
    tree = lxml.etree.parse(xmlfile)
    root = tree.getroot()
    # Extract spots
    spots = root.find("Spots")
    print("N Spots:", len(spots))

    rows = []
    # Loop through spots
    for spot in spots:
        row = {}
        spotid = -1
        adduct = ""
        name = ""
        mz = -1
        spectrum = None
        specstring = ""
        for elem in spot:
            # Extract SpotID
            if elem.tag == "SpotId":
                spotid = elem.text
                # print("SpotID:", spotid)
                continue

            # If "Reference" element exists, strip
            if elem.tag == "Reference":
                adduct = elem.find("Adduct").text
                name = elem.find("Name").text
                mz = elem.find("Mz").text
                # print("Reference:", adduct, name, mz)
                continue

            spectrum = elem.find("MatchedSpectrum")
            if spectrum is not None:
                if len(spectrum) == 0:
                    continue
                for subelem in spectrum:
                    # If Peak is "Reference", strip
                    if subelem.tag == "MatchedSpectrumPeak":
                        for subsubelem in subelem:
                            if subsubelem.tag == "Reference":
                                # Extract Mz and Intensity
                                mz = subsubelem.find("Mz").text
                                inten = subsubelem.find("Intensity").text
                                specstring += f"{mz}:{inten} "



        if spotid != -1 and specstring != "":
            row["SpotID"] = spotid
            row["Adduct"] = adduct
            row["Name"] = name
            row["Mz"] = mz
            row["MS/MS spectrum ref"] = specstring
            rows.append(row)

    columns = ["SpotID", "Adduct", "Name", "Mz", "MS/MS spectrum ref"]
    df = pd.DataFrame(rows, columns=columns)

    return df


topdir = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\FinalResults"
os.chdir(topdir)

# posfile = "PosIDsC1.csv"
# posdf = load_df(posfile)
# refposdf = parse_xml("Pos.xml")
# posdf = add_ref_from_xml(posdf, refposdf)
# posdf.to_csv("PosIDsC1_Ref.csv", index=False)

negfile = "NegIDsC1.csv"
negdf = load_df(negfile)
refnegdf = parse_xml("Neg.xml")
negdf = add_ref_from_xml(negdf, refnegdf)
negdf.to_csv("NegIDsC1_Ref.csv", index=False)

# poslibfile = "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-C1.msp"
# neglibfile = "MSDIAL-TandemMassSpectralAtlas-VS69-Neg-C1.msp"
# libfile = neglibfile
# outlibfile = libfile.replace(".msp", "_Ref.msp")
# df, outdata = add_ref_spectra(df, libfile)
#
# # Save outdata to new msp file
# save_as_msp(outdata, outlibfile, mode="w")
# # Save dfpos to new csv file
# df.to_csv(outcsvfile, index=False)




