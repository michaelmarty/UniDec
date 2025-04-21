import os
from unidec.modules.thermo_reader import ThermoImporter
import unidec.tools as ud
import pandas as pd

proteins = ["e protein", "vp4", "vpu"]
detergents = ["ldao", "c8e4", "dpc", "triton"]

def decode_file_name(path):
    protein = ""
    fname = str(path).lower()
    for p in proteins:
        if p in fname:
            protein = p

    detergent = ""
    for d in detergents:
        if d in fname:
            detergent = d

    return protein, detergent

# path = "Z:\Group Share\Julia\\viroporins in detergent\Vpu\detergent screen\C8E4\\230215 jat vpu c8e4 15uM.RAW"

topdir = "Z:\\Group Share\\Julia\\viroporins in detergent"
os.chdir(topdir)

files = ud.match_files_recursive(topdir, ".RAW")

df = pd.DataFrame(files, columns=["Path"])

print(files)


header_queries = ["Extended Trapping"]

if False:

    for i, f in enumerate(files):
        TDI = ThermoImporter.ThermoDataImporter(f)

        p, d = decode_file_name(f)
        df.loc[i, "Protein"] = p
        df.loc[i, "Detergent"] = d

        sourcefrag = TDI.get_sid_voltage()
        df.loc[i, "Source Fragmentation"] = sourcefrag

        for h in header_queries:
            try:
                value = TDI.msrun.header[h]
                print(h, value)
                # Save back to df

                df.loc[i, h] = float(value)
            except:
                pass

        try:
            df.loc[i, "Pressure"] = float(TDI.msrun.uhv_pressure)
        except:
            pass

    # Drop rows where pressure is ""
    df = df[df["Pressure"] != ""]
    # Drop rows where detergent and protein are ""
    df = df[df["Detergent"] != ""]
    df = df[df["Protein"] != ""]
    # Drop rows where pressure is nan
    df = df[df["Pressure"].notna()]

    # sort by protein first and then detergent
    df = df.sort_values(by=["Protein", "Detergent"])

    df.to_excel("Params.xlsx")
else:
    df = pd.read_excel("Params.xlsx")

    #print min and max of each param
    for c in df.columns:
        if c not in ["Path", "Protein", "Detergent"]:
            print(c, df[c].min(), df[c].max(), df[c].mean(), df[c].mode(), df[c].median())


