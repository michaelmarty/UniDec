import os
import pandas as pd

resfile = "C:\\Users\\mm96978\\Downloads\\UT-Austin_Dr Marty_7600_7500+_Data_May2025\\UT-Austin_Dr Marty_7600_7500+_Data_May2025\\Sample H CID Matches.xlsx"
df = pd.read_excel(resfile)


datalist = []
for i, r in df.iterrows():
    name = r["Name"]
    ltype = r["Ontology"]
    rtmid = r["RT (min)"]

    rtstart = r["RT left(min)"]
    rtend = r["RT right (min)"]
    rtwidth = rtend - rtstart

    q1 = r["Precursor m/z"]

    ms2spec = r["MSMS spectrum"]
    ms2data = []
    for mzi in ms2spec.split(" "):
        mz, i = mzi.split(":")
        ms2data.append((float(mz), float(i)))

    # Find highest intensity m/z in ms2data that isn't the q1 m/z
    q3 = 0
    max_intensity = 0
    for mz, intensity in ms2data:
        if mz != q1 and intensity > max_intensity:
            max_intensity = intensity
            q3 = mz



    row = {
        "Molecule Name": name,
        "Molecule List Name": ltype,
        "Explicit Retention Time": rtmid,
        "Explicit Retention Time Window": rtwidth,
        "Precursor m/z": q1,
        "Product m/z": q3,
        "Precursor Charge": 1,  # Assuming charge state is 1, adjust if needed
        "Product Charge": 1,  # Assuming charge state is 1, adjust if needed
    }
    # Add row to DataFrame
    datalist.append(row)

# Convert list of dictionaries to DataFrame
df = pd.DataFrame(datalist)

# Save to CSV
df.to_csv(os.path.join(os.path.dirname(resfile), "MRMHR_TransitionListPos.csv"), index=False)






