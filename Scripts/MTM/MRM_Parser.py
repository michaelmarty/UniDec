import os
import pymzml
import pandas as pd

file = "C:\\Users\\mm96978\\Downloads\\UT-Austin_Dr Marty_7600_7500+_Data_May2025\\UT-Austin_Dr Marty_7600_7500+_Data_May2025\\20250304_7500_Samples\\7500 Pos LX MRM_Samples-H Pos.mzML"

d = pymzml.run.Reader(file, skip_chromatogram=False)

numc = d.get_chromatogram_count()

datalist = []
for c in d:
    s = c.ID
    if "SRM" in s:
        q1 = 0
        q3 = 0
        rtstart = 0
        rtend = 0
        rtmid = 0
        rtwidth = 0
        name = ""
        type = ""


        splits = s.split(" ")
        for split in splits:
            if "Q1" in split:
                q1 = float(split.split("=")[1])
            if "Q3" in split:
                q3 = float(split.split("=")[1])
            if "start" in split:
                rtstart = float(split.split("=")[1])
            if "end" in split:
                rtend = float(split.split("=")[1])
            if "name" in split:
                name = split.split("=")[1]
                type = name.split("(")[0]


        rtmid = (rtend + rtstart)/2.
        rtwidth = rtend - rtstart

        # Save to row on PD dataframe
        # row = {
        #     "Name": name,
        #     "RT Midpoint": rtmid,
        #     "RT Width": rtwidth,
        #     "Q1": q1,
        #     "Q3": q3,
        #     "RT Start": rtstart,
        #     "RT End": rtend,
        # }

        row = {
            "Molecule Name": name,
            "Molecule List Name": type,
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
df.to_csv(os.path.join(os.path.dirname(file), "MRM_TransitionListPos.csv"), index=False)






