import os
import pandas as pd
import numpy as np

def fix_header(df, hlen=5):
    # Merge multi-row headers into single row
    columns = df.columns
    new_columns = []
    for col in columns:
        # Get first hlen rows for this column
        header_parts = []
        for i in range(hlen):
            part = str(df.iloc[i][col]).strip()
            if part and part != 'nan':
                header_parts.append(part)
        new_col = "_".join(header_parts)
        new_columns.append(new_col)
    df.columns = new_columns
    # Drop the first hlen rows used for headers
    df = df.iloc[hlen:, :].reset_index(drop=True)
    # Drop completely empty columns
    df = df.dropna(axis=1, how='all')
    return df

def read_lipidblast_db(lbfile, classes=None):
    xl = pd.ExcelFile(lbfile)
    all_sheets = xl.sheet_names
    print("Sheets in Excel file:", all_sheets)

    data_frames = {}
    for sheet in all_sheets:
        if classes is not None and sheet not in classes:
            continue
        df = xl.parse(sheet)
        df = fix_header(df, hlen=5)
        data_frames[sheet] = df
    return data_frames

def match_ref_to_row(refmsms, lbrow, tol=0.01):
    # Convert refmsms string to array
    ref_array = np.array([m.split(":") for m in refmsms.strip().split(" ")]).astype(float)
    ref_mz = ref_array[:, 0]
    # ref_int = ref_array[:, 1]
    # ref_int_thresh = intthresh * np.max(ref_int)
    # b1 = ref_int >= ref_int_thresh
    # ref_mz = ref_mz[b1]

    matches = []
    # Find any in row that match ref_mz within tol
    for mz in ref_mz:
        # Search any in lbrow that matches mz within tol
        match = ""
        for i, val in enumerate(lbrow):
            try:
                lb_mz = float(val)
                if abs(lb_mz - mz) <= tol:
                    # Match is column name of librow at i
                    match = lbrow.index[i]
            except:
                continue
        if match:
            print("Matched mz", mz, "to", match)
        matches.append(match)
    return matches

def match_lipidblast_spectra(datadf, lipidblast_data, tol=0.01):
    matches = []
    for index, row in datadf.iterrows():
        lclass = row["Ontology"]
        # Find sheet that matches lclass
        if lclass not in lipidblast_data:
            matches.append("")
            continue
        lb_df = lipidblast_data[lclass]
        if "Name" not in lb_df.columns:
            print("Name column not found in LipidBlast data for class", lclass)
            matches.append("")
            continue

        name = row["Metabolite name"]
        # Split name and replace _ with -
        if "|" in name:
            name = name.split("|")[1]
        name = name.replace("_", "-")
        # print("Matching", row["Metabolite name"], "in class", lclass, "with", len(lb_df), "entries.")
        # Match name in lb_df
        matched_rows = lb_df[lb_df["Name"] == name]
        if matched_rows.empty:
            matched_rows = lb_df[lb_df["Short Name"] == name]
            if matched_rows.empty:
                matches.append("")
                continue

        print("Matched rows for", name, len(matched_rows))
        refmsms = row["Ref Spec"]
        adduct = row["Adduct type"]
        print("Adduct:", adduct)
        matchvalues = match_ref_to_row(refmsms, matched_rows.iloc[0], tol=tol)
        matches.append(", ".join([m for m in matchvalues if m]))

    datadf["LipidBlast_Matches"] = matches
    return datadf



datafile = r"C:\Python\pyRetip\notebooks\lipidomics\PosIDsC1_Ref.csv"
datadf = pd.read_csv(datafile)


lbfile = r"C:\Python\pyRetip\notebooks\lipidomics\LipidBlast DBs fork GlyceroLipids VS12\LipidBlast DBs fork GlyceroLipids VS12.xlsm"

print("Reading LipidBlast database...")
lipidblast_data = read_lipidblast_db(lbfile, classes=None)
print("Finished reading LipidBlast database.")
print(lipidblast_data.keys())

matched_df = match_lipidblast_spectra(datadf, lipidblast_data)
# Write to CSV
outputfile = r"C:\Python\pyRetip\notebooks\lipidomics\PosIDsC1_Ref_LipidBlastMatches.csv"
matched_df.to_csv(outputfile, index=False)
print("Wrote matched results to", outputfile)