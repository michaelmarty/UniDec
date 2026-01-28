import pandas as pd

# ----------------------------
# USER INPUTS
# ----------------------------
input_xlsx = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\NewColumn_Areas_1,2,3,4,5,6.xlsx"
output_xlsx = r"Z:\Group Share\Annika\Stellar\EC\AqpZ New Data\Normalized_NewColumn_Areas_1,2,3,4,5,6.xlsx"

# Class → Internal standard name (EXACT match to Row Labels)
CLASS_IS = {
    "CL": "CL 14:0-14:0-14:0-14:0",
    "PE": "PE 15:0-18:1(+[2]H7)",
    "PG": "PG 15:0-18:1(+[2]H7)"
}

# ----------------------------
# LOAD WORKBOOK
# ----------------------------
xls = pd.ExcelFile(input_xlsx)
writer = pd.ExcelWriter(output_xlsx, engine="openpyxl")

# ----------------------------
# PROCESS SHEETS
# ----------------------------
for sheet in xls.sheet_names:
    df = pd.read_excel(xls, sheet_name=sheet)

    if "Row Labels" not in df.columns or len(df.columns) < 2:
        df.to_excel(writer, sheet_name=sheet, index=False)
        continue

    area_col = df.columns[1]
    df[area_col] = pd.to_numeric(df[area_col], errors="coerce")

    # Extract lipid class from Row Labels (text before first space)
    df["Class"] = df["Row Labels"].str.split().str[0]

    # Pre-calculate IS areas for this sheet
    is_area_by_class = {}

    for lipid_class, is_name in CLASS_IS.items():
        match = df.loc[df["Row Labels"] == is_name, area_col]
        if not match.empty:
            is_area_by_class[lipid_class] = match.values[0]
        else:
            print(f"⚠️ {lipid_class} IS not found in sheet {sheet}")

    # ----------------------------
    # CLASS-SPECIFIC IS NORMALIZATION
    # ----------------------------
    def normalize_row(row):
        cls = row["Class"]
        if cls in is_area_by_class and pd.notna(is_area_by_class[cls]):
            return row[area_col] / is_area_by_class[cls]
        return float("nan")

    df["IS_normalized_class"] = df.apply(normalize_row, axis=1)

    # ----------------------------
    # TOTAL NORMALIZATION (AFTER IS)
    # ----------------------------
    total = df["IS_normalized_class"].sum(skipna=True)
    df["IS_then_total_norm_%"] = (df["IS_normalized_class"] / total) * 100

    # Save
    df.to_excel(writer, sheet_name=sheet, index=False)

# ----------------------------
# FINALIZE
# ----------------------------
writer.close()
print("✅ Class-specific IS normalization complete")
