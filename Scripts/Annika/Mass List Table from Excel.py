import pandas as pd

# ===== FILE PATHS =====
input_file = "Z:\Group Share\Annika\Stellar\IS_Information.xlsx"   # your Excel file
sheet_name = "UltimateSPLASH One"              # <-- CHANGE to the sheet you want
output_csv = "Z:\\Group Share\\Annika\\Stellar\\UltimateSPLASH_One_masslist.csv"

# ===== LOAD DATA FROM SPECIFIC SHEET =====
df = pd.read_excel(input_file, sheet_name=sheet_name)

# Clean column names just in case
df.columns = df.columns.str.strip()

output_rows = []

for _, row in df.iterrows():
    compound_name = f"{row['Class']} {row['Lipid']}"

    # Positive polarity entry
    if pd.notna(row['m/z (+)']):
        output_rows.append({
            "Compound": compound_name,
            "m/z": row['m/z (+)'],
            "z": 1,
            "Polarity": "Positive",
            "Retention Time (min)": "",
            "RT Window (min)": ""
        })

    # Negative polarity entry
    if pd.notna(row['m/z (-)']):
        output_rows.append({
            "Compound": compound_name,
            "m/z": row['m/z (-)'],
            "z": 1,
            "Polarity": "Negative",
            "Retention Time (min)": "",
            "RT Window (min)": ""
        })

# ===== SAVE OUTPUT =====
out_df = pd.DataFrame(output_rows)
out_df.to_csv(output_csv, index=False)

print(f"Mass list saved to: {output_csv}")
