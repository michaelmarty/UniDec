import pandas as pd

# ==========================
# USER INPUTS
# ==========================
transition_list_path = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Combined_TransitionList_Clean.csv"
isolation_list_path = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\HEK_IsoList.csv"
output_path = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Combined_Cleantransition_list_updated_RT.csv"

# ==========================
# LOAD FILES
# ==========================
trans = pd.read_csv(transition_list_path)
iso = pd.read_csv(isolation_list_path)

# Clean names
trans["Molecule Name"] = trans["Molecule Name"].astype(str).str.strip()
iso["Compound"] = iso["Compound"].astype(str).str.strip()

# Extract base lipid name (before "|") for BOTH lists
trans["Base Lipid"] = trans["Molecule Name"].str.split("|").str[0]
iso["Base Lipid"] = iso["Compound"].str.split("|").str[0]

# Ensure numeric
iso["Retention Time (min)"] = pd.to_numeric(
    iso["Retention Time (min)"], errors="coerce"
)
iso["RT Window (min)"] = pd.to_numeric(
    iso["RT Window (min)"], errors="coerce"
)

# ==========================
# BUILD LOOKUP
# ==========================
rt_lookup = (
    iso
    .dropna(subset=["Base Lipid", "Retention Time (min)"])
    .set_index("Base Lipid")[["Retention Time (min)", "RT Window (min)"]]
)

# ==========================
# APPLY RTs
# ==========================
updated = 0

for compound, rt_vals in rt_lookup.iterrows():
    mask = trans["Base Lipid"] == compound

    if mask.any():
        trans.loc[mask, "Explicit Retention Time"] = rt_vals["Retention Time (min)"]
        trans.loc[mask, "Explicit Retention Time Window"] = rt_vals["RT Window (min)"]
        updated += mask.sum()
    else:
        print(f"⚠️ No transitions found for: {compound}")

# ==========================
# SAVE
# ==========================
trans.to_csv(output_path, index=False)

print(f"\nUpdated {updated} transition rows")
print(f"Saved to: {output_path}")

# ==========================
# CHECK MISSING
# ==========================
trans_names = set(trans["Base Lipid"].str.lower().str.strip())
iso_names = set(iso["Base Lipid"].str.lower().str.strip())

missing = iso_names - trans_names
print("Compounds in isolation list not found in transition list:")
print(missing)
