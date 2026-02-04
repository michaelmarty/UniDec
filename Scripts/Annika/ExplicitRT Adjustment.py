import pandas as pd
import re

# ==========================
# USER INPUTS
# ==========================
transition_list_path = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Peak Masters & xmls\Combined_TransitionList_Clean.csv"
isolation_list_path  = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Peak Masters & xmls\AdjustedERTs_IsoList_All.csv"
output_path          = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Peak Masters & xmls\Combined_Clean_AdjustedRT_TransList.csv"

# ==========================
# LOAD FILES
# ==========================
trans = pd.read_csv(transition_list_path)
iso   = pd.read_csv(isolation_list_path)

# ==========================
# CLEAN NAMES
# ==========================
trans["Molecule Name"] = trans["Molecule Name"].astype(str).str.strip()
iso["Molecule Name"]   = iso["Molecule Name"].astype(str).str.strip()

# ==========================
# EXTRACT BASE LIPID
# (LEFT of | + optional _A/_B at END)
# ==========================
def extract_base_lipid(name):
    name = str(name).strip()

    # Left side of pipe
    base = name.split("|", 1)[0].strip()

    # Isomer suffix at very end (_A, _B, etc)
    m = re.search(r'(_[A-Z])\s*$', name)
    if m:
        base = f"{base}{m.group(1)}"

    return base

# MUST happen BEFORE any dropna / lookup
trans["Base Lipid"] = trans["Molecule Name"].apply(extract_base_lipid)
iso["Base Lipid"]   = iso["Molecule Name"].apply(extract_base_lipid)

# ==========================
# ENSURE NUMERIC
# ==========================
iso["Explicit Retention Time"] = pd.to_numeric(
    iso["Explicit Retention Time"], errors="coerce"
)
iso["Explicit Retention Time Window"] = pd.to_numeric(
    iso["Explicit Retention Time Window"], errors="coerce"
)

# ==========================
# BUILD LOOKUP
# ==========================
rt_lookup = (
    iso
    .dropna(subset=["Base Lipid", "Explicit Retention Time"])
    .set_index("Base Lipid")[["Explicit Retention Time", "Explicit Retention Time Window"]]
)

# ==========================
# APPLY RTs
# ==========================
updated = 0

for compound, rt_vals in rt_lookup.iterrows():
    mask = trans["Base Lipid"] == compound

    if mask.any():
        trans.loc[mask, "Explicit Retention Time"] = rt_vals["Explicit Retention Time"]
        trans.loc[mask, "Explicit Retention Time Window"] = rt_vals["Explicit Retention Time Window"]
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
# CHECK MISSING COMPOUNDS
# ==========================
trans_names = set(trans["Base Lipid"].str.lower().str.strip())
iso_names   = set(iso["Base Lipid"].str.lower().str.strip())

missing = iso_names - trans_names
print("\nCompounds in isolation list not found in transition list:")
print(missing)

# ==========================
# SANITY CHECK: ISOMERS
# ==========================
def strip_isomer(name):
    return re.sub(r'_[A-Z]$', '', name)

trans_iso_counts = (
    trans["Base Lipid"]
    .dropna()
    .apply(strip_isomer)
    .value_counts()
)

iso_iso_counts = (
    iso["Base Lipid"]
    .dropna()
    .apply(strip_isomer)
    .value_counts()
)

problem_lipids = []

for lipid, trans_count in trans_iso_counts.items():
    iso_count = iso_iso_counts.get(lipid, 0)

    if trans_count > 1 and iso_count <= 1:
        problem_lipids.append((lipid, trans_count, iso_count))

if problem_lipids:
    print("\n⚠️ Potential isomer RT collisions detected:")
    for lipid, t_cnt, i_cnt in problem_lipids:
        print(
            f"  {lipid}: {t_cnt} transition isomers "
            f"but {i_cnt} isolation entries"
        )
else:
    print("\n✅ No isomer RT conflicts detected")
