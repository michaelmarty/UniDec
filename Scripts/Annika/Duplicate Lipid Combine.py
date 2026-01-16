import pandas as pd
import numpy as np

# Load the CSV
file = "Z:\\Group Share\\Annika\\Stellar\\Untargeted DDA\\HEK\\Stellar_MassList_HEKPos2.csv"
df = pd.read_csv(file)

# --------------------------------------------------
# Step 1: Condense true duplicates (Compound + m/z exact)
# --------------------------------------------------
df_dedup = (
    df.groupby(["Compound", "m/z"], as_index=False)
      .agg({
          "Product Charge": "first",
          "Explicit Retention Time Window": "max",
          "Explicit Retention Time": "first",
          "Precursor Adduct": "first"
      })
)

# --------------------------------------------------
# Step 2: Fuzzy merge compounds within ±1 m/z
# --------------------------------------------------
def fuzzy_merge(group, mz_col="m/z", tolerance=1.0):
    """
    Merge rows with the same Compound within ±tolerance m/z
    """
    group = group.sort_values(mz_col).reset_index(drop=True)
    merged_rows = []
    i = 0
    while i < len(group):
        # Start a cluster
        cluster = [group.loc[i]]
        j = i + 1
        while j < len(group) and abs(group.loc[j, mz_col] - group.loc[i, mz_col]) <= tolerance:
            cluster.append(group.loc[j])
            j += 1
        # Merge cluster
        merged = {
            "Explicit Retention Time Window": max([row["Explicit Retention Time Window"] for row in cluster]),
            "Explicit Retention Time": np.mean([row["Explicit Retention Time"] for row in cluster]),
            "Product Charge": cluster[0]["Product Charge"],  # keep first
            "Compound": cluster[0]["Compound"],  # same compound
            "m/z": np.mean([row["m/z"] for row in cluster]),  # avg m/z
            "Precursor Adduct": cluster[0]["Precursor Adduct"]  # keep first
        }
        merged_rows.append(merged)
        i = j
    return pd.DataFrame(merged_rows)

# Apply fuzzy merge per compound
df_fuzzy = pd.concat([fuzzy_merge(group) for _, group in df_dedup.groupby("Compound")], ignore_index=True)

# --------------------------------------------------
# Step 3: Optional combine all compounds with identical m/z (as before)
# --------------------------------------------------
df_combined = (
    df_fuzzy.groupby("m/z", as_index=False)
            .agg({
                "Compound": lambda x: " \\ ".join(sorted(set(x))),
                "Product Charge": "first",
                "Explicit Retention Time Window": "max",
                "Explicit Retention Time": "first",
                "Precursor Adduct": "first"
            })
)

# Reorder columns to match Stellar expectations
df_combined = df_combined[
    [
        "Explicit Retention Time Window",
        "Explicit Retention Time",
        "Product Charge",
        "Compound",
        "m/z",
        "Precursor Adduct"
    ]
]

# Save output
output_path = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\\HEK\stellar_mass_table_unique_mz_fuzzy2.csv"
df_combined.to_csv(output_path, index=False)

print(f"Saved fuzzy-merged mass table to:\n{output_path}")
