import pandas as pd

# ===== FILE PATHS =====
file1 = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Combined_IsolationList2.csv"
file2 = r"list2.csv"
output = "combined_skyline_list.csv"

# ===== READ FILES =====
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Optional: enforce consistent column order/names
columns = [
    "Compound",
    "m/z",
    "z",
    "Polarity",
    "Retention Time (min)",
    "RT Window (min)"
]

df1 = df1[columns]
df2 = df2[columns]

# ===== COMBINE + REMOVE DUPLICATES =====
combined = pd.concat([df1, df2], ignore_index=True)
combined = combined.drop_duplicates()

# ===== SAVE RESULT =====
combined.to_csv(output, index=False)

print(f"Combined list saved to: {output}")
