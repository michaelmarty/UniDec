import pandas as pd


def keep_highest_area_adduct(df):
    df["Total Area"] = pd.to_numeric(df["Total Area"], errors="coerce")
    df["Precursor Mz"] = pd.to_numeric(df["Precursor Mz"], errors="coerce")
    df["Molecule Retention Time"] = pd.to_numeric(df["Molecule Retention Time"], errors="coerce")

    df = df.dropna(subset=[
        "Molecule",
        "Precursor Adduct",
        "Precursor Mz",
        "Molecule Retention Time",
        "Total Area"
    ])

    idx = df.groupby("Molecule")["Total Area"].idxmax()
    return df.loc[idx].reset_index(drop=True)


def combine_overlapping_mz_rt(df, mz_tolerance=0.3, rt_tolerance=0.3):
    """
    Combine rows ONLY if:
      - |Δm/z| <= mz_tolerance
      - |ΔRT| <= rt_tolerance (minutes)
    """
    df = df.rename(columns={
        "Molecule": "Compound",
        "Product Charge": "z",
        "Precursor Adduct": "Adduct",
        "Precursor Mz": "m/z",
        "Molecule Retention Time": "Retention Time (min)"
    })

    df = df.sort_values("m/z").reset_index(drop=True)
    used = set()
    combined_rows = []

    for i, row in df.iterrows():
        if i in used:
            continue

        mz_center = row["m/z"]
        rt_center = row["Retention Time (min)"]
        group_indices = [i]

        for j in range(i + 1, len(df)):
            if j in used:
                continue

            mz_diff = abs(df.loc[j, "m/z"] - mz_center)
            rt_diff = abs(df.loc[j, "Retention Time (min)"] - rt_center)

            if mz_diff <= mz_tolerance and rt_diff <= rt_tolerance:
                group_indices.append(j)

        used.update(group_indices)
        group = df.loc[group_indices]

        combined_rows.append({
            "Compound": " & ".join(sorted(set(group["Compound"].astype(str)))),
            "Retention Time (min)": group["Retention Time (min)"].median(),
            "z": group["z"].iloc[0],
            "Adduct": " & ".join(sorted(set(group["Adduct"].astype(str)))),
            "m/z": group["m/z"].median(),
            "Total Area": group["Total Area"].sum()
        })

    return pd.DataFrame(combined_rows)


# -----------------------
# USER INPUT
# -----------------------
input_csv = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Skyline Exports\New_FullHEK_PCs_LPCs.csv"
output_csv = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\Skyline Exports\Fixed_New_FullHEK_PCs_LPCs.csv"

mz_tolerance = 0.3      # m/z tolerance
rt_tolerance = 0.25    # RT tolerance in minutes (≈15 sec)

# -----------------------
# RUN
# -----------------------
df = pd.read_csv(input_csv)

df_filtered = keep_highest_area_adduct(df)
df_combined = combine_overlapping_mz_rt(
    df_filtered,
    mz_tolerance=mz_tolerance,
    rt_tolerance=rt_tolerance
)

df_combined = df_combined.sort_values("m/z").reset_index(drop=True)
df_combined.to_csv(output_csv, index=False)

print(
    f"Input rows: {len(df)}, "
    f"After highest-area adduct filter: {len(df_filtered)}, "
    f"After m/z + RT combining: {len(df_combined)}"
)
