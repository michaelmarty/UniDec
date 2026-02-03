import pandas as pd

def combine_overlapping_molecules(
    input_csv,
    output_csv,
    mz_tolerance=0.3,
    rt_tolerance=0.3
):
    df = pd.read_csv(input_csv)

    # Clean column names
    df.columns = df.columns.str.strip()

    # Sort so grouping is deterministic
    df = df.sort_values(
        by=["Precursor Mz", "Molecule Retention Time"]
    ).reset_index(drop=True)

    combined_rows = []
    used = set()

    for i, row in df.iterrows():
        if i in used:
            continue

        group_indices = [i]

        for j in range(i + 1, len(df)):
            if j in used:
                continue

            other = df.loc[j]

            mz_close = abs(row["Precursor Mz"] - other["Precursor Mz"]) <= mz_tolerance
            rt_close = abs(row["Molecule Retention Time"] - other["Molecule Retention Time"]) <= rt_tolerance
            same_charge = row["Product Charge"] == other["Product Charge"]
            same_adduct = row["Precursor Adduct"] == other["Precursor Adduct"]

            if mz_close and rt_close and same_charge and same_adduct:
                group_indices.append(j)

        used.update(group_indices)
        group = df.loc[group_indices]

        combined_rows.append({
            "Molecule": " & ".join(sorted(group["Molecule"].astype(str).unique())),
            "Product Charge": group["Product Charge"].iloc[0],
            "Precursor Adduct": group["Precursor Adduct"].iloc[0],
            "Precursor Mz": group["Precursor Mz"].mean(),
            "Molecule Retention Time": group["Molecule Retention Time"].mean()
        })

    out_df = pd.DataFrame(combined_rows)

    # =========================
    # ADDITIONS START HERE
    # =========================

    RT_WINDOW_VALUE = 2.0  # <-- set your desired RT window (min)

    out_df = out_df.rename(columns={
        "Molecule": "Compound",
        "Product Charge": "z",
        "Precursor Mz": "m/z",
        "Molecule Retention Time": "Retention Time (min)"
    })

    out_df["RT Window (min)"] = RT_WINDOW_VALUE

    # =========================
    # ADDITIONS END HERE
    # =========================

    out_df.to_csv(output_csv, index=False)

    print(f"Combined {len(df)} rows → {len(out_df)} rows")
    print(f"Saved to: {output_csv}")


if __name__ == "__main__":
    input_csv = r"Z:\Group Share\Annika\Stellar\HEK\tMS2\HEK_IsoList.csv"
    output_csv = r"Z:\\Group Share\\Annika\\Stellar\\HEK\\tMS2\\NoDups_HEK_IsoList.csv"

    combine_overlapping_molecules(
        input_csv,
        output_csv,
        mz_tolerance=0.3,
        rt_tolerance=0.7
    )
