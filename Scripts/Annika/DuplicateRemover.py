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
        by=["m/z", "Retention Time (min)"]
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

            mz_close = abs(row["m/z"] - other["m/z"]) <= mz_tolerance
            rt_close = abs(row["Retention Time (min)"] - other["Retention Time (min)"]) <= rt_tolerance
            same_charge = row["z"] == other["z"]


            if mz_close and rt_close and same_charge:
                group_indices.append(j)

        used.update(group_indices)
        group = df.loc[group_indices]

        combined_rows.append({
            "Compound": " & ".join(sorted(group["Compound"].astype(str).unique())),
            "z": group["z"].iloc[0],
            "m/z": group["m/z"].mean(),
            "Retention Time (min)": group["Retention Time (min)"].mean()
        })

    out_df = pd.DataFrame(combined_rows)


    out_df.to_csv(output_csv, index=False)

    print(f"Combined {len(df)} rows → {len(out_df)} rows")
    print(f"Saved to: {output_csv}")


if __name__ == "__main__":
    input_csv = r"Z:\Group Share\Annika\IQ-X\AcquireX and Neo\AcquireX\Combined_IsolationList_AcquireX.csv"
    output_csv = r"Z:\Group Share\Annika\IQ-X\AcquireX and Neo\AcquireX\NoDups_Combined_IsolationList_AcquireX.csv"

    combine_overlapping_molecules(
        input_csv,
        output_csv,
        mz_tolerance=0.3,
        rt_tolerance=0.7
    )
