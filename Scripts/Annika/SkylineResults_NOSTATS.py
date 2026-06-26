import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# USER SETTINGS
# ============================================================

# Input spreadsheet.
# This can be .csv or .xlsx
input_file = r"Z:\Group Share\Annika\Stellar\EC\B2\Molecule Transition Results.csv"

# Output files
output_table = "fold_change_results.csv"
output_plot = "fold_change_plot.pdf"

# Your replicate names exactly as they appear in the "Replicate Name" column
control_replicate = "HEK Extr"
analyte_replicate = "B2 Discs"

# Transition mode:
# "Products" = exclude precursor transitions
# "Precursors" = only precursor transitions
# "All" = use all transitions
transition_mode = "Products"

# Whether to normalize each lipid by its internal standard first
normalize_by_internal_standard = True

# Whether to remove internal standards after IS normalization
drop_internal_standards = True

# Whether to normalize each sample by total signal after IS normalization
normalize_by_mtic = True

# Optional: only plot the top N lipids by absolute log2 fold change.
# Set to None to plot everything.
top_n_to_plot = 40


# ============================================================
# DATA LOADING
# ============================================================

def load_input_file(path):
    """
    Load a CSV or Excel file.
    """
    if path.lower().endswith(".csv"):
        return pd.read_csv(path)
    elif path.lower().endswith(".xlsx") or path.lower().endswith(".xls"):
        return pd.read_excel(path)
    else:
        raise ValueError("Input file must be .csv, .xlsx, or .xls")


# ============================================================
# CLEANING FUNCTIONS
# ============================================================

def clean_classes(df):
    """
    Clean common lipid class naming differences.
    """
    df = df.copy()

    df["Molecule List Name"] = df["Molecule List Name"].astype(str)

    df["Molecule List Name"] = df["Molecule List Name"].str.replace("DAG", "DG", regex=False)
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("TAG", "TG", regex=False)
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("SE 27:1", "CE", regex=False)
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("SHexCer", "HexCer", regex=False)
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("Hex2Cer", "HexCer", regex=False)

    return df


def clean_adducts(df):
    """
    Clean common labeled internal-standard adduct names so they match unlabeled forms.
    """
    df = df.copy()

    adduct_col = "Precursor Adduct" if "Precursor Adduct" in df.columns else "Adduct"

    replacements = {
        "[M7H2+H]": "[M+H]",
        "[M7H2-H]": "[M-H]",
        "[M5H2+H]": "[M+H]",
        "[M5H2-H]": "[M-H]",
        "[M5H2+CH3COO]": "[M+CH3COO]",
        "[M5H2-CH3COO]": "[M-CH3COO]",
        "[M9H2+H]": "[M+H]",
        "[M9H2-H]": "[M-H]",
        "[M7H2+H-H2O]": "[M+H-H2O]",
        "[M9H2+CH3COO]": "[M+CH3COO]",
        "[M9H2-CH3COO]": "[M-CH3COO]",
    }

    df[adduct_col] = df[adduct_col].replace(replacements)

    return df


# ============================================================
# TRANSITION SUMMING
# ============================================================

def sum_transitions(df, mode="Products"):
    """
    Sum transition areas so there is one row per molecule/adduct,
    with one column per replicate.

    Expected input columns:
        Molecule List Name
        Molecule
        Replicate Name
        Area
        Adduct or Precursor Adduct
        Fragment Ion
        RT or Retention Time
        Precursor Mz, optional
    """

    required_cols = [
        "Molecule List Name",
        "Molecule",
        "Replicate Name",
        "Area",
    ]

    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    adduct_col = "Precursor Adduct" if "Precursor Adduct" in df.columns else "Adduct"
    if adduct_col not in df.columns:
        raise ValueError("Missing required column: Adduct or Precursor Adduct")

    rt_col = "RT" if "RT" in df.columns else "Retention Time"
    if rt_col not in df.columns:
        raise ValueError("Missing required column: RT or Retention Time")

    if "Fragment Ion" not in df.columns:
        raise ValueError("Missing required column: Fragment Ion")

    df = df.copy()
    df = clean_classes(df)
    df = clean_adducts(df)

    df["Area"] = pd.to_numeric(df["Area"], errors="coerce").fillna(0)

    replicates = list(df["Replicate Name"].dropna().unique())

    rows = []

    group_cols = ["Molecule List Name", "Molecule", adduct_col]

    for group_values, subdf in df.groupby(group_cols):
        class_name, molecule, adduct = group_values

        if mode == "Products":
            use_df = subdf[
                ~subdf["Fragment Ion"].astype(str).str.contains(
                    "precursor", case=False, na=False
                )
            ]
        elif mode == "Precursors":
            use_df = subdf[
                subdf["Fragment Ion"].astype(str).str.contains(
                    "precursor", case=False, na=False
                )
            ]
        elif mode == "All":
            use_df = subdf
        else:
            raise ValueError("mode must be 'Products', 'Precursors', or 'All'")

        newrow = {
            "Molecule List Name": class_name,
            "Molecule": molecule,
            "Adduct": adduct,
            "RT": subdf[rt_col].mean(),
        }

        if "Precursor Mz" in subdf.columns:
            newrow["MZ"] = subdf["Precursor Mz"].mean()

        for rep in replicates:
            rep_sum = use_df.loc[use_df["Replicate Name"] == rep, "Area"].sum()
            newrow[rep] = rep_sum

        rows.append(newrow)

    outdf = pd.DataFrame(rows)

    return outdf, replicates


# ============================================================
# INTERNAL STANDARD NORMALIZATION
# ============================================================

def get_internal_standard_rows(df, class_name, adduct):
    """
    Find internal standards for a given lipid class and adduct.

    Internal standards are assumed to have Molecule names starting with "IS ".
    """
    subdf = df[
        (df["Molecule List Name"] == class_name)
        & (df["Adduct"] == adduct)
    ]

    is_rows = subdf[subdf["Molecule"].astype(str).str.startswith("IS ")]

    return is_rows


def normalize_by_is(df, replicates):
    """
    Normalize each lipid by its matching internal standard.

    Matching is done by:
        1. Same Molecule List Name
        2. Same Adduct
        3. Closest RT if multiple internal standards are available
    """

    normdf = df.copy()

    for rep in replicates:
        normdf[rep] = pd.to_numeric(normdf[rep], errors="coerce").fillna(0)

    for idx, row in normdf.iterrows():
        molecule = str(row["Molecule"])

        # Do not normalize IS rows here; they will be dropped later if desired.
        if molecule.startswith("IS "):
            continue

        class_name = row["Molecule List Name"]
        adduct = row["Adduct"]
        rt = row["RT"]

        is_rows = get_internal_standard_rows(normdf, class_name, adduct)

        if len(is_rows) == 0:
            print(
                f"WARNING: No internal standard found for "
                f"{molecule}, class {class_name}, adduct {adduct}. "
                f"Setting values to 0."
            )
            for rep in replicates:
                normdf.at[idx, rep] = 0
            continue

        # If multiple internal standards exist, choose the closest RT.
        if len(is_rows) > 1:
            is_rows = is_rows.copy()
            is_rows["RT_Diff"] = (is_rows["RT"] - rt).abs()
            is_row = is_rows.sort_values("RT_Diff").iloc[0]
        else:
            is_row = is_rows.iloc[0]

        for rep in replicates:
            is_value = float(is_row[rep])
            lipid_value = float(row[rep])

            if is_value <= 0:
                print(
                    f"WARNING: Internal standard value is zero for "
                    f"{molecule}, replicate {rep}. Setting value to 0."
                )
                normdf.at[idx, rep] = 0
            else:
                normdf.at[idx, rep] = lipid_value / is_value

    return normdf


# ============================================================
# TIC / mTIC NORMALIZATION
# ============================================================

def normalize_tic(df, replicates):
    """
    Normalize each replicate column by its total signal.
    """
    normdf = df.copy()

    for rep in replicates:
        total = normdf[rep].sum()

        if total <= 0:
            print(f"WARNING: Total signal for {rep} is zero. Leaving unchanged.")
            continue

        normdf[rep] = normdf[rep] / total

    return normdf


# ============================================================
# FOLD CHANGE CALCULATION
# ============================================================

def calculate_fold_change(df, control_col, analyte_col):
    """
    Calculate fold change and log2 fold change for one control and one analyte.
    No pseudocount is used.
    """

    outdf = df.copy()

    outdf[control_col] = pd.to_numeric(outdf[control_col], errors="coerce")
    outdf[analyte_col] = pd.to_numeric(outdf[analyte_col], errors="coerce")

    outdf = outdf.dropna(subset=[control_col, analyte_col])

    # Keep only positive values because log2 fold change requires positive data.
    outdf = outdf[
        (outdf[control_col] > 0)
        & (outdf[analyte_col] > 0)
    ].copy()

    outdf["Control"] = outdf[control_col]
    outdf["Analyte"] = outdf[analyte_col]

    outdf["FC"] = outdf["Analyte"] / outdf["Control"]
    outdf["Log2FC"] = np.log2(outdf["FC"])

    outdf = outdf.replace([np.inf, -np.inf], np.nan)
    outdf = outdf.dropna(subset=["FC", "Log2FC"])

    return outdf


# ============================================================
# PLOTTING
# ============================================================

def make_fold_change_plot(
    fcdf,
    molecule_col="Molecule",
    fc_col="Log2FC",
    title="Analyte vs Control Fold Change",
    top_n=None,
    output_file=None,
):
    """
    Make a simple log2 fold-change bar plot.
    """

    plotdf = fcdf.copy()

    if top_n is not None:
        plotdf = plotdf.reindex(
            plotdf[fc_col].abs().sort_values(ascending=False).index
        ).head(top_n)

    plotdf = plotdf.sort_values(fc_col)

    plt.figure(figsize=(12, 6))
    plt.bar(plotdf[molecule_col], plotdf[fc_col])

    plt.axhline(0, linestyle="--", linewidth=1)
    plt.ylabel("Log2 Fold Change")
    plt.xlabel("Molecule")
    plt.title(title)

    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()

    if output_file is not None:
        plt.savefig(output_file, bbox_inches="tight")

    plt.show()


# ============================================================
# MAIN SCRIPT
# ============================================================

def main():
    print("Loading input file...")
    rawdf = load_input_file(input_file)

    print("Summing transitions...")
    summed_df, replicates = sum_transitions(
        rawdf,
        mode=transition_mode
    )

    print("Replicates found:")
    print(replicates)

    if control_replicate not in replicates:
        raise ValueError(
            f"Control replicate '{control_replicate}' was not found. "
            f"Available replicates are: {replicates}"
        )

    if analyte_replicate not in replicates:
        raise ValueError(
            f"Analyte replicate '{analyte_replicate}' was not found. "
            f"Available replicates are: {replicates}"
        )

    normdf = summed_df.copy()

    if normalize_by_internal_standard:
        print("Normalizing by internal standard...")
        normdf = normalize_by_is(normdf, replicates)

    if drop_internal_standards:
        print("Dropping internal standards...")
        normdf = normdf[
            ~normdf["Molecule"].astype(str).str.startswith("IS ")
        ].copy()

    if normalize_by_mtic:
        print("Normalizing by TIC/mTIC...")
        normdf = normalize_tic(normdf, replicates)

    print("Calculating fold change...")
    fcdf = calculate_fold_change(
        normdf,
        control_col=control_replicate,
        analyte_col=analyte_replicate
    )

    print(f"Saving results table to {output_table}...")
    fcdf.to_csv(output_table, index=False)

    print(f"Saving plot to {output_plot}...")
    make_fold_change_plot(
        fcdf,
        title=f"{analyte_replicate} vs {control_replicate}",
        top_n=top_n_to_plot,
        output_file=output_plot
    )

    print("Done.")


if __name__ == "__main__":
    main()