"""Export 90-100 kDa E. coli proteins ordered by residue content."""

from pathlib import Path

import pandas as pd


FASTA_DIR = Path(r"C:\Data\IsoNN\fastas")
INPUT_PATH = FASTA_DIR / "ecoli.xlsx"
OUTPUT_PATH = FASTA_DIR / "ecoli_filter.xlsx"
CYSTEINE_OUTPUT_PATH = FASTA_DIR / "ecoli_filter_cysteine.xlsx"
MIN_MASS_DA = 90_000
MAX_MASS_DA = 100_000
OUTER_MEMBRANE_TERM = "cell outer membrane [GO:0009279]"


def add_residue_content(proteins, residue, residue_name):
    """Add count and percent columns for a residue to protein records."""
    sequences = proteins["Sequence"].fillna("").astype(str).str.upper()
    proteins = proteins.copy()
    proteins[residue_name] = sequences.str.count(residue)
    proteins[f"{residue_name} (%)"] = 100 * proteins[residue_name] / sequences.str.len()
    return proteins


def add_methionine_content(proteins):
    """Add methionine count and percent columns to protein records."""
    return add_residue_content(proteins, "M", "Methionine")


def add_cysteine_content(proteins):
    """Add cysteine count and percent columns to protein records."""
    return add_residue_content(proteins, "C", "Cysteine")


def filter_by_mass_and_residue(
    input_path, output_path, residue, residue_name
):
    """Write 90-100 kDa proteins with a residue content in ascending order."""
    proteins = pd.read_excel(input_path)
    required_columns = {
        "Mass",
        "Sequence",
        "Gene Ontology (cellular component)",
    }
    missing_columns = required_columns.difference(proteins.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Input workbook is missing required columns: {missing}")

    proteins = add_residue_content(proteins, residue, residue_name)
    filtered = proteins.loc[
        proteins["Mass"].between(MIN_MASS_DA, MAX_MASS_DA, inclusive="both")
    ].copy()
    filtered = filtered.sort_values(
        [f"{residue_name} (%)", residue_name, "Mass"],
        kind="stable",
        na_position="last",
    )
    filtered.to_excel(output_path, index=False)
    return filtered


def filter_by_mass_and_methionine(input_path=INPUT_PATH, output_path=OUTPUT_PATH):
    """Write 90-100 kDa proteins with methionine content in ascending order."""
    return filter_by_mass_and_residue(
        input_path, output_path, "M", "Methionine"
    )


def filter_by_mass_and_cysteine(
    input_path=INPUT_PATH, output_path=CYSTEINE_OUTPUT_PATH
):
    """Write 90-100 kDa proteins with cysteine content in ascending order."""
    return filter_by_mass_and_residue(input_path, output_path, "C", "Cysteine")


def average_residue_content(proteins, residue_name):
    """Return average residue percentage for all and outer-membrane proteins."""
    whole_dataset_average = proteins[f"{residue_name} (%)"].mean()
    outer_membrane_proteins = proteins.loc[
        proteins["Gene Ontology (cellular component)"]
        .fillna("")
        .str.contains(OUTER_MEMBRANE_TERM, regex=False)
    ]
    outer_membrane_average = outer_membrane_proteins[f"{residue_name} (%)"].mean()
    return whole_dataset_average, outer_membrane_average, len(outer_membrane_proteins)


def average_methionine_content(proteins):
    """Return average methionine percentage for all and outer-membrane proteins."""
    return average_residue_content(proteins, "Methionine")


def average_cysteine_content(proteins):
    """Return average cysteine percentage for all and outer-membrane proteins."""
    return average_residue_content(proteins, "Cysteine")


if __name__ == "__main__":
    filtered_proteins = filter_by_mass_and_methionine()
    print(f"Wrote {len(filtered_proteins)} proteins to {OUTPUT_PATH}")
    all_proteins = add_methionine_content(pd.read_excel(INPUT_PATH))
    (
        whole_dataset_average,
        outer_membrane_average,
        outer_membrane_count,
    ) = average_methionine_content(all_proteins)
    print(f"Whole dataset average methionine content: {whole_dataset_average:.2f}%")
    print(
        f"Outer-membrane protein average methionine content "
        f"({outer_membrane_count} proteins): {outer_membrane_average:.2f}%"
    )
    filtered_cysteines = filter_by_mass_and_cysteine()
    print(f"\nWrote {len(filtered_cysteines)} proteins to {CYSTEINE_OUTPUT_PATH}")
    all_proteins = add_cysteine_content(pd.read_excel(INPUT_PATH))
    (
        whole_dataset_average,
        outer_membrane_average,
        outer_membrane_count,
    ) = average_cysteine_content(all_proteins)
    print(f"Whole dataset average cysteine content: {whole_dataset_average:.2f}%")
    print(
        f"Outer-membrane protein average cysteine content "
        f"({outer_membrane_count} proteins): {outer_membrane_average:.2f}%"
    )