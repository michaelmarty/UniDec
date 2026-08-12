"""Report human and mouse proteins near 95 kDa with low methionine content."""

import math

from fasta_analysis import FASTA_DIR, protein_mass
from fasta_analysis_outliers import protein_name_and_accession, read_fasta_records


PROTEOME_FASTAS = {"E. coli": FASTA_DIR / "ecoli.fasta"}
TARGET_MASS_KDA = 95
MASS_TOLERANCE_KDA = 5
TOP_N = 20


def find_low_methionine_proteins(fasta_path):
    min_mass = TARGET_MASS_KDA - MASS_TOLERANCE_KDA
    max_mass = TARGET_MASS_KDA + MASS_TOLERANCE_KDA
    proteins = []

    for header, sequence in read_fasta_records(fasta_path):
        mass_kda = protein_mass(sequence) / 1000
        if math.isnan(mass_kda) or not min_mass <= mass_kda <= max_mass:
            continue

        protein_name, accession = protein_name_and_accession(header)
        proteins.append(
            {
                "name": protein_name,
                "accession": accession,
                "mass_kda": mass_kda,
                "methionine_percent": 100 * sequence.count("M") / len(sequence),
            }
        )

    return sorted(
        proteins,
        key=lambda protein: (
            protein["methionine_percent"],
            abs(protein["mass_kda"] - TARGET_MASS_KDA),
            protein["name"],
        ),
    )


def average_zinc_finger_methionine(fasta_path):
    """Return mean methionine percentage for 90-100 kDa zinc finger proteins."""
    min_mass = TARGET_MASS_KDA - MASS_TOLERANCE_KDA
    max_mass = TARGET_MASS_KDA + MASS_TOLERANCE_KDA
    methionine_percentages = []

    for header, sequence in read_fasta_records(fasta_path):
        protein_name, _ = protein_name_and_accession(header)
        if "zinc finger protein" not in protein_name.lower():
            continue

        mass_kda = protein_mass(sequence) / 1000
        if math.isnan(mass_kda) or not min_mass <= mass_kda <= max_mass:
            continue

        methionine_percentages.append(
            100 * sequence.count("M") / len(sequence)
        )

    if not methionine_percentages:
        return float("nan"), 0
    return (
        sum(methionine_percentages) / len(methionine_percentages),
        len(methionine_percentages),
    )


def print_top_proteins(proteins, proteome):
    min_mass = TARGET_MASS_KDA - MASS_TOLERANCE_KDA
    max_mass = TARGET_MASS_KDA + MASS_TOLERANCE_KDA
    print(
        f"Top {TOP_N} methionine-poor {proteome.lower()} proteins from "
        f"{min_mass:g}-{max_mass:g} kDa"
    )
    print(f"{'Rank':>4}  {'Protein':<64} {'Mass (kDa)':>10} {'Methionine (%)':>15}")
    print("-" * 99)

    for rank, protein in enumerate(proteins[:TOP_N], start=1):
        protein_label = f"{protein['name']} [{protein['accession']}]"
        print(
            f"{rank:>4}  {protein_label:<64.64} "
            f"{protein['mass_kda']:>10.2f} "
            f"{protein['methionine_percent']:>15.2f}"
        )


if __name__ == "__main__":
    for index, (proteome, fasta_path) in enumerate(PROTEOME_FASTAS.items()):
        if index:
            print()
        print_top_proteins(
            find_low_methionine_proteins(fasta_path),
            proteome,
        )
        average_methionine, protein_count = average_zinc_finger_methionine(
            fasta_path
        )
        mass_range = (
            f"{TARGET_MASS_KDA - MASS_TOLERANCE_KDA:g}-"
            f"{TARGET_MASS_KDA + MASS_TOLERANCE_KDA:g} kDa"
        )
        if protein_count:
            print(
                f"\nAverage methionine content for all {protein_count} "
                f"{proteome} zinc finger proteins from {mass_range}: "
                f"{average_methionine:.2f}%"
            )
        else:
            print(
                f'\nNo {proteome} proteins containing "Zinc finger protein" '
                f"were found from {mass_range}."
            )
