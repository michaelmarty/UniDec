"""Report methionine-rich human proteins near 35 kDa."""

import math

from fasta_analysis import FASTA_DIR, protein_mass


PROTEOME_FASTAS = {
    "Human": FASTA_DIR / "human.fasta",
    "Mouse": FASTA_DIR / "mouse.fasta",
}
TARGET_MASS_KDA = 35
MASS_TOLERANCE_KDA = 5
TOP_N = 10


def read_fasta_records(path):
    """Yield (header, sequence) pairs from a multi-line FASTA file."""
    header = None
    sequence_parts = []

    with path.open(encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts).upper()
                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line)

    if header is not None:
        yield header, "".join(sequence_parts).upper()


def protein_name_and_accession(header):
    """Extract a readable protein name and accession from a UniProt header."""
    identifier, _, description = header.partition(" ")
    identifier_parts = identifier.split("|")
    accession = identifier_parts[1] if len(identifier_parts) > 1 else identifier
    protein_name = description.split(" OS=", maxsplit=1)[0] or identifier
    return protein_name, accession


def find_methionine_rich_proteins(fasta_path):
    min_mass = TARGET_MASS_KDA - MASS_TOLERANCE_KDA
    max_mass = TARGET_MASS_KDA + MASS_TOLERANCE_KDA
    proteins = []

    for header, sequence in read_fasta_records(fasta_path):
        mass_kda = protein_mass(sequence) / 1000
        if math.isnan(mass_kda) or not min_mass <= mass_kda <= max_mass:
            continue

        protein_name, accession = protein_name_and_accession(header)
        methionine_count = sequence.count("M")
        proteins.append(
            {
                "name": protein_name,
                "accession": accession,
                "mass_kda": mass_kda,
                "methionine_percent": 100 * methionine_count / len(sequence),
            }
        )

    return sorted(
        proteins,
        key=lambda protein: (
            -protein["methionine_percent"],
            abs(protein["mass_kda"] - TARGET_MASS_KDA),
            protein["name"],
        ),
    )


def average_olfactory_receptor_methionine(fasta_path):
    """Return the mean methionine percentage for named olfactory receptors."""
    methionine_percentages = []

    for header, sequence in read_fasta_records(fasta_path):
        protein_name, _ = protein_name_and_accession(header)
        if "olfactory receptor" in protein_name.lower():
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
        f"Top {TOP_N} methionine-rich {proteome.lower()} proteins from "
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
            find_methionine_rich_proteins(fasta_path),
            proteome,
        )
        average_methionine, receptor_count = (
            average_olfactory_receptor_methionine(fasta_path)
        )
        print(
            f"\nAverage methionine content for all {receptor_count} "
            f"{proteome.lower()} olfactory receptors: "
            f"{average_methionine:.2f}%"
        )
