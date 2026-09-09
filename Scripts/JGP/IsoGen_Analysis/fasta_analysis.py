from pathlib import Path
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FASTA_DIR = Path(r"C:\Data\IsoNN\fastas")
os.chdir(FASTA_DIR)
PROTEOMES = {
    "Yeast": "yeast.fasta",
    "Human": "human.fasta",
    "Mouse": "mouse.fasta",
    "E. coli": "ecoli.fasta",
}
OUTPUT_PATH = FASTA_DIR / "proteome_cysteine_methionine_vs_mass.png"
NORMALIZED_OUTPUT_PATH = (
    FASTA_DIR / "proteome_cysteine_methionine_fraction_vs_mass.png"
)

# Average residue masses (Da). Adding H2O converts the sum of residue masses to
# the neutral mass of the intact protein.
RESIDUE_MASSES = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}
WATER_MASS = 18.0153


def read_fasta(path):
    """Return complete protein sequences from a possibly multi-line FASTA."""
    sequences = []
    sequence_parts = []

    with path.open(encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sequence_parts:
                    sequences.append("".join(sequence_parts).upper())
                    sequence_parts = []
            else:
                sequence_parts.append(line)

    if sequence_parts:
        sequences.append("".join(sequence_parts).upper())
    print(len(sequences), "sequences read from", path)
    return sequences


def protein_mass(sequence):
    """Calculate average neutral protein mass; return NaN for ambiguous residues."""
    try:
        return WATER_MASS + sum(RESIDUE_MASSES[residue] for residue in sequence)
    except KeyError:
        return np.nan


def load_proteomes():
    records = []
    for proteome, filename in PROTEOMES.items():
        for sequence in read_fasta(FASTA_DIR / filename):
            sequence_length = len(sequence)
            cysteine_count = sequence.count("C")
            methionine_count = sequence.count("M")
            records.append(
                {
                    "Proteome": proteome,
                    "Mass (kDa)": protein_mass(sequence) / 1000,
                    "Cysteine": cysteine_count,
                    "Methionine": methionine_count,
                    "Cysteine (%)": 100 * cysteine_count / sequence_length,
                    "Methionine (%)": 100 * methionine_count / sequence_length,
                }
            )
    return pd.DataFrame(records).dropna(subset=["Mass (kDa)"])


def bin_by_mass(data, bin_width=2, max_mass=100, min_mass=10):
    """Calculate residue-count means in fixed-width protein-mass bins."""
    bins = np.arange(min_mass, max_mass + bin_width, bin_width)
    binned_data = data.loc[data["Mass (kDa)"] <= max_mass].copy()
    binned_data["Mass bin"] = pd.cut(
        binned_data["Mass (kDa)"], bins=bins, include_lowest=True
    )
    averages = (
        binned_data.groupby(["Proteome", "Mass bin"], observed=True)[
            ["Cysteine", "Methionine", "Cysteine (%)", "Methionine (%)"]
        ]
        .mean()
        .reset_index()
    )
    averages["Mass (kDa)"] = averages["Mass bin"].map(
        lambda interval: interval.mid
    )
    return averages


def plot_residues_vs_mass(averages):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharex=True)

    for axis, residue in zip(axes, ["Cysteine", "Methionine"]):
        for proteome in PROTEOMES:
            proteome_data = averages[averages["Proteome"] == proteome]
            axis.plot(
                proteome_data["Mass (kDa)"],
                proteome_data[residue],
                linewidth=2,
                label=proteome,
            )
        all_proteome_average = (
            averages.pivot(
                index="Mass (kDa)", columns="Proteome", values=residue
            )
            .dropna(subset=list(PROTEOMES))
            .mean(axis=1)
        )
        axis.plot(
            all_proteome_average.index,
            all_proteome_average,
            color="black",
            linewidth=2.5,
            linestyle="--",
            label="4-proteome average",
        )
        axis.set_title(f"{residue} (C)" if residue == "Cysteine" else "Methionine (M)")
        axis.set_xlabel("Protein mass (kDa)")
        axis.set_ylabel(f"Average number of {residue.lower()} residues")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_PATH.with_suffix(".pdf"), bbox_inches="tight")


def plot_normalized_residues_vs_mass(averages):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharex=True)

    for axis, residue in zip(axes, ["Cysteine", "Methionine"]):
        normalized_column = f"{residue} (%)"
        for proteome in PROTEOMES:
            proteome_data = averages[averages["Proteome"] == proteome]
            axis.plot(
                proteome_data["Mass (kDa)"],
                proteome_data[normalized_column],
                linewidth=2,
                label=proteome,
            )
        all_proteome_average = (
            averages.pivot(
                index="Mass (kDa)",
                columns="Proteome",
                values=normalized_column,
            )
            .dropna(subset=list(PROTEOMES))
            .mean(axis=1)
        )
        axis.plot(
            all_proteome_average.index,
            all_proteome_average,
            color="black",
            linewidth=2.5,
            linestyle="--",
            label="4-proteome average",
        )
        axis.set_title(f"{residue} (C)" if residue == "Cysteine" else "Methionine (M)")
        axis.set_xlabel("Protein mass (kDa)")
        axis.set_ylabel(f"Average {residue.lower()} content (% of residues)")
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(NORMALIZED_OUTPUT_PATH, dpi=300, bbox_inches="tight")
    fig.savefig(NORMALIZED_OUTPUT_PATH.with_suffix(".pdf"), bbox_inches="tight")


if __name__ == "__main__":
    proteome_data = load_proteomes()
    mass_binned_averages = bin_by_mass(proteome_data)
    plot_residues_vs_mass(mass_binned_averages)
    plot_normalized_residues_vs_mass(mass_binned_averages)
    plt.show()
