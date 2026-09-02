import os
import pandas as pd
import numpy as np

def get_unique_enzymes(df, column="UniProtKB AC(s)"):
    # Get unique entries. Will need to collapse -1 and -2 isoforms into the same
    unique_entries = df[column].dropna().unique()
    unique_enzymes = set()
    for entry in unique_entries:
        # Split by hyphen and strip whitespace
        proteins = [protein.strip() for protein in entry.split("-")]
        for protein in proteins:
            if protein[0] == "P" or protein[0] == "Q":  # Only consider entries starting with P or Q
                unique_enzymes.add(protein)
    # Sort the unique enzymes for better readability
    unique_enzymes = sorted(unique_enzymes)
    return unique_enzymes


def filter_by_reaction_text(df, reaction_text, rtcol="Reaction text"):
    # Filter the dataframe based on the reaction text
    filtered_df = df[df[rtcol].str.contains(reaction_text, case=False, na=False)]
    return filtered_df


filepath = r"Y:\Data\SwissLipids\human_enzymes.xlsx"

df = pd.read_excel(filepath)

df = filter_by_reaction_text(df, "choline")

unique_enzymes = get_unique_enzymes(df)

print(f"Found {len(unique_enzymes)} unique enzymes.")