"""Sequence coverage view for :func:`isodec.match_fragments` results."""

import math

import isogen
from isogen.protein_mods import _parse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ION_COLORS = {
    "a": "#d97706", "b": "#dc2626", "c": "#7c3aed",
    "x": "#0891b2", "y": "#2563eb", "z": "#16a34a",
}


def plot_fragment_matches(ax, sequence, pks, residues_per_line=30):
    """Draw matched cleavage marks on an existing Matplotlib axis.

    ``pks`` must have been processed by ``match_fragments`` for ``sequence``.
    N-terminal ions appear above each sequence row; C-terminal ions below it.
    The row's left number is the one-based index of its first residue.
    """
    if residues_per_line < 1:
        raise ValueError("residues_per_line must be positive")
    residues, modifications, _, _ = _parse(sequence)
    modified = set()
    for _, _, positions in modifications:
        if isinstance(positions, frozenset):
            modified.update(positions)
        elif positions == "N-term":
            modified.add(0)
        elif positions == "C-term":
            modified.add(len(residues) - 1)
    table = pks.fragment_matches
    lines = max(1, math.ceil(len(residues) / residues_per_line))
    ax.clear()

    ion_types = [column[:-6] for column in table if column.endswith("_match")
                 and table[column].notna().any()]
    upper = [ion for ion in ion_types if ion[0] in "abc"]
    lower = [ion for ion in ion_types if ion[0] in "xyz"]
    # Keep related series together while giving modified ions their own lane.
    upper.sort(key=lambda ion: ("abc".index(ion[0]), ion))
    lower.sort(key=lambda ion: ("xyz".index(ion[0]), ion))
    row_height = 1.6 + 0.23 * (len(upper) + len(lower))

    for line in range(lines):
        start = line * residues_per_line
        chunk = residues[start:start + residues_per_line]
        baseline = line * row_height + 0.65 + 0.23 * len(upper)
        ax.text(-1.5, baseline, str(start + 1), ha="right", va="center",
                color="#64748b", fontsize=9)
        for position, residue in enumerate(chunk):
            ax.text(position + 0.5, baseline, residue, ha="center", va="center",
                    fontfamily="monospace", fontsize=12, fontweight="bold",
                    color="#dc2626" if start + position in modified else "black")

        for position in range(1, len(chunk) + 1):
            cleavage = start + position
            if cleavage not in table.index:
                continue
            for series, side in ((upper, -1), (lower, 1)):
                for lane, ion in enumerate(series):
                    if math.isnan(table.at[cleavage, f"{ion}_match"]):
                        continue
                    y = baseline + side * (0.45 + lane * 0.23)
                    color = ION_COLORS[ion[0]]
                    # The vertical tip points to the actual peptide bond.
                    x0, x1 = ((position - 0.8, position) if side < 0
                              else (position, position + 0.8))
                    ax.plot((x0, x1), (y, y), color=color, lw=1.8,
                            solid_capstyle="round")
                    ax.plot((position, position), (y, y - side * 0.17),
                            color=color, lw=1.8)

    ax.set_xlim(-2.5, residues_per_line + 1)
    ax.set_ylim(lines * row_height, -0.3)
    ax.set_axis_off()
    ax.set_title("Fragment matches  |  sequence coverage {:.1%}  |  matched peaks {:.1f}%".format(
        pks.sequence_coverage, pks.fragment_match_percent), loc="left", fontsize=12)
    if ion_types:
        ax.legend([Line2D([0], [0], color=ION_COLORS[ion[0]], lw=2)
                   for ion in ion_types], ion_types, ncol=len(ion_types),
                  loc="lower center", bbox_to_anchor=(0.5, -0.03),
                  frameon=False, fontsize=9)
    return ax


def show_fragment_matches(sequence, pks, residues_per_line=30):
    """Open a standalone Matplotlib window with the sequence coverage view."""
    residues = isogen.strip_proforma(sequence)
    if residues_per_line < 1:
        raise ValueError("residues_per_line must be positive")
    lines = max(1, math.ceil(len(residues) / residues_per_line))
    fig, ax = plt.subplots(figsize=(13, max(4, lines * 1.1 + 1.5)))
    plot_fragment_matches(ax, sequence, pks, residues_per_line)
    fig.tight_layout()
    plt.show()
    return fig
