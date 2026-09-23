from types import SimpleNamespace

from matplotlib.figure import Figure
import isogen
import pandas as pd

from isodec import match_fragments
from isodec.fragment_view import plot_fragment_matches


def test_fragment_view_marks_wrapped_cleavages_and_both_termini():
    sequence = "PEPTIDE"
    masses = isogen.calc_pep_fragments(sequence, ion_types="by")
    pks = SimpleNamespace(peaks=[SimpleNamespace(monoiso=masses[label])
                                 for label in ("b3", "y3")])
    match_fragments(pks, sequence, ion_types="by")
    ax = Figure().subplots()

    plot_fragment_matches(ax, sequence, pks, residues_per_line=4)

    assert [text.get_text() for text in ax.texts if text.get_text() in ("1", "5")] == ["1", "5"]
    assert len(ax.lines) == 4  # Each matched cleavage has a horizontal and vertical mark.
    assert {line.get_color() for line in ax.lines} == {"#dc2626", "#2563eb"}
    assert "33.3%" in ax.get_title(loc="left")  # Two matched sites among six possible cleavages.


def test_fragment_view_colors_modified_residue_red():
    sequence = "S[Acetylation]HHS"
    pks = SimpleNamespace(peaks=[])
    match_fragments(pks, sequence)
    ax = Figure().subplots()

    plot_fragment_matches(ax, sequence, pks, residues_per_line=2)

    letters = [text for text in ax.texts if text.get_text() in "SH"]
    assert [(text.get_text(), text.get_color()) for text in letters] == [
        ("S", "#dc2626"), ("H", "black"), ("H", "black"), ("S", "black")
    ]


def test_wrapped_z_mark_moves_left_while_c_mark_stays_right():
    table = pd.DataFrame(index=range(1, 8), columns=["c_match", "z'_match"])
    table.loc[4] = [1.0, 1.0]
    pks = SimpleNamespace(fragment_matches=table, sequence_coverage=1 / 7,
                          fragment_match_percent=100)
    ax = Figure().subplots()

    plot_fragment_matches(ax, "PEPTIDEK", pks, residues_per_line=4)

    c_line, _, z_line, _ = ax.lines
    assert tuple(c_line.get_xdata()) == (3.3, 4)
    assert tuple(z_line.get_xdata()) == (0, 0.7)
    assert z_line.get_ydata()[0] > c_line.get_ydata()[0]
