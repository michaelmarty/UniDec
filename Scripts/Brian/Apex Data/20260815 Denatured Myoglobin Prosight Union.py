#!/usr/bin/env python3
"""
Union sequence coverage across ProSight Lite exports.

Collapses matched fragment ions to backbone cleavage SITES, so b/y ions (HCD,
IRMPD) and c/z ions (ETD, EThcD) combine correctly. Reports per-condition
coverage, the union, and pairwise orthogonality. Also writes ProSight-style
fusion fragment maps as publication-ready PNG and SVG files.

Optional: define DOMAINS to get per-domain coverage (e.g. MSP vs mCherry barrel).

Usage:
    python prosight_union.py                 # uses configured files below
    python prosight_union.py *.csv           # or override them on the command line
"""

import sys
import glob
import re
from pathlib import Path
from itertools import combinations

import pandas as pd

# ----------------------------------------------------------------------------
# CONFIG -- edit these
# ----------------------------------------------------------------------------

SEQUENCE = (
    # Paste the EXACT sequence you used in ProSight Lite, including tag/TEV
    # and with the same residue numbering. One string, no spaces or newlines.
    "GLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELGFQG"
)

# Exact ProSight Lite fragment exports to combine when no CLI paths are given.
INPUT_FILES = [
    r"Z:\Group Share\BHT\Data\Apex Data\20260815\IRMPD Fragments.csv",
    r"Z:\Group Share\BHT\Data\Apex Data\20260815\HCD Fragments.csv",
    r"Z:\Group Share\BHT\Data\Apex Data\20260815\IR-ETD Fragments.csv",
    r"Z:\Group Share\BHT\Data\Apex Data\20260815\IR-EThcD Fragments.csv",
    r"Z:\Group Share\BHT\Data\Apex Data\20260815\HCD Fragments.csv",
    r"Z:\Group Share\BHT\Data\Apex Data\20260815\EThcD Fragments.csv",
]

# Fallback glob used only when INPUT_FILES is empty.
INPUT_GLOB = "*.csv"

# ProSight-style graphical output names. Results are written to the directory
# containing the input fragment files, independent of the working directory.
MAP_OUTPUT_NAMES = ("fusion_fragment_map.png", "fusion_fragment_map.svg")
RESIDUES_PER_ROW = 30
MAP_DPI = 300

# Optional domain boundaries, 1-based inclusive residue numbers in SEQUENCE.
# Set to {} to skip. Example:
# DOMAINS = {"tag": (1, 30), "MSP": (31, 250), "mCherry": (251, 480)}
DOMAINS = {}

# ----------------------------------------------------------------------------

NTERM_IONS = {"a", "b", "c"}
CTERM_IONS = {"x", "y", "z"}

# ProSight ion-family colors: a/x green, b/y blue, and c/z red.
ION_COLORS = {
    "a": "#16a34a",
    "b": "#1646ff",
    "c": "#ef1f2d",
    "x": "#16a34a",
    "y": "#1646ff",
    "z": "#ef1f2d",
}

ION_PATTERN = re.compile(r"^\s*([abcxyz])\s*[-_]?\s*(\d+)", re.IGNORECASE)


def find_col(df, candidates):
    """Case-insensitive fuzzy column lookup."""
    lowered = {c.lower().strip(): c for c in df.columns}
    for cand in candidates:
        if cand in lowered:
            return lowered[cand]
    for cand in candidates:
        for low, orig in lowered.items():
            if cand in low:
                return orig
    return None


def parse_ions(df):
    """Yield (ion_letter, ion_number) pairs from a ProSight-style table.

    Handles both a single combined label column (e.g. 'b12', 'c 12') and
    separate type/number columns.
    """
    label_col = find_col(df, ["fragment", "ion", "label", "name", "match"])
    type_col = find_col(df, ["ion type", "iontype", "type", "series"])
    num_col = find_col(df, ["ion number", "ionnumber", "number", "index",
                            "position", "residue"])

    if type_col and num_col:
        for t, n in zip(df[type_col], df[num_col]):
            t = str(t).strip().lower()
            if not t:
                continue
            letter = t[0]
            try:
                num = int(float(n))
            except (TypeError, ValueError):
                continue
            if letter in NTERM_IONS | CTERM_IONS:
                yield letter, num
        return

    if label_col is None:
        raise ValueError(
            "Could not identify fragment columns. Found: "
            + ", ".join(map(str, df.columns))
        )

    for val in df[label_col].astype(str):
        m = ION_PATTERN.match(val)
        if m:
            yield m.group(1).lower(), int(m.group(2))


def sites_from_file(path, seq_len):
    """Return all cleavage sites and the sites for each ion series."""
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep, engine="python")
    df.columns = [str(c).strip() for c in df.columns]

    sites = set()
    series_sites = {ion: set() for ion in "abcxyz"}
    dropped = 0
    for letter, num in parse_ions(df):
        if letter in NTERM_IONS:
            site = num
        else:
            site = seq_len - num
        if 1 <= site <= seq_len - 1:
            sites.add(site)
            series_sites[letter].add(site)
        else:
            dropped += 1

    if dropped:
        print(f"    note: {dropped} ion(s) fell outside 1..L-1 and were skipped")
    return sites, series_sites


def pct(sites, total_sites):
    return 100.0 * len(sites) / total_sites if total_sites else 0.0


def domain_coverage(sites, start, end):
    """Coverage within a domain, counting internal bonds only."""
    bonds = range(start, end)          # sites start..end-1
    total = len(list(bonds))
    hit = sum(1 for s in bonds if s in sites)
    return hit, total, (100.0 * hit / total if total else 0.0)


def write_fragment_maps(seq, series_sites, condition_names, output_directory):
    """Write ProSight-style union maps with ion-series colors and glyphs."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams["font.family"] = "Arial"

    residues_per_row = RESIDUES_PER_ROW
    if residues_per_row < 1:
        raise ValueError("RESIDUES_PER_ROW must be at least 1")

    row_count = (len(seq) + residues_per_row - 1) // residues_per_row
    fig_height = max(3.0, 1.1 + 0.50 * row_count)
    fig, ax = plt.subplots(figsize=(17, fig_height))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    sequence_font = max(10.0, min(18.0, 540.0 / residues_per_row))
    number_color = "#a8a8a8"

    for row_index in range(row_count):
        start = row_index * residues_per_row
        end = min(start + residues_per_row, len(seq))
        y = row_count - row_index - 1

        left_label = "N" if start == 0 else str(start + 1)
        right_label = "C" if end == len(seq) else str(end)
        ax.text(-1.15, y, left_label, ha="center", va="center",
                color=number_color, fontsize=sequence_font - 2)
        ax.text(residues_per_row + 0.65, y, right_label,
                ha="center", va="center", color=number_color,
                fontsize=sequence_font - 2)

        for column, residue in enumerate(seq[start:end]):
            ax.text(column, y, residue, ha="center", va="center",
                    color="black", family="Arial",
                    weight="bold", fontsize=sequence_font)

    def draw_glyph(site, ion, upper):
        """Draw a terminally oriented cleavage flag at one sequence bond."""
        row_index = (site - 1) // residues_per_row
        y = row_count - row_index - 1
        x = (site - 1) % residues_per_row + 0.5
        color = ION_COLORS[ion]
        line_width = 2.6

        # ProSight uses a diagonal shoulder for b/y glyphs and a square elbow
        # for c/z glyphs. Their elbow coordinates align at +/-0.29, while the
        # blue vertical stem stops sooner at the opposite end. Each glyph is a
        # continuous path so its corner cannot acquire an anti-aliasing seam.
        if ion in {"b", "y"}:
            if upper:
                xs = [x - 0.19, x, x]
                ys = [y + 0.36, y + 0.29, y - 0.19]
            else:
                xs = [x, x, x + 0.19]
                ys = [y + 0.19, y - 0.29, y - 0.36]
        elif upper:
            xs = [x - 0.19, x, x]
            ys = [y + 0.29, y + 0.29, y - 0.29]
        else:
            xs = [x, x, x + 0.19]
            ys = [y + 0.29, y - 0.29, y - 0.29]

        ax.plot(xs, ys, color=color, linewidth=line_width,
                solid_capstyle="butt", solid_joinstyle="miter", zorder=3)

    # Overlapping ion types share the exact cleavage coordinate. Draw square
    # red c/z glyphs first and the slightly shorter blue b/y glyphs afterward,
    # allowing both shapes to remain visible without horizontal staggering.
    for site in range(1, len(seq)):
        ions = [ion for ion in "czaxby" if site in series_sites[ion]]
        for ion in ions:
            draw_glyph(site, ion, upper=(ion in NTERM_IONS))

    union = set().union(*series_sites.values())
    total_sites = len(seq) - 1
    coverage = pct(union, total_sites)
    methods = ", ".join(
        name.removesuffix(" Fragments") for name in condition_names
    )
    fig.suptitle("Fusion Fragmentation Map", x=0.035, y=0.975,
                 ha="left", fontsize=15, weight="bold")
    fig.text(0.035, 0.910,
             f"{len(union)}/{total_sites} cleavage sites ({coverage:.1f}%)"
             f"  |  {methods}",
             ha="left", va="top", fontsize=9.5, color="#555555")
    fig.text(0.965, 0.910,
             "upper: N-terminal   lower: C-terminal",
             ha="right", va="top", fontsize=9.5, color="#777777")

    legend_handles = [
        Line2D([0], [0], color=ION_COLORS[ion], linewidth=3,
               label=label)
        for ion, label in (("a", "a/x"), ("b", "b/y"), ("c", "c/z"))
    ]
    fig.legend(handles=legend_handles, loc="upper right",
               bbox_to_anchor=(0.965, 0.985), ncol=3, frameon=False,
               handlelength=1.6, columnspacing=1.2, fontsize=9.5)

    ax.set_xlim(-1.9, residues_per_row + 1.4)
    ax.set_ylim(-0.45, row_count - 0.60)
    ax.axis("off")
    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.87)

    outputs = []
    for output_name in MAP_OUTPUT_NAMES:
        output = output_directory / output_name
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=MAP_DPI, bbox_inches="tight",
                    pad_inches=0.04, facecolor="white")
        outputs.append(output.resolve())
    plt.close(fig)
    return outputs


def main():
    seq = "".join(SEQUENCE.split()).upper()
    if not seq:
        sys.exit("Set SEQUENCE at the top of the script first.")

    L = len(seq)
    total_sites = L - 1
    print(f"Sequence length: {L} residues  ({total_sites} cleavage sites)\n")

    configured_files = INPUT_FILES or sorted(glob.glob(INPUT_GLOB))
    files = [Path(p) for p in (sys.argv[1:] or configured_files)]
    if not files:
        sys.exit(f"No input files matched {INPUT_GLOB!r}")
    output_directory = files[0].resolve().parent

    per_condition = {}
    union_series = {ion: set() for ion in "abcxyz"}
    for path in files:
        print(f"Reading {path.name}")
        try:
            sites, series_sites = sites_from_file(path, L)
        except Exception as exc:
            print(f"    SKIPPED: {exc}")
            continue
        per_condition[path.stem] = sites
        for ion, ion_sites in series_sites.items():
            union_series[ion].update(ion_sites)
        print(f"    {len(sites)} sites -> {pct(sites, total_sites):.1f}% coverage")

    if not per_condition:
        sys.exit("Nothing parsed successfully.")

    union = set().union(*per_condition.values())

    print("\n" + "=" * 60)
    print("PER-CONDITION COVERAGE")
    print("=" * 60)
    for name, sites in sorted(per_condition.items(),
                              key=lambda kv: -len(kv[1])):
        print(f"  {name:<32} {len(sites):>5} sites  {pct(sites, total_sites):>6.1f}%")

    print("-" * 60)
    print(f"  {'UNION':<32} {len(union):>5} sites  {pct(union, total_sites):>6.1f}%")

    best = max(per_condition.values(), key=len)
    print(f"\n  Gain over best single condition: "
          f"+{pct(union, total_sites) - pct(best, total_sites):.1f} percentage points")

    # Unique contribution: sites found by exactly one condition
    print("\n" + "=" * 60)
    print("UNIQUE CONTRIBUTION (sites no other condition found)")
    print("=" * 60)
    for name, sites in sorted(per_condition.items()):
        others = set().union(
            *[s for n, s in per_condition.items() if n != name]
        ) if len(per_condition) > 1 else set()
        uniq = sites - others
        print(f"  {name:<32} {len(uniq):>5} sites  {pct(uniq, total_sites):>6.1f}%")

    # Pairwise overlap -- low Jaccard = orthogonal chemistry
    if len(per_condition) > 1:
        print("\n" + "=" * 60)
        print("PAIRWISE OVERLAP (Jaccard; lower = more orthogonal)")
        print("=" * 60)
        for a, b in combinations(sorted(per_condition), 2):
            sa, sb = per_condition[a], per_condition[b]
            j = len(sa & sb) / len(sa | sb) if (sa | sb) else 0.0
            print(f"  {a} vs {b}: {j:.2f}  "
                  f"({len(sa & sb)} shared, {len(sa | sb)} combined)")

    # Per-domain breakdown
    if DOMAINS:
        print("\n" + "=" * 60)
        print("PER-DOMAIN COVERAGE")
        print("=" * 60)
        header = f"  {'condition':<28}" + "".join(
            f"{d:>14}" for d in DOMAINS
        )
        print(header)
        rows = list(sorted(per_condition.items())) + [("UNION", union)]
        for name, sites in rows:
            line = f"  {name:<28}"
            for dname, (start, end) in DOMAINS.items():
                hit, tot, p = domain_coverage(sites, start, end)
                line += f"{p:>13.1f}%"
            print(line)
        print("\n  (domain bonds counted internally; boundary bonds excluded)")

    # Write the union site list for plotting / figure work
    out = output_directory / "union_sites.csv"
    rows = []
    for site in range(1, L):
        rows.append({
            "site": site,
            "residue_before": seq[site - 1],
            "residue_after": seq[site],
            "in_union": site in union,
            **{f"{ion}_ion": (site in union_series[ion]) for ion in "abcxyz"},
            **{name: (site in s) for name, s in sorted(per_condition.items())},
        })
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nWrote per-site matrix to {out.resolve()}")

    map_outputs = write_fragment_maps(
        seq, union_series, per_condition, output_directory
    )
    for map_output in map_outputs:
        print(f"Wrote fusion fragment map to {map_output}")


if __name__ == "__main__":
    main()
