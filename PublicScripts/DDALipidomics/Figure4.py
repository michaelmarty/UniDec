"""Plot CL (panel A, left) and TG (panel B, right), with five dataset networks each.

Reuse SI 5-6's scoring, eligibility rules, and per-lipid normalized colors.
Build networks from accepted, non-provisional spectra only, before calculating
strength. Defaults: CL [M-H]- and TG [M+Na]+ across all
five datasets. Save one PNG (600 dpi) and PDF directly beside the input in the DDA Pipeline folder.
The 6.5 x 6.75-inch figure fits a page with Arial text at its final print size.
Use --input, --output-dir, --datasets, or --no-show as needed.
"""

import argparse
from pathlib import Path
import sys

import matplotlib as mpl
if __name__ == "__main__" and "--no-show" in sys.argv:
    mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Make the package available when running this script directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import networkx as nx
from PublicScripts.Lipids.LipidFunctions import network_plot
from Figure3 import datasets as DATASETS
from SI_8_9 import (INPUT_FILE, MIN_LIPIDS, NORMALIZED, NETWORK_GROUPS,
                    recompute_adduct_networks, prepare_strength, comparable_points)

CONFIRMED_ONLY = True
SHOW_NODE_LABELS = False
DATASET_DISPLAY_NAMES = {'OTOT': 'OT/OT', 'OTIT': 'OT/IT', 'ITIT': 'IT/IT', 'Stellarboth': 'Stellar±'}


PANELS = [("CL", "[M-H]-"), ("EtherPC", "[M+H]+")]
OUTPUT_BASENAME = "Figure4"
FIGURE_SIZE = (6.5, 6.75)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def prepare_network_context(source, confirmed_only=CONFIRMED_ONLY):
    """Use the distribution pipeline unchanged; retain its feature graphs for display."""
    graphs = {}
    scored, _ = recompute_adduct_networks(source, graph_store=graphs, accepted_only=confirmed_only)
    points = prepare_strength(scored, all_ids=not confirmed_only, allow_empty=True)
    scored["Accepted"] = scored["Overall Match Status"].isin([True, "TRUE", "True", "true", 1, "1"])
    scored["Known Status"] = scored["Accepted"] | scored["Overall Match Status"].isin(
        [False, "FALSE", "False", "false", 0, "0"])
    return scored, points, graphs


def build_class_networks(frame, lipid_class, datasets, confirmed_only=CONFIRMED_ONLY,
                          adduct=None, context=None):
    scored, points, source_graphs = (prepare_network_context(frame, confirmed_only)
                                     if context is None else context)
    selected = points.loc[points["Class Group"].str.casefold().eq(lipid_class.strip().casefold())].copy()
    if adduct:
        selected = selected.loc[selected["Adduct type"].eq(adduct.strip())]
    # SI 6-7 excludes original-class/adduct groups without positive strengths.
    selected = selected.loc[selected["Normalization maximum"].gt(0)]
    networks = []
    canonical_datasets = {name.casefold(): name for name in scored["Dataset"].dropna().unique()}
    datasets = [canonical_datasets.get(name.strip().casefold(), name.strip()) for name in datasets]
    selected, _ = comparable_points(selected, MIN_LIPIDS, dataset_order=datasets)
    for polarity, mode_points in selected.groupby("Polarity", sort=False):
        counts = mode_points.groupby("Dataset")["Lipid ID"].nunique()
        graphs = []
        for dataset in datasets:
            dataset_points = mode_points.loc[mode_points["Dataset"].eq(dataset)]
            graph = nx.Graph()
            for key, group in dataset_points.groupby(NETWORK_GROUPS, sort=False):
                original = source_graphs[key]
                by_identity = group.set_index("Lipid ID")
                included = []
                attributes = {}
                labels = {}
                for node, attrs in original.nodes(data=True):
                    row = scored.loc[attrs["Input Row"]]
                    accepted = row["Accepted"]
                    provisional = str(row.get("Validation Status", "")).strip().casefold() == "provisional"
                    identity = row["Ontology"].casefold() + "|" + row["Metabolite name"].split("|", 1)[0].strip().casefold()
                    if provisional or not row["Known Status"] or (confirmed_only and not accepted):
                        continue
                    if identity not in by_identity.index:
                        continue
                    included.append(node)
                    point = by_identity.loc[identity]
                    attributes[node] = {NORMALIZED: float(point[NORMALIZED]),
                                        "Lipid Median Class Strength": float(point["Class Strength"]),
                                        "Normalization maximum": float(point["Normalization maximum"]),
                                        "Lipid ID": identity}
                    labels[node] = str(attrs["Input Row"]) + "~" + row["Metabolite name"]
                display = original.subgraph(included).copy()
                nx.set_node_attributes(display, attributes)
                display = nx.relabel_nodes(display, labels)
                graph = nx.compose(graph, display)
            graph.graph["Title"] = (f"{DATASET_DISPLAY_NAMES.get(dataset, dataset)}\n"
                                    f"{counts[dataset]} lipids ({graph.number_of_nodes()} spectra), "
                                    f"{graph.number_of_edges():,} links")
            graph.graph["Dataset"] = dataset
            graph.graph["Eligible Lipids"] = int(counts[dataset])
            graphs.append(graph)
        networks.append((polarity, graphs))
    return networks


def plot_cl_tg(frame, datasets=DATASETS, output_dir=None):
    """Draw CL/TG panels using accepted-only feature graphs and SI 5-6 colors."""
    if not datasets:
        raise ValueError("Select at least one dataset.")
    context = prepare_network_context(frame, confirmed_only=True)
    prepared = []
    for lipid_class, adduct in PANELS:
        groups = build_class_networks(frame, lipid_class, datasets,
                                               confirmed_only=True, adduct=adduct, context=context)
        if len(groups) != 1:
            raise ValueError(f"{lipid_class} {adduct} does not meet the SI 5-6 eligibility rules.")
        prepared.append(groups[0])

    fig = plt.figure(figsize=FIGURE_SIZE)
    outer = fig.add_gridspec(1, 2, left=0.035, right=0.98, bottom=0.16,
                             top=0.99, wspace=0.14)
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    for panel, ((lipid_class, adduct), (polarity, graphs)) in enumerate(zip(PANELS, prepared)):
        inner = outer[panel].subgridspec(2, 1, height_ratios=[0.08, 1], hspace=0.22)
        header = fig.add_subplot(inner[0, 0])
        columns = min(2, len(datasets))
        rows = (len(datasets) + columns - 1) // columns
        network_grid = inner[1, 0].subgridspec(rows, columns * 2, hspace=0.65, wspace=0.25)
        header.set_axis_off()
        header.text(0, 0.5, chr(65 + panel), fontsize=12, fontweight="bold", va="center")
        header.text(0.55, 0.5, f"{lipid_class} | {adduct}",
                    fontsize=10, fontweight="bold", ha="center", va="center")
        for index, graph in enumerate(graphs):
            column = 1 if columns == 2 and index == len(graphs) - 1 and len(graphs) % 2 else 2 * (index % columns)
            ax = fig.add_subplot(network_grid[index // columns, column:column + 2])
            network_plot(graph, axis=ax, colortarget=NORMALIZED)
            dataset = graph.graph["Dataset"]
            ax.set_title(DATASET_DISPLAY_NAMES.get(dataset, dataset),
                         fontsize=9, fontweight="bold", pad=23)
            ax.text(0.5, 1.02,
                    f"{graph.graph['Eligible Lipids']} lipids / {graph.number_of_nodes()} spectra\n"
                    f"{graph.number_of_edges():,} links",
                    fontsize=7, ha="center", va="bottom", transform=ax.transAxes)
            values = np.array([attrs[NORMALIZED] for _, attrs in graph.nodes(data=True)])
            for collection in ax.collections:
                if isinstance(collection, mpl.collections.PathCollection):
                    collection.set_sizes([18])
                    collection.set_edgecolor("white")
                    collection.set_linewidth(0.25)
                    collection.set_cmap("viridis")
                    collection.set_norm(norm)
                    collection.set_array(values)
                elif isinstance(collection, mpl.collections.LineCollection):
                    collection.set_alpha(0.40)
                    collection.set_linewidths(np.array(collection.get_linewidths()) * 0.8)
            for label in ax.texts:
                if label.get_transform() != ax.transAxes:
                    label.set_visible(SHOW_NODE_LABELS)
                    label.set_fontsize(5)
            ax.set_axis_off()

    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="viridis"),
                            cax=fig.add_axes([0.15, 0.075, 0.70, 0.022]),
                            orientation="horizontal")
    colorbar.set_label("Within-dataset normalized class strength", fontsize=9, labelpad=5)
    colorbar.ax.tick_params(labelsize=8, length=3)
    colorbar.outline.set_linewidth(0.5)
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        for extension in ["png", "pdf"]:
            path = output_dir / f"{OUTPUT_BASENAME}.{extension}"
            fig.savefig(path, dpi=600)
            print(f"Saved {path}")
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT_FILE)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--datasets", nargs="+", default=DATASETS)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)
    source = (pd.read_csv(args.input) if args.input.suffix.lower() == ".csv"
              else pd.read_excel(args.input, sheet_name="Sheet1"))
    figure = plot_cl_tg(source, args.datasets, args.output_dir or args.input.parent / "Final Figures")
    if args.no_show:
        plt.close(figure)
    else:
        plt.show()


if __name__ == "__main__":
    main()
