import pandas as pd
import os
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams

matplotlib.use('WxAgg')

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 0.75
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 14
rcParams['font.sans-serif'] = "Arial"


def make_graph(df, m_name="G", s_name=None):
    if s_name is None:
        s_name=m_name+"E"
    if m_name == "S":
        r = 4
    else:
        r=2
    print(df)
    els = {}
    edges = [("WT", "R224A"), ("WT", "W14A"), ("WT", "DM2"), ("R224A", "DM2"), ("W14A", "DM2"), ("R224A", "W14A")]
    G = nx.Graph()
    weights = []
    for i in range(0, 6):
        e1 = edges[i][0]
        e2 = edges[i][1]
        mutant = e1 + "_" + e2
        row = df[df["Mutant"] == mutant]
        m = float(row[m_name])
        s = np.abs(float(row[s_name]))

        G.add_edge(e1, e2, weight=m)
        label = str(round(m, r)) + "\n± " + str(round(s, r))
        els.update({edges[i]: label})
        sm = np.sqrt(3)
        t = 3.182
        cil = m - s #* t/sm
        clh = m + s #* t/sm
        if cil < 0 < clh:
            weights.append(0)
        else:
            weights.append(m)
    return G, els, weights


path = "Z:\\Group Share\\FAM\\FAM_HSJ_AqpZ CL\\Final8\\"
os.chdir(path)
df = pd.read_excel("CombinedResults2.xlsx")
names=[ "H", "TS", "GA"]#["G", "GA", "S", "H", "TS"]
for name in names:
    #name = "G"
    outstring = "Square_"+name+"_"
    if name == "S":
        minmax = [-0.01, 0.01]
    elif name == "H" or name == "TS":
        minmax = [-3, 3]
    else:
        minmax=[-0.75, 0.75]
    if name in ["H", "S", "GA", "TS"]:
        temps = [25]
    else:
        temps = [15, 20, 25, 30, 35]
    for j in range(0, len(temps)):
        fig = plt.figure(figsize=(6, 8))
        t = temps[j]
        for i in range(0, 6):
            df2 = df[df["# Lnum"] == i]
            df2 = df2[df2["Temp"] == t]

            G, els, weights = make_graph(df2, m_name=name)

            p = {"WT": [0, 1], "DM2": [1, 0], "R224A": [1, 1], "W14A": [0, 0]}
            ax1 = plt.subplot(3, 2, i + 1, aspect="equal")
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            nx.draw(G, p)
            ns = 1800
            fs = 12
            bbox = {"boxstyle": "round, pad=0.1", "ec": "None", "fc": "None", "alpha": 0.5}
            # weights = list(nx.get_edge_attributes(G, 'weight').values())
            nx.draw_networkx_nodes(G, p, margins=0.15, node_color="w", node_shape="o", node_size=ns, edgecolors="k")
            nx.draw_networkx_edges(G, p, width=32, edge_color="k")
            nx.draw_networkx_edges(G, p, width=30, edge_color=weights, edge_cmap=plt.cm.bwr,
                                   node_size=ns, edge_vmin=minmax[0], edge_vmax=minmax[1])
            nx.draw_networkx_labels(G, p, font_size=fs)
            nx.draw_networkx_edge_labels(G, p, edge_labels=els, font_size=fs, bbox=bbox, font_color="k")

            plt.title("Lipid " + str(i) + "→" + str(i + 1))
        plt.tight_layout()
        figfile = outstring + str(t)
        fig.suptitle('ΔΔ' + name + " "+ str(t) + ' °C', fontsize=16)
        plt.savefig(figfile + ".png")
        plt.savefig(figfile + ".pdf")

plt.show()
