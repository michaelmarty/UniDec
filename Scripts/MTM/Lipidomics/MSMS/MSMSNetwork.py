import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplcursors
import matchms as mms
import networkx as nx
from Scripts.MTM.Lipidomics.LipidFunctions import detect_outliers_mahalanobis

mpl.use("WxAgg")

def class_network_analysis(df):
    graphs = []
    outdf = pd.DataFrame()
    for c in df['Ontology'].unique():
        # Get subset for class
        subdf = df[df['Ontology'] == c].reset_index(drop=True)

        if len(subdf) < 3:
            outdf = pd.concat([outdf, subdf])
            continue

        G, newdf = network_analysis(subdf)
        graphs.append(G)
        outdf = pd.concat([outdf, newdf])
    outdf = outdf.reset_index(drop=True)
    return outdf, graphs

def robust_z(values_dict):
    # returns dict of z-scores
    vals = np.array(list(values_dict.values()), dtype=float)
    med = np.median(vals)
    mad = np.median(np.abs(vals - med))
    eps = 1e-9
    z = 0.6745 * (vals - med) / (mad + eps)
    return dict(zip(values_dict.keys(), z))

def create_msms_score_matrix(df, mz_tol=0.01):
    # Create a score matrix based on MS/MS similarity
    spectra = []
    nlspectra = []
    for i, row in df.iterrows():
        if pd.isna(row['MS/MS spectrum']):
            raise ValueError(f"Row {i} has no MS/MS spectrum data.")
        text = row['MS/MS spectrum']
        data = text.split(" ")
        data2 = np.array([d.split(":") for d in data]).astype(float)
        mz = row['Average Mz']
        metadata = {"precursor_mz": mz}
        nldata = data2.copy()
        nldata[:, 0] = mz - data2[:, 0]
        # Sort nldata by m/z
        nldata = nldata[np.argsort(nldata[:, 0])]
        # Take only nldata > 10
        nldata = nldata[nldata[:, 0] > 10]
        spec = mms.Spectrum(data2[:, 0], data2[:, 1], metadata=metadata)
        spectra.append(spec)
        nlspec = mms.Spectrum(nldata[:, 0], nldata[:, 1], metadata=metadata)
        nlspectra.append(nlspec)

    score_matrix = mms.calculate_scores(references=spectra, queries=spectra,
                                        similarity_function=mms.similarity.CosineGreedy(),
                                        is_symmetric=True, ).scores.to_array()
    score_matrix = np.array(score_matrix.tolist())[:, :, 0]

    nlscore_matrix = mms.calculate_scores(references=nlspectra, queries=nlspectra,
                                          similarity_function=mms.similarity.CosineGreedy(),
                                          is_symmetric=True).scores.to_array()
    nlscore_matrix = np.array(nlscore_matrix.tolist())[:, :, 0]

    score_matrix += nlscore_matrix
    score_matrix /= 2.0
    return score_matrix

def network_analysis(df, mz_tol=0.01):
    # Create score matrix
    score_matrix = create_msms_score_matrix(df, mz_tol=mz_tol)

    # Create graph from score matrix
    G = nx.Graph()
    G.add_weighted_edges_from([ (i, j, score_matrix[i,j]) for i in range(score_matrix.shape[0]) for j in range(i, score_matrix.shape[1]) if score_matrix[i,j] > 0.1 ])

    # Add node labels
    labels = {i: str(i) + "~"+ df.iloc[i]['Metabolite name'] for i in range(len(df))}
    nx.relabel_nodes(G, labels, copy=False)

    # Drop self-loops
    G.remove_edges_from(nx.selfloop_edges(G))

    # Add distance attribute for shortest-path calculations
    for u, v, data in G.edges(data=True):
        data['distance'] = 1.0 - data['weight']

    # Measure a variety of features of the graph
    strength = dict(G.degree(weight='weight'))  # s_i
    degree = dict(G.degree())  # unweighted degree
    clust_w = nx.clustering(G, weight='weight')  # weighted clustering coeff
    betw = nx.betweenness_centrality(G, weight='distance')  # uses shortest-path distances
    close = nx.closeness_centrality(G, distance='distance')  # distance-based closeness
    pr = nx.pagerank(G, weight='weight')  # weighted PageRank
    tri = nx.triangles(G)  # triangle count per node
    ev = nx.eigenvector_centrality(G, weight='weight', tol=1e-04)  # eigenvector centrality

    # Network features dict
    n_features = {
        'strength': strength,
        'degree': degree,
        'clust_w': clust_w,
        'betw': betw,
        'close': close,
        'pr': pr,
        'tri': tri,
        'ev': ev,
    }

    z_features = {
        'z_strength': robust_z(strength),
        'z_degree': robust_z(degree),
        'z_clust_w': robust_z(clust_w),
        'z_betw': robust_z(betw),
        'z_close': robust_z(close),
        'z_pr': robust_z(pr),
        'z_tri': robust_z(tri),
        'z_ev': robust_z(ev),
    }

    # Simple aggregate anomaly score A_i = sum |z_ik|
    nodes = list(G.nodes())
    A = {}
    for n in nodes:
        A[n] = sum(zf[n] for zf in n_features.values())
    # Add anomaly score as z_features
    z_features['Anomaly_Score'] = A
    z_features['z_Anomaly_Score'] = robust_z(A)
    outdf = df.copy()

    # Set these as columns in df
    for i, row in df.iterrows():
        name = str(i) + "~" + row['Metabolite name']
        # Apply to outdf
        for feat_name, feat_dict in n_features.items():
            outdf.at[i, feat_name] = float(feat_dict.get(name, 0.0))
        for zfeat_name, zfeat_dict in z_features.items():
            outdf.at[i, zfeat_name] = float(zfeat_dict.get(name, 0.0))

        if name in G:
            # Apply to nodes in G
            attrs = {}
            for feat_name, feat_dict in n_features.items():
                attrs[feat_name] = float(feat_dict.get(name, 0.0))
            for zfeat_name, zfeat_dict in z_features.items():
                attrs[zfeat_name] = float(zfeat_dict.get(name, 0.0))
            attrs['Adduct type'] = row['Adduct type']
            G.nodes[name].update(attrs)

    # outdf = detect_outliers_mahalanobis(outdf, features=list(z_features.keys()), verbose=True, sort=False)
    # # Embed mahalanobis scores into graph nodes
    # for i, row in outdf.iterrows():
    #     name = str(i) + "~" + row['Metabolite name']
    #     if name in G:
    #         attrs = {}
    #         attrs['mahalanobis'] = row['mahalanobis']
    #         attrs['outlier'] = row['outlier']
    #         G.nodes[name].update(attrs)
    #     else:
    #         print(f"Node {name} not in graph.", G.nodes())
    return G, outdf




os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")
posfile = "PosIDsRTP2_HeadgroupTailCheck.csv"
negfile = "NegIDsRTP2_HeadgroupTailCheck.csv"
posdf = pd.read_csv(posfile)
negdf = pd.read_csv(negfile)

posdf, posgraphs = class_network_analysis(posdf)
negdf, neggraphs = class_network_analysis(negdf)
graphs = posgraphs + neggraphs

# Save updated dataframes
posdf.to_csv("PosIDsRTP2_HeadgroupTailCheck_Network.csv", index=False)
negdf.to_csv("NegIDsRTP2_HeadgroupTailCheck_Network.csv", index=False)

# Get global min and max degree for color normalization
colortarget = "z_Anomaly_Score"
all_degrees = []
for G in graphs:
    all_degrees.extend([G.nodes[n].get(colortarget, 0.0) for n in G.nodes()])
global_min_degree = np.min(all_degrees)
global_max_degree = np.max(all_degrees)

for G in graphs:
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    pos = nx.spring_layout(G, k=0.5, iterations=50)

    # edges on the specific Axes
    weights = [G[u][v]['weight'] for u, v in G.edges()]
    widths = [max(0.1, float(np.round(w, 2))) for w in weights]
    nx.draw_networkx_edges(G, pos, width=widths, ax=ax)

    # nodes as PathCollection on the specific Axes
    node_list = list(G.nodes())
    sizes = []
    for n in node_list:
        if G.nodes[n].get(colortarget, 0.0) < -2:
            sizes.append(100)  # square
        else:
            sizes.append(500)  # circle
    nc = nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_size=sizes, ax=ax)
    # If colortarget < -2 set the shape to square for that node
    # Set color as degree centrality
    degrees = np.array([G.nodes[n].get(colortarget, 0.0) for n in node_list])
    norm = mpl.colors.Normalize(vmin=np.amin(degrees), vmax=np.amax(degrees))
    cmap = mpl.cm.viridis
    node_colors = cmap(norm(degrees))
    nc.set_facecolor(node_colors)

    # labels on the same Axes
    nx.draw_networkx_labels(G, pos, labels={n: n.split("~")[1] for n in node_list}, font_size=8, ax=ax)

    # prepare label mapping
    labels_text = {node: (
        f"{node.split("~")[1]}\n"
        # f"Betweenness: {G.nodes[node].get('Betweenness', 0.0):.4f}\n"
        f"Score: {G.nodes[node].get(colortarget, 0.0):.4f}\n"
        f"Adduct: {G.nodes[node].get('Adduct type', 'N/A')}\n"
        # f"Closeness: {G.nodes[node].get('Closeness', 0.0):.4f}"
    ) for node in node_list}

    # attach cursor to this collection
    cursor = mplcursors.cursor(nc, hover=True)


    @cursor.connect("add")
    def on_add(sel, node_list=node_list, labels_text=labels_text, nc=nc):
        # robustly extract a single integer index
        idx = getattr(sel, "index", None)
        if idx is None:
            idx = getattr(getattr(sel, "target", None), "index", None)

        # if index is array-like, take first valid element
        if isinstance(idx, (np.ndarray, list, tuple)):
            if len(idx) == 0:
                idx = None
            else:
                try:
                    idx = int(np.atleast_1d(idx)[0])
                except Exception:
                    idx = None

        # fallback: nearest point by mouse coordinates
        if idx is None:
            try:
                offsets = nc.get_offsets()
                mx, my = sel.mouseevent.xdata, sel.mouseevent.ydata
                if mx is None or my is None:
                    return
                d = np.sum((offsets - np.array([mx, my])) ** 2, axis=1)
                idx = int(np.argmin(d))
            except Exception:
                return

        # validate index bounds
        if idx < 0 or idx >= len(node_list):
            return

        node = node_list[int(idx)]
        sel.annotation.set_text(labels_text[node])
plt.show()

