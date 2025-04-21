import pandas as pd
import re
import numpy as np
import distance
from sklearn.cluster import AffinityPropagation
import networkx as nx
import matplotlib.pyplot as plt
from sklearn import metrics
import matplotlib.cm as cm


file = "D:\\Data\\Lipidomics\\OpenMS\\FilteredFixedAll.csv"
file = "D:\Data\Lipidomics\Exported_TLs\\220816_MS1_Ecoli_lipid_ND_res_outputs_1_filtered.csv"

df = pd.read_csv(file)
print(df.keys())
names = df["Full Name"].to_numpy()
classes = df["Ontology"]
df2 = df[classes=="PG"]
#df2 = df2[df2["Tail Unsaturation"]==0]
#print(df2["Full Name"])

colors = cm.rainbow(df2["Tail Unsaturation"].to_numpy()/4)

mz = df2["mass"]
RT = df2["Average Rt(min)"]
cRT = RT /df2["Tail Lengths"]
cRT2 = RT
md = df2["kmassd"]

plt.subplot(121)
plt.scatter(mz, cRT2, color=colors)
plt.subplot(122)
plt.scatter(mz, md, color=colors)
plt.show()


exit()

lev_similarity = -1*np.array([[distance.levenshtein(w1,w2) for w1 in names] for w2 in names])

print(lev_similarity)

affprop = AffinityPropagation(affinity="precomputed", damping=0.9)
affprop.fit(lev_similarity)
cluster_centers_indices = affprop.cluster_centers_indices_
labels = affprop.labels_

n_clusters_ = len(cluster_centers_indices)

print("Estimated number of clusters: %d" % n_clusters_)

for cluster_id in np.unique(affprop.labels_):
    exemplar = names[affprop.cluster_centers_indices_[cluster_id]]
    cluster = np.unique(names[np.nonzero(affprop.labels_==cluster_id)])
    cluster_str = ", ".join(cluster)
    print(" - *%s:* %s" % (exemplar, cluster_str))

from itertools import cycle

plt.close("all")
plt.figure(1)
plt.clf()

colors = cycle("bgrcmykbgrcmykbgrcmykbgrcmyk")
for k, col in zip(range(n_clusters_), colors):
    class_members = labels == k
    cluster_center = lev_similarity[cluster_centers_indices[k]]
    plt.plot(lev_similarity[class_members, 0], lev_similarity[class_members, 1], col + ".")
    #plt.plot(cluster_center[0],cluster_center[1],"o",markerfacecolor=col, markeredgecolor="k",markersize=14,)
    for x in lev_similarity[class_members]:
        plt.plot([cluster_center[0], x[0]], [cluster_center[1], x[1]], col)

plt.title("Estimated number of clusters: %d" % n_clusters_)
plt.show()