import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# Set backend to Agg
mpl.use('WxAgg')

files = ["1108-10", "1108-11", "1108-12", "1116-10", "1116-11", "1116-12"]

id_rawprsms = [2148,2300,2203,2409,2381,2427]
id_filteredprsms = [997,1047,988,1141,1130,1169]
id_rawmatchedfrags = [45.67,45.92,45.92,51.55,51.97,51.74]
id_filteredmatchedfrags = [54.68,55.84,55.97,63.62,63.25,64.43]
id_rawevals = [0.000128,0.000124,0.000116,0.000134,0.000112,0.00012]
id_filteredevals = [0.0000336,0.0000495,0.0000385,0.0000398,0.0000331,0.0000159]


tfd_rawprsms = [2550,2763,2606,2955,2996,3003]
tfd_filteredprsms = [1291,1381,1273,1528,1560,1577]
tfd_rawmatchedfrags = [38.84,39.01,39.11,43.12,43.19,43.7]
tfd_filteredmatchedfrags = [45.27,46.23,46.73,51.94,51.9,52.4]
tfd_rawevals = [0.000207,0.000191,0.00019,0.000205,0.0002,0.0002]
tfd_filteredevals = [0.000119,0.000116,0.000844,0.0000985,0.0000745,0.0000951]


rawprsms_diffs = np.array(id_rawprsms) - np.array(tfd_rawprsms)
filteredprsms_diffs = np.array(id_filteredprsms) - np.array(tfd_filteredprsms)

n = len(files)
index = np.arange(n)
bar_width = 0.35

fig,ax = plt.subplots()

bars1 = ax.bar(index - bar_width/2, id_filteredmatchedfrags, bar_width, label="IsoDec Raw Matched Fragments")
bars2 = ax.bar(index + bar_width/2, tfd_filteredmatchedfrags, bar_width, label="TopFD Raw Matched Fragments")
ax.set_xlabel("File")
ax.set_ylabel("# Matched Fragments")
ax.set_xticks(index)
ax.set_xticklabels(files)
ax.legend()
ax.set_ylim(0,80)
plt.show()