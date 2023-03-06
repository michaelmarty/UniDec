import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def convolve(ungal):
    ungalcy = np.convolve(ungal[:, 1], ungal[:, 1])
    ungalcx = np.arange(2 * np.amin(ungal[:, 0]), 2 * np.amax(ungal[:, 0]) + 1)
    return np.transpose([ungalcx, ungalcy])


os.chdir("C:\Data\Luis Genentech\\MHC2")


file2 = "MHC2_DPA_tryp"
#file2 = "MHC2_DPA_v8de"
file2 = file2 + "_glycan_list.csv"

df = pd.read_csv(file2)

print(df.keys())

dfs1 = df[df["S1"] > 0]
dfs2 = df[df["S2"] > 0]
dfs3 = df[df["S3"] > 0]
dfs4 = df[df["S4"] > 0]

lens = [len(dfs1), len(dfs2), len(dfs3), len(dfs4)]

indexes = np.array(list(np.ndindex(tuple(lens))))

probs = []
nhexs = []
ngals = []
ngn = []
Fucs = []
nsas = []
for i in indexes:
    r1 = dfs1.iloc[i[0]].fillna(1)
    r2 = dfs2.iloc[i[1]].fillna(1)
    r3 = dfs3.iloc[i[2]].fillna(1)
    r4 = dfs4.iloc[i[3]].fillna(1)
    rs = [r1, r2, r3, r4]

    p = r1["S1"] / 100 * r2["S2"] / 100 * r3["S3"] / 100 * r4["S4"] / 100

    probs.append(float(p))
    nhexs.append(np.sum([r["NHex"] for r in rs]))
    ngals.append(np.sum([r["NGal"] for r in rs]))
    ngn.append(np.sum([r["NGn"] for r in rs]))
    Fucs.append(np.sum([r["Fuc"] for r in rs]))
    nsas.append(np.sum([r["NSa"] for r in rs]))

probs = np.array(probs)
nhexs = np.array(nhexs)
nhexs[nhexs<0] = 0
ngals = np.array(ngals)
ngn = np.array(ngn)
Fucs = np.array(Fucs)
nsas = np.array(nsas)

ugal = np.unique(ngals)
ugal = np.arange(np.amin(ugal), np.amax(ugal)+1)

ugn = np.unique(ngn)
ugn = np.arange(np.amin(ugn), np.amax(ugn)+1)

ufuc = np.unique(Fucs)
ufuc = np.arange(np.amin(ufuc), np.amax(ufuc)+1)

usa = np.unique(nsas)
usa = np.arange(np.amin(usa), np.amax(usa)+1)

#unhex = np.transpose([np.unique(nhexs), [np.sum(probs[np.array(nhexs) == u]) for u in np.unique(nhexs)]])
ungal = np.transpose([ugal, [np.sum(probs[np.array(ngals) == u]) for u in ugal]])
ungn = np.transpose([ugn, [np.sum(probs[np.array(ngn) == u]) for u in ugn]])
unfuc = np.transpose([ufuc, [np.sum(probs[np.array(Fucs) == u]) for u in ufuc]])
unsa = np.transpose([usa, [np.sum(probs[np.array(nsas) == u]) for u in usa]])

if False:
    ungal = convolve(ungal)
    ungn = convolve(ungn)
    unfuc = convolve(unfuc)
    unsa = convolve(unsa)

#plt.plot(unhex[:, 0], unhex[:, 1], label="Hex")
alpha=0.4
plt.bar(unsa[:, 0], unsa[:, 1], label="SA", alpha=alpha)
plt.plot(unsa[:, 0], unsa[:, 1], alpha=alpha)
plt.bar(ungn[:, 0], ungn[:, 1], label="GlcNac", alpha=alpha)
plt.plot(ungn[:, 0], ungn[:, 1],  alpha=alpha)
plt.bar(ungal[:, 0], ungal[:, 1], label="Gal", alpha=alpha)
plt.plot(ungal[:, 0], ungal[:, 1],  alpha=alpha)
plt.bar(unfuc[:, 0], unfuc[:, 1], label="Fuc", alpha=alpha)
plt.plot(unfuc[:, 0], unfuc[:, 1],  alpha=alpha)

outdata = np.array([unsa, ungn, ungal, unfuc], dtype=object)
np.save( file2[:-4]+"_stats.npy",outdata)

plt.legend()
plt.title(file2[:-4])
plt.savefig(file2[:-4]+".png")
plt.savefig(file2[:-4]+".pdf")
plt.show()
