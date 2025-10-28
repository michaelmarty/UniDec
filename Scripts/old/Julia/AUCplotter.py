import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

dir = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\M2"
os.chdir(dir)

file = "AUC data_MTM.xlsx"

df = pd.read_excel(file)

print(df)

dets = ["c8e4", "LDAO"]

reps = [1,2,3]

ax1 = plt.subplot(121)
ax2 = plt.subplot(122)
axes = [ax1, ax2]

colors=["maroon", "r", "firebrick", "b", "royalblue" , "darkblue"]

for i, d in enumerate(dets):
    for j, r in enumerate(reps):
        label = d + " " + str(r)

        svals = df[label+" s"]
        ivals = df[label+" i"]

        ivals /= np.amax(ivals)

        b1 = svals < 7.5
        ax=axes[i]
        ax.plot(svals[b1], ivals[b1], color=colors[3*i+j])

        if i ==0:
            ax.set_ylabel("Normalized Intensity")

        ax.set_xlabel("Sedimentation Coefficient (S)")
        ax.set_title(d.upper() + " pH 5")

plt.savefig("AUCfig.pdf")
plt.savefig("AUCfig.png")
plt.show()



