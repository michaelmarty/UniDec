import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

dir = "Z:\\Group Share\\Skippy\\Projects\\HT\\Good data sets\\Bgal bit len"

lens = [1, 3, 5]
colors = ["r", "g", "b"]

plt.figure(figsize=(3, 2), layout='constrained')
#fig, axes = plt.subplots(2, 1, sharey=True, figsize=(6.5,3.25))
#plt.ylabel("Intensity (# of Ions)", loc='center')

for len in lens:
    top_folder = os.path.join(dir, str(len) + " bit")

    chromfiles = ud.match_files_recursive(top_folder, "chrom.txt")
    cdata = [np.loadtxt(file) for file in chromfiles]

    files = ud.match_files_recursive(top_folder, "TIC_decon.txt")
    print(files)
    datas = [np.loadtxt(file) for file in files]

    if len == 1:
        datas = cdata

    plt.subplot(212)
    for d in datas:
        plt.plot(d[:, 0], d[:, 1], color=colors[lens.index(len)], label=str(len) + " bit")
        print(d.shape)


    plt.subplot(211)
    for d in cdata:
        plt.plot(d[:, 0], d[:, 1], color=colors[lens.index(len)], label=str(len) + " bit")
        print(d.shape)

plt.subplot(211)
plt.xlim(0, 41)
ax = plt.gca()
ax.set_xticklabels([])
ax.spines[['right', 'top']].set_visible(False)
plt.ylabel("(# of Ions)", labelpad=9, loc='bottom')
# add legend using only every third entry
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[::3], labels[::3], loc='upper left')
#plt.legend(loc='upper left')
#plt.axis('off')
#ax.set_xticklabels([])

plt.subplot(212)
plt.xlim(8, 12)
ax = plt.gca()
ax.spines[['right', 'top']].set_visible(False)
plt.xlabel("Time (min)")
plt.ylabel("Intensity", labelpad=4, loc='top')
plt.show()
