import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from unidec.IsoDec.encoding import *

# Set backend to Agg
mpl.use('WxAgg')
os.chdir("C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures")
path = "Z:\\Group Share\\JGP\\PXD019247\\20170105_L_MaD_ColC4_Ecoli20161108-12_01.pkl"

tr, te, zd = encode_phase_file(path)

plt.figure()

# turn off axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
z = 9
for t in tr:
    emat = t[0]
    centroids = t[1]
    zval = t[2]

    if zval == z:
        centroids[:, 1] /= np.amax(centroids[:, 1])

        for c in centroids:
            plt.vlines(c[0], 0, c[1], color="k")
        plt.yticks([])
        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.savefig("training_data9.png", dpi=300, transparent=True)
        plt.savefig("training_data9.pdf")

        xaxis = np.arange(np.amin(centroids[:, 0]), np.amax(centroids[:, 0]), 0.001)
        yaxis = np.sin(xaxis*(2 * np.pi)*z / 1.0033)/2 - 0.5

        plt.figure()
        np.savetxt("emat.txt", emat)
        plt.imshow(emat, aspect="auto", cmap="viridis", interpolation="none")
        # Set ytick marks to be the charge states
        plt.yticks(range(0, len(emat), 10), range(1, len(emat) + 1, 10))
        plt.xlabel("Phase Bin")
        plt.ylabel("Charge")
        plt.savefig("phaseimage.png", dpi=300, transparent=True)
        plt.savefig("phaseimage.pdf")
        plt.show()
        exit()
        plt.figure()
        # find all places where sine crosses zero
        zeros = np.where(np.diff(np.sign(yaxis+0.5)))[0]

        # find all places where sine crosses zero and the derivative is positive
        pos = np.where(np.diff(np.sign(yaxis+0.5)) > 0)[0]
        # find all places where sine crosses zero and the derivative is negative
        neg = np.where(np.diff(np.sign(yaxis+0.5)) < 0)[0]

        # Final halfway between pos and neg
        p2 = np.array([pos, neg]).mean(axis=0).astype(int)
        # find halfway between pos and neg again
        p3 = np.array([pos[1:], neg[:-1]]).mean(axis=0).astype(int)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.plot(xaxis, yaxis, color="b")
        for zero in pos:
            plt.vlines(xaxis[zero], -1, 1, color="r")

        for p in p2:
            plt.vlines(xaxis[p], -1, 1, color="purple")

        for p in p3:
            plt.vlines(xaxis[p], -1, 1, color="yellow")

        for p in neg:
            plt.vlines(xaxis[p], -1, 1, color="green")

        plt.xlim(602, 603)
        plt.xlabel("m/z")
        plt.title("+" + str(z))
        plt.savefig("phase2.png", dpi=300, transparent=True)
        plt.savefig("phase2.pdf")
        plt.show()
        break


exit()
zvals = range(1, 31)

for z in zvals:
    plt.subplot(10, 3, z)
    plt.ylabel("+" + str(z))
    plt.xticks([])
    plt.yticks([])
    # turn off axis
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    for t in tr:
        emat = t[0]
        centroids = t[1]
        zval = t[2]

        if zval == z:
            plt.imshow(emat, aspect="auto", cmap="viridis", interpolation="none")
            break

plt.savefig("encoding.png", dpi=300, transparent=True)
plt.savefig("encoding.pdf")
plt.tight_layout()
plt.show()
