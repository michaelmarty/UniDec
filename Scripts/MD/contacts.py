import os
import matplotlib.pyplot as plt
from Scripts.MD.tools import *
import h5py
import sys

sys.path.append("/mnt/c/Python/unidec")

topdir = "C:\\VMD\\AmtB_Nanodiscs"
topdir = "Z:\\mtmarty\\VMD"
dirs = ["POPG_E3_4NH2_all_atom", "POPC_E3_4NH2_all_atom", "PGPC_E3_4NH2_all_atom"]
dirs = ["PGPC_E3_4NH2_all_atom"]
colors = ["b", "b", "g"]

edges = np.arange(0, 20, 0.1)
mids = np.array([np.mean(edges[x - 1:x + 1]) for x in range(1, len(edges))])

plt.figure()
for i, d in enumerate(dirs):
    path = os.path.join(topdir, d, "analysis.hdf5")
    hdf = h5py.File(path)

    agroup = hdf.get("rmsf")
    resnames = get_dataset(agroup, "resnames")
    barray = np.array([b"POPC" in r for r in resnames])
    parray = barray[barray]

    agroup = hdf.get("all_atoms")
    md = get_dataset(agroup, "min_dists")
    dat = np.ravel(md)
    hist = np.histogram(dat, bins=edges)[0] / float(len(dat))
    plt.subplot(311)
    cutoffs, peaks = find_cutoff(mids, hist, window=10, threshold=0.05)
    counts, percents = cutoff_count(md, cutoffs, barray)
    print("All Atom:",counts, percents)
    mval = np.amax(hist)
    #for p in peaks:
    #    plt.vlines(p, 0, mval, color=colors[i])
    for j, cutoff in enumerate(cutoffs):
        plt.vlines(cutoff, 0, mval)
        count = int(round(counts[j]))
        plt.text(cutoff, mval * 1.05, str(count))
    plt.plot(mids, hist, color=colors[i])

    agroup = hdf.get("all_noh")
    md = get_dataset(agroup, "min_dists")
    #print(md.shape)
    dat = np.ravel(md)
    hist = np.histogram(dat, bins=edges)[0] / float(len(dat))
    cutoffs, peaks = find_cutoff(mids, hist)
    counts, percents = cutoff_count(md, cutoffs, barray)
    print("All Atom: NoH", counts, percents)
    mval = np.amax(hist)
    #for p in peaks:
    #    plt.vlines(p, 0, mval, color=colors[i])
    for j, cutoff in enumerate(cutoffs):
        if cutoff < 10:
            plt.vlines(cutoff, 0, mval, linestyle="--")
            count = int(round(counts[j]))
            plt.text(cutoff, mval * 1.05, str(count))
    # plt.subplot(311)
    plt.plot(mids, hist, color=colors[i], linestyle="--")

    agroup = hdf.get("lipid_cation")
    md = get_dataset(agroup, "min_dists")
    dat = np.ravel(md)
    hist = np.histogram(dat, bins=edges)[0] / float(len(dat))
    plt.subplot(312)
    cutoffs, peaks = find_cutoff(mids, hist, window=10)
    counts, percents = cutoff_count(md, cutoffs, parray)
    print("Lipid Cation:", counts, percents)
    mval = np.amax(hist)
    # for p in peaks:
    #    plt.vlines(p, 0, mval, color=colors[i])
    for j, cutoff in enumerate(cutoffs):
        if cutoff < 10:
            plt.vlines(cutoff, 0, mval, linestyle="-")
            count = int(round(counts[j]))
            plt.text(cutoff, mval * 1.05, str(count))
    plt.plot(mids, hist, color=colors[i])

    agroup = hdf.get("lipid_anion")
    md = get_dataset(agroup, "min_dists")
    dat = np.ravel(md)
    hist = np.histogram(dat, bins=edges)[0] / float(len(dat))
    #plt.subplot(312)
    cutoffs, peaks = find_cutoff(mids, hist, window=10)
    counts, percents = cutoff_count(md, cutoffs, barray)
    print("Lipid Anion:", counts, percents)
    mval = np.amax(hist)
    # for p in peaks:
    #    plt.vlines(p, 0, mval, color=colors[i])
    for j, cutoff in enumerate(cutoffs):
        if cutoff < 10:
            plt.vlines(cutoff, 0, mval, linestyle="--")
            count = int(round(counts[j]))
            plt.text(cutoff, mval * 0.8, str(count))
    plt.plot(mids, hist, color=colors[i], linestyle="--")


    agroup = hdf.get("lipid_hdon")
    md = get_dataset(agroup, "min_dists")
    dat = np.ravel(md)
    hist = np.histogram(dat, bins=edges)[0] / float(len(dat))
    plt.subplot(313)

    cutoffs, peaks = find_cutoff(mids, hist, window=10)
    counts, percents = cutoff_count(md, cutoffs, barray)
    print("Lipid H Donor:", counts, percents)
    mval = np.amax(hist)
    # for p in peaks:
    #    plt.vlines(p, 0, mval, color=colors[i])
    for j, cutoff in enumerate(cutoffs):
        if cutoff < 10:
            plt.vlines(cutoff, 0, mval, linestyle="-")
            count = int(round(counts[j]))
            plt.text(cutoff, mval * 1.05, str(count))

    plt.plot(mids, hist, color=colors[i])

    agroup = hdf.get("lipid_hacc")
    md = get_dataset(agroup, "min_dists")
    dat = np.ravel(md)
    hist = np.histogram(dat, bins=edges)[0] / float(len(dat))
    #plt.subplot(313)

    cutoffs, peaks = find_cutoff(mids, hist, window=10)
    counts, percents = cutoff_count(md, cutoffs, barray)
    print("Lipid H Acceptor:", counts, percents)
    mval = np.amax(hist)
    # for p in peaks:
    #    plt.vlines(p, 0, mval, color=colors[i])
    for j, cutoff in enumerate(cutoffs):
        if cutoff < 10:
            plt.vlines(cutoff, 0, mval, linestyle="--")
            count = int(round(counts[j]))
            plt.text(cutoff, mval * 1.05, str(count))
    plt.plot(mids, hist, color=colors[i], linestyle="--")

    agroup = hdf.get("lipid_heads")
    md = get_dataset(agroup, "min_dists")
    dat = np.ravel(md)
    hist = np.histogram(dat, bins=edges)[0] / float(len(dat))
    #plt.subplot(313)

    cutoffs, peaks = find_cutoff(mids, hist, window=10)
    counts, percents = cutoff_count(md, cutoffs, barray)
    print("Head Groups:", counts, percents)
    mval = np.amax(hist)
    # for p in peaks:
    #    plt.vlines(p, 0, mval, color=colors[i])
    for j, cutoff in enumerate(cutoffs):
        if cutoff < 10:
            plt.vlines(cutoff, 0, mval, linestyle="-.")
            count = int(round(counts[j]))
            plt.text(cutoff, mval * 1.05, str(count))
    plt.plot(mids, hist, color=colors[i], linestyle="-.")

    hdf.close()
plt.show()
