from random import randint
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.runtime import IsoDecRuntime
import matplotlib.pyplot as plt
from Scripts.Joe.gaussian import *

from UniDecImporter import ImporterFactory
from unidec.modules.isotopetools import calc_averagine_isotope_dist
import numpy as np
from unidec.IsoDec.plots import plot_pks_ax
from unidec.IsoDec import datatools as dt




def create_mix_data(z1, z2, ratio=0.1, noise=0.01, deltamz=0.1, mz=1000):
    if z1 > z2:
        z1, z2 = z2, z1
    mass1 = z1 * mz
    mass2 = z2 * mz

    # Create isotopic distributions
    iso1 = calc_averagine_isotope_dist(mass1, charge=z1, crop=True)
    iso2 = calc_averagine_isotope_dist(mass2, charge=z2, crop=True)

    # Normalize intensities
    iso1[:, 1] = iso1[:, 1] / iso1[:, 1].max()
    iso2[:, 1] = iso2[:, 1] / iso2[:, 1].max() * ratio

    # Add noise
    noise1 = noise * np.random.randn(len(iso1))
    noise2 = noise * np.random.randn(len(iso2))

    # Shift second distribution
    iso2[:, 0] += deltamz

    # Combine distributions
    isodist = np.vstack((iso1, iso2))

    # Resort isodist by m/z
    isodist = isodist[isodist[:, 0].argsort()]

    # If peaks are within 10 ppm, merge them
    newisodist = []
    for d in isodist:
        if len(newisodist) == 0 or abs(d[0] - newisodist[-1][0]) > 0.00001:
            newisodist.append(d)
        else:
            newisodist[-1][1] += d[1]
    isodist = np.array(newisodist)

    return isodist

def decon(isodist):
    eng = IsoDecRuntime()
    eng.config.verbose = False
    pks = eng.batch_process_spectrum(isodist, centroided=True)
    #print(pks)

    parray=[]
    for p in pks:
        parray.append([p.mz, p.peakint, p.z])
    parray = np.array(parray)

    # Sort parray by last column
    if len(parray) > 1:
        parray = parray[parray[:, 2].argsort()]
    else:
        print("None or only 1 peak detected")

    return parray, pks

def check_decon(z1, z2, ratio, deltamz, mz, parray):
    if len(parray) > 2:
        print("Too many peaks detected")
        return False
    if len(parray) < 2:
        print("Too few peaks detected")
        return False
    if parray[0, 2] != z1:
        print("Incorrect first charge state", parray[0, 2], z1)
        return False
    if parray[1, 2] != z2:
        print("Incorrect second charge state", parray[1, 2], z2)
        return False
    meas_ratio = parray[1, 1] / parray[0, 1]
    if abs(meas_ratio - ratio) > 0.1:
        print("Incorrect ratio", ratio, meas_ratio)
        return False
    return True

def create_single_plot(z1, z2, ratio, deltamz, mz):
    isodist = create_mix_data(z1, z2, ratio=ratio, noise=0.01, deltamz=deltamz, mz=mz)
    parray, pks = decon(isodist)
    result = check_decon(z1, z2, ratio, deltamz, mz, parray)
    fig, ax = plt.subplots()
    plot_pks_ax(pks, ax, centroids=isodist, ccolor="k")
    if len(parray)> 0:
        # plt.ylim(-1, 1)
        plt.yticks(np.arange(-1, 1.1, 0.5), fontsize=12)
        plt.hlines(0, np.amin(isodist[:, 0]-0.5), np.amax(isodist[:, 0]+0.5), color="k")
    else:
        # plt.ylim(0, 1)
        plt.yticks(np.arange(0, 1.1, 0.5), fontsize=12)

    plt.xticks(np.arange(0, 1005, 1), fontsize=12)
    plt.xlabel("m/z", fontsize=16)
    plt.xlim(np.amin(isodist[:, 0]-0.25), np.amax(isodist[:, 0]+0.25))

    plt.tight_layout()
    if False:
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\mix_tester2.png",
                    dpi=300)
        plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\mix_tester2.pdf",
                    dpi=300)
    return isodist



def overlap_search(direction, value, dat):
    if value > dat[1,-1] or value < dat[1,0]:
        return
    curr_point = dat[value, 1]
    print(dat[value, 0], dat[value, 1])
    if direction == "left":
        left_check = dat[value-1, 1]
        if left_check > curr_point:
            print("found overlap between", curr_point, " and ", left_check)
            plt.vlines(dat[value-1, 0], 0, np.max(dat[:, 1]), color="darkred", linewidth = 2.5)
        else:
            overlap_search("left", value-1, dat)
    if direction == "right":
        right_check = dat[value+1, 1]
        if right_check > curr_point:
            print("found overlap between", curr_point, " and ", right_check)
            plt.vlines(dat[value-1, 0], 0, -np.max(dat[:, 1]), color="black", linewidth = 2.5)
        else:
            overlap_search("right", value+1, dat)




if __name__ == '__main__':
    # fig, axs = plt.subplots(2, 3, figsize=(6, 4))
    # axs = axs.flatten()
    #
    # for i, r in enumerate(ratios):
    #     create_single_plot(zvals[0], zvals[1], r, delta_mz, mz)
    #     exit()
    #     isodist = create_mix_data(zvals[0], zvals[1], ratio=r, noise=noise_level, deltamz=delta_mz, mz=mz)
    #     parray, pks = decon(isodist)
    #     result = check_decon(zvals[0], zvals[1], r, delta_mz, mz, parray)
    #     ax = axs[i]
    #     plot_pks_ax(pks, ax, centroids=isodist, title=str(result), ccolor="k")
    #     # ax.set_ylim(-1, 1)
    #
    #     labfont = 12
    #     if i > 2:
    #         ax.set_xlabel("m/z", fontsize=labfont)
    #     else:
    #         ax.set_xlabel("")
    #         # clear ticks
    #         ax.set_xticks([])
    #     if i % 3 == 0:
    #         ax.set_ylabel("Intensity", fontsize=labfont)
    #         ax.set_yticks([])
    #     else:
    #         ax.set_ylabel("")
    #         ax.set_yticks([])

    if False:
        ratios = [0.1, 0.5, 0.9, 1, 2, 10]
        zvals = [6, 12]
        noise_level = 0.01
        delta_mz = 0
        mz = 1000

        ratio = 0.75
        dat = create_single_plot(zvals[0], zvals[1], ratio, delta_mz, mz)
        highest = np.amax(dat[:, 1])
        index = np.where(dat == highest)[0]
        overlap_search("left", index, dat)
        overlap_search("right", index, dat)
        plt.tight_layout()
        plt.show()

    if True:
        importer = ImporterFactory.create_importer("Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw")
        dat = importer.grab_centroid_data(1)
        eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False)
        eng.config.activescan = 1
        eng.config.activescanrt = -1
        eng.config.activescanorder = 2

        color_wheel = ['r','g','b','cyan','lime','purple','magenta', 'olive', 'teal']
        eng.batch_process_spectrum(dat, centroided=True)
        mz = []
        z = []
       # plt.ylim([min(reg[:, 1]) - 0.1 * abs(min(reg[:, 1])), max(reg[:, 1]) + 0.1 * abs(max(reg[:, 1]))])
        position = {}
        ranges = {}
        for i in dat:
            plt.vlines(i[0], 0, -i[1], color = 'b')
        highest_max = 0
        for i in eng.pks.peaks:
            max = 0
            for j in i.centroids:
                if j[1] > max:
                    max = j[1]
            #this will get me the starting positions of each cluster
            position[i.mz] = max
            ranges[i.mz] = [i.centroids[0][0], i.centroids[-1][0]]
            if max > highest_max:
                highest_max = max
        for item in position.items():
            random = randint(0, 8)
            plt.vlines(item[0], 0, item[1], color = 'k')
            plt.fill_between([ranges[item[0]][0], ranges[item[0]][1]], 0, highest_max, color = color_wheel[random], alpha = 0.1)

            bin_start, bin_end = ranges[item[0]]
            bin_indices = np.where((dat[:, 0] >= bin_start) & (dat[:, 0] <= bin_end))[0]

            if len(bin_indices) > 0:
                bin_dat = dat[bin_indices]
                mz_value = item[0]
                index = np.where(bin_dat[:, 0] == mz_value)[0]

                if len(index) > 0:
                    index = index[0]
                    overlap_search("right", index, bin_dat)
                    overlap_search("left", index, bin_dat)
                else:
                    print(f"Warning: m/z {mz_value} not found in bin data!")

            else:
                print("bin fail")

        plt.show()






    # plt.plot(mz_values, z_values, color = 'b')
    # # for peak in peaks:
    # #     print(peak)
    # #     start = np.where(dat[:, 0] == peak[1])
    # #     print("start", start)
    # #     overlap_search("left", start, dat)
    # #     overlap_search("right", start, dat)

