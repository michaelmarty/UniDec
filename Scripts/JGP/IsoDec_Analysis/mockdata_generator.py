from unidec.modules.isotopetools import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from unidec.IsoDec.datatools import *
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
import unidec.tools as ud
from unidec.IsoDec.match import *

def generate_gaussian(x_axis, mu, sigma, height=1):
    return height * np.exp(-np.power(x_axis - mu, 2.) / (2 * np.power(sigma, 2.)))

def get_zdist(all_charges, peakcharge, sigma):
    zdist = generate_gaussian(all_charges, peakcharge, sigma)
    #Normalize the distribution to max out at 1
    zdist = zdist / np.max(zdist)
    #Set anything below 0.1 to 0
    zdist[zdist < 0.1] = 0
    return zdist

def add_isodist(mzs, intensities, mass, charge, intensity,
                min_mz, max_mz, resolution,
                annotated_masses=None, annotated_mzs=None, annotated_zs=None):
    dist = fast_calc_averagine_isotope_dist(mass, charge)

    initial_length = len(mzs)
    #Scale the distribution to the intensity
    dist[:,1] = dist[:,1] * intensity
    #Check if both the dist min mz and max mz are within the mz range
    if dist[0,0] > max_mz or dist[-1,0] < min_mz:
        return mzs, intensities

    if annotated_masses is not None and annotated_mzs is not None and annotated_zs is not None:
        annotated_masses.append(mass)
        #Find the mz at the maximum intensity
        apex_mz = dist[np.argmax(dist[:,1]),0]
        annotated_mzs.append(apex_mz)
        annotated_zs.append(charge)


    #Now see if we can add each peak in this distribution to the data
    for i in range(dist.shape[0]):
        if dist[i,0] < min_mz or dist[i,0] > max_mz:
            continue
        if len(mzs) == 0:
            mzs.append(dist[i,0])
            intensities.append(dist[i,1])
        else:
            nearest_idx = fastnearest(np.array(mzs), dist[i,0])
            delta_mz = np.abs(dist[i,0] - mzs[nearest_idx])
            theo_res = dist[i,0] / (np.abs(dist[i,0]-mzs[nearest_idx]))

            if theo_res > resolution or delta_mz == 0:
                intensities[nearest_idx] += dist[i,1]
            else:
                #Check if the current mz is greater than the nearest mz
                if dist[i,0] > mzs[nearest_idx] and nearest_idx < len(mzs)-1:
                    #Insert the new mz and intensity at the nearest index + 1
                    mzs.insert(nearest_idx+1, dist[i,0])
                    intensities.insert(nearest_idx+1, dist[i,1])
                elif dist[i,0] > mzs[nearest_idx] and nearest_idx == len(mzs)-1:
                    mzs.append(dist[i,0])
                    intensities.append(dist[i,1])
                else:
                    mzs.insert(nearest_idx, dist[i,0])
                    intensities.insert(nearest_idx, dist[i,1])

    final_length = len(mzs)
    if initial_length == final_length:
        print("No peaks added for mass: ", mass, "charge: ", charge)
    return mzs, intensities


def generate_mockdata(mass_count, min_mz, max_mz, resolution, all_charges):
    annotated_masses = []
    annotated_mzs = []
    annotated_zs = []

    # I want a list of N random masses between 1000 and 35000
    masses = np.random.uniform(1000, 35000, mass_count)
    charges = np.zeros_like(masses)
    dist_intensities = np.random.uniform(1000, 1000000, mass_count)
    zdist_widths = np.random.uniform(1, 5, mass_count)

    # Go through each mass and select a random charge state that will fall within the mz range
    for i, mass in enumerate(masses):
        max_charge = np.ceil(mass / min_mz)
        min_charge = np.floor(mass / max_mz)

        if max_charge > all_charges[-1]:
            max_charge = all_charges[-1]
        if min_charge < all_charges[0]:
            min_charge = all_charges[0]

        if min_charge == max_charge or min_charge > max_charge:
            charges[i] = min_charge
        else:
            charges[i] = np.random.randint(min_charge, max_charge + 1)

    mzs = []
    intensities = []

    for i in range(mass_count):
        # Create an isotopic distribution for each charge state of each mass
        zdist = get_zdist(all_charges, charges[i], zdist_widths[i])
        for z, z_int in enumerate(zdist):
            if z_int == 0 or z == 0:
                continue
            else:
                mzs, intensities = add_isodist(mzs, intensities, masses[i], all_charges[z], dist_intensities[i] * z_int,
                                               min_mz, max_mz, resolution, annotated_masses, annotated_mzs,
                                               annotated_zs)

    if len(mzs) == 0:
        print("No peaks generated")
        return None, None, None, None

    # Zip the mzs and intensities together
    centroids = np.array(list(zip(mzs, intensities)))
    # Sort the centroids by mz
    centroids = centroids[centroids[:, 0].argsort()]

    # for i in range(len(mzs)):
    #     plt.plot([mzs[i], mzs[i]], [0, intensities[i]], 'k-')
    # plt.show()

    return centroids, annotated_masses, annotated_mzs, annotated_zs

def process_mockdata(centroids, eng):

    eng.batch_process_spectrum(centroids, centroided=True)

    if False:
        matched_experimental, allz, unmatchedz = compare_experimental_peaks_to_annotated(eng.pks, annotated_mzs, annotated_zs, ppm_tol=5)

        matched_exp = [eng.pks.peaks[i] for i in matched_experimental]
        matchedpeaks_collection = MatchedCollection()
        matchedpeaks_collection.add_peaks(matched_exp)

        unmatched_exp = [eng.pks.peaks[i] for i in range(len(eng.pks.peaks)) if i not in matched_experimental]
        unmatchedpeaks_collection = MatchedCollection()
        unmatchedpeaks_collection.add_peaks(unmatched_exp)

        plot_pks(unmatchedpeaks_collection, centroids=centroids, forcecolor="r", show=False, plotmass=False,
                 nocentroids=True, labelz=True)
        plot_pks(matchedpeaks_collection, centroids=centroids, ccolor='k', forcecolor="g", show=False, plotmass=False,
                 nocentroids=True, labelz=True)
        plt.show()

    return eng.pks

def compare_experimental_peaks_to_annotated(experimental_mc,
                                            annotated_mzs, annotated_zs,
                                            ppm_tol=20, verbose=False):
    annotated_mzs_matched = []
    #Go through each peak in the experimental peaks and see if it matches an annotated peak
    #Order the annotated_masses, annotated_mzs, and annotated_zs by mz
    annotated_mzs = np.array(annotated_mzs)
    annotated_zs = np.array(annotated_zs)
    sorted_indices = np.argsort(annotated_mzs)
    annotated_mzs = annotated_mzs[sorted_indices]
    annotated_zs = annotated_zs[sorted_indices]

    matched_experimental = []

    for i in range(len(experimental_mc.peaks)):
        nearest_idx = fastnearest(annotated_mzs, experimental_mc.peaks[i].mz)
        allowed_mzs = [experimental_mc.peaks[i].mz]
        shift = 1.0033 / experimental_mc.peaks[i].z
        for j in range(1,3):
            allowed_mzs.append(experimental_mc.peaks[i].mz + j*shift)
            allowed_mzs.append(experimental_mc.peaks[i].mz - j*shift)

        if any([ud.within_ppm(x, annotated_mzs[nearest_idx], ppm_tol) for x in allowed_mzs]):
            #Now check if the charge state is the same
            if experimental_mc.peaks[i].z == annotated_zs[nearest_idx]:
                annotated_mzs_matched.append(nearest_idx)
                matched_experimental.append(i)
            else:
                if verbose:
                    print("Charge state mismatch: ", experimental_mc.peaks[i].mz, experimental_mc.peaks[i].z, annotated_mzs[nearest_idx], annotated_zs[nearest_idx])

    if verbose:
        print("Total annotated peaks: ", len(annotated_mzs))
        print("Total annotated peaks matched: ", len(annotated_mzs_matched))
        print("Percentage of annotated peaks matched: ", len(annotated_mzs_matched)/len(annotated_mzs)*100)

        print("Total experimental peaks: ", len(experimental_mc.peaks))
        print("Total experimental peaks matched: ", len(annotated_mzs_matched))
        print("Percentage of experimental peaks matched: ", len(annotated_mzs_matched)/len(experimental_mc.peaks)*100)

    unmatched_mzs = [annotated_mzs[i] for i in range(len(annotated_mzs)) if i not in annotated_mzs_matched]
    unmatched_zs = [annotated_zs[i] for i in range(len(annotated_zs)) if i not in annotated_mzs_matched]

    allz_dictionary = {}
    unmatchedz_dictionary = {}

    for z in annotated_zs:
        if z in allz_dictionary:
            allz_dictionary[z] += 1
        else:
            allz_dictionary[z] = 1

    for z in unmatched_zs:
        if z in unmatchedz_dictionary:
            unmatchedz_dictionary[z] += 1
        else:
            unmatchedz_dictionary[z] = 1

    #Order the dictionaries by key
    allz_dictionary = dict(sorted(allz_dictionary.items()))
    unmatchedz_dictionary = dict(sorted(unmatchedz_dictionary.items()))

    return matched_experimental, allz_dictionary, unmatchedz_dictionary

def sum_dictionaries(existing_zdict, current_zdict):
    for z in current_zdict:
        if z in existing_zdict:
            existing_zdict[z] += current_zdict[z]
        else:
            existing_zdict[z] = current_zdict[z]
    return existing_zdict

def generate_mockdata_equalz(mass_count, resolution, all_charges):
    ann_masses = []
    ann_mzs = []
    ann_zs = []


    # I want a list of N random masses between 1000 and 35000
    masses = np.random.uniform(10000, 35000, mass_count)
    dist_intensities = np.random.uniform(1000, 1000000, mass_count)

    print("Masses: ", masses)
    print("Intensities: ", dist_intensities)

    mzs = []
    intensities = []

    for i in range(mass_count):
        for j in all_charges:

            mzs, intensities = add_isodist(mzs, intensities, masses[i], j, dist_intensities[i],
                                               0, 38000, resolution, ann_masses, ann_mzs,
                                               ann_zs)

    if len(mzs) == 0:
        print("No peaks generated")
        return None, None, None, None


    # max_intensity = np.max(intensities)
    # #add random noise from the minimum centroid mz - 10 to the maximum centroid mz + 10
    # noise_mzs = np.random.uniform(mzs[0]-10, mzs[-1]+10, 10000)
    # noise_intensities = np.random.uniform(0, max_intensity*0.01, 10000)
    # #If any of the noise mzs are within the mz range of the centroids, remove them
    # for i in range(len(noise_mzs)):
    #     nearest_idx = fastnearest(np.array(mzs), noise_mzs[i])
    #     delta_mz = np.abs(noise_mzs[i] - mzs[nearest_idx])
    #     theo_res = noise_mzs[i] / (np.abs(noise_mzs[i] - mzs[nearest_idx]))
    #     if theo_res > resolution and delta_mz != 0:
    #         intensities[nearest_idx] += noise_intensities[i]
    #     else:
    #         #Insert the mz and intensity at the proper location
    #         if mzs[nearest_idx] < noise_mzs[i]:
    #             mzs.insert(nearest_idx+1, noise_mzs[i])
    #             intensities.insert(nearest_idx+1, noise_intensities[i])
    #         else:
    #             mzs.insert(nearest_idx, noise_mzs[i])
    #             intensities.insert(nearest_idx, noise_intensities[i])

    # for i in range(len(mzs)):
    #     plt.plot([mzs[i], mzs[i]], [0, intensities[i]], 'k-')
    # plt.show()

    # Zip the mzs and intensities together
    centroids = np.array(list(zip(mzs, intensities)))
    # Sort the centroids by mz
    centroids = centroids[centroids[:, 0].argsort()]

    return centroids, ann_masses, ann_mzs, ann_zs

def generate_mockdata_singlez(charge, min_mz=200, max_mz=2000):
    min_mass = (min_mz+20) * charge
    max_mass = (max_mz-20) * charge

    mass = np.random.uniform(min_mass, max_mass)
    intensity = np.random.uniform(1000, 1000000)
    dist = fast_calc_averagine_isotope_dist(mass, charge)
    mzs = []
    intensities = []
    for i in range(dist.shape[0]):
        if dist[i,0] < min_mz or dist[i,0] > max_mz:
            continue
        mzs.append(dist[i,0])
        intensities.append(dist[i,1]*intensity)
    centroids = np.array(list(zip(mzs, intensities)))
    if len(centroids) < 3:
        return []

    centroids = centroids[centroids[:, 0].argsort()]


    return centroids



# Set backend to Agg
mpl.use('WxAgg')
centroids = None
if True:
    min_mz = 0
    max_mz = 37000
    resolution = 240000
    mass_count = 1
    all_charges = np.arange(1,50)

    annotated_masses = []


    loop_count = 10000

    existing_allz_dict = {}
    existing_unmatchedz_dict = {}

    for i in all_charges:
        existing_allz_dict[i] = 0
        existing_unmatchedz_dict[i] = 0

    eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False)
    eng.config.activescan = 1
    eng.config.activescanrt = -1
    eng.config.activescanorder = 2

    eng.config.knockdown_rounds = 5
    eng.config.css_thresh = 0.7
    eng.config.matchtol = 5
    eng.predmode = 0
    start_time = time.perf_counter()
    for i in range(loop_count):
        if i % 10 == 0:
            print("Iteration: ", i)
            current_time = time.perf_counter()
            print("Time elapsed: ", current_time - start_time)
            start_time = current_time

        # centroids, annotated_masses, annotated_mzs, annotated_zs = (
        #     generate_mockdata(mass_count, min_mz, max_mz, resolution, all_charges))
        # centroids, annotated_masses, annotated_mzs, annotated_zs = (
        #         generate_mockdata_equalz(mass_count, resolution, all_charges)
        # )
        for z in all_charges:
            centroids = generate_mockdata_singlez(z)
            if len(centroids) >= 3:
                annotated_mzs = []
                annotated_zs = []

                annotated_zs.append(z)
                annotated_mzs.append(centroids[np.argmax(centroids[:,1]),0])

                mc = process_mockdata(centroids,eng)
                matched_experimental, allz_dictionary, unmatchedz_dictionary = (
                    compare_experimental_peaks_to_annotated(mc, annotated_mzs, annotated_zs, ppm_tol=5))

                existing_allz_dict = sum_dictionaries(existing_allz_dict, allz_dictionary)
                existing_unmatchedz_dict = sum_dictionaries(existing_unmatchedz_dict, unmatchedz_dictionary)

                mc.peaks = []
            else:
                print("Insufficient peaks!")


    fracunmatched_zdict = {}
    for z in existing_allz_dict.keys():
        unmatched = existing_unmatchedz_dict[z] if z in existing_unmatchedz_dict else 0
        total = existing_allz_dict[z] if z in existing_allz_dict else 0
        if total == 0:
            fracunmatched_zdict[z] = 0
        else:
            fracunmatched_zdict[z] = unmatched / total

    #Order dictionaries by key
    existing_allz_dict = dict(sorted(existing_allz_dict.items()))
    existing_unmatchedz_dict = dict(sorted(existing_unmatchedz_dict.items()))

    if True:
        plt.plot(list(existing_allz_dict.keys()), list(existing_allz_dict.values()), 'bo-', label="All Synthetic Clusters")
        plt.plot(list(existing_unmatchedz_dict.keys()), list(existing_unmatchedz_dict.values()), 'ro-', label="Unmatched Clusters")
        plt.legend(bbox_to_anchor=(1, 1.2))
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xlabel("Charge State", fontsize=18)
        plt.xticks(fontsize=14)
        plt.ylabel("Count", fontsize=18)
        plt.yticks(fontsize=14)
        #Set the y minimum to 0
        plt.ylim(0, max(list(existing_allz_dict.values()))*1.05)
        plt.xlim(0,50)
        plt.tight_layout()
        if True:
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\synthdata_zdist_nooverlap_newmodel.png", dpi=300)
            plt.savefig("D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\synthdata_zdist_nooverlap_newmodel.pdf", dpi=300)
        plt.show()

    if True:
        plt.plot(fracunmatched_zdict.keys(), fracunmatched_zdict.values(), 'ro-', clip_on=False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.xlabel("Charge State", fontsize=18)
        plt.xticks(fontsize=14)
        plt.ylabel("Fraction of\nUnmatched Clusters", ha='center', fontsize=18)
        plt.yticks(fontsize=14)
        # Set the y minimum to 0
        plt.ylim(0, 1)
        plt.xlim(0, 50)
        plt.tight_layout()
        if True:
            plt.savefig(
                "D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\synthdata_fdist_nooverlap_newmodel.png",
                dpi=300)
            plt.savefig(
                "D:\\Johnny\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\synthdata_fdist_nooverlap_newmodel.pdf",
                dpi=300)
        plt.show()
    exit()

    np.savetxt("Z:\\Group Share\\JGP\\MockData\\mockdata_centroids.csv", centroids, delimiter=",")
    np.savetxt("Z:\\Group Share\\JGP\\MockData\\annotated_mzs.csv", annotated_mzs, delimiter=",")
    np.savetxt("Z:\\Group Share\\JGP\\MockData\\annotated_zs.csv", annotated_zs, delimiter=",")

if centroids is None:
    centroids = np.loadtxt("Z:\\Group Share\\JGP\\MockData\\mockdata_centroids.csv", delimiter=",")
    annotated_mzs = np.loadtxt("Z:\\Group Share\\JGP\\MockData\\annotated_mzs.csv", delimiter=",")
    annotated_zs = np.loadtxt("Z:\\Group Share\\JGP\\MockData\\annotated_zs.csv", delimiter=",")

mzrange = [1114,1116]
mzrange = None

if mzrange is not None:
    centroids = ud.datachop(centroids, mzrange[0], mzrange[1])

eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=True)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

eng.config.knockdown_rounds = 5
eng.config.css_thresh = 0.7
eng.config.matchtol = 5

eng.batch_process_spectrum(centroids, centroided=True)

matched_experimental =  compare_experimental_peaks_to_annotated(eng.pks, annotated_mzs, annotated_zs, ppm_tol=5)

matched_exp = [eng.pks.peaks[i] for i in matched_experimental]
matchedpeaks_collection = MatchedCollection()
matchedpeaks_collection.add_peaks(matched_exp)

unmatched_exp = [eng.pks.peaks[i] for i in range(len(eng.pks.peaks)) if i not in matched_experimental]
unmatchedpeaks_collection = MatchedCollection()
unmatchedpeaks_collection.add_peaks(unmatched_exp)

plot_pks(unmatchedpeaks_collection, centroids=centroids, forcecolor="r", show=False, plotmass=False, nocentroids=True, labelz=True)
plot_pks(matchedpeaks_collection, centroids=centroids, ccolor='k', forcecolor="g", show=False, plotmass=False, nocentroids=True, labelz=True)
plt.show()







