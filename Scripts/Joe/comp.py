from UniDecImporter.ImporterFactory import ImporterFactory
from Scripts.Joe.gaussian import plot_data
import matplotlib.pyplot as plt
import numpy as np



def main():

    # plt.figure(figsize=(10, 6))
    # path = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.mzML"
    # importer = ImporterFactory.create_importer(path)
    # dat = importer.get_single_scan(1)
    # plot_data(dat, label='Control', color='g', subplot_num=222, centroided=False)
    # plt.show()
    plot_gaussian()


def gaussian(x, mu, sigma, height=1):
    return height * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sigma, 2.)))


def plot_gaussian():
    x_axis = np.linspace(400, 2000, 100000)
    unresolved_masses = [57650, 63420, 96500, 234570]
    unresolved_maxcharges = [104, 72, 83, 146]
    unresolved_distwidths = [2.5, 2.5, 2.5, 2.25]
    unresolved_heights = [10000, 14000, 15000, 20000]
    resolution = 100000
    all_charges = np.arange(1, 150)

    intensities = np.zeros_like(x_axis)

    plt.figure(figsize=(18, 12))

    for i in range(len(unresolved_masses)):

        for charge in unresolved_maxcharges:
            zdist = gaussian(all_charges, charge, unresolved_distwidths[i])
            for j in range(len(zdist)):
                if zdist[j] < 0.01:
                    continue
                mz = unresolved_masses[i] / all_charges[j]

                sigma = unresolved_masses[i] / resolution

                gauss = gaussian(x_axis, mz, sigma, unresolved_heights[i] * zdist[j])

                intensities += gauss

    plt.plot(x_axis, intensities)
    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    plt.show()


main()