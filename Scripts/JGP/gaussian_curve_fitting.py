import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from UniDecImporter.ImporterFactory import ImporterFactory
import unidec.tools as ud
import unidec.IsoDec.datatools as dt
import matplotlib as mpl

# Set backend to Agg
mpl.use('WxAgg')


def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))


def multi_gaussian(x, *params):
    # Number of Gaussians is determined by the number of parameters
    n_gaussians = len(params) // 3
    result = np.zeros_like(x)

    # Each Gaussian has 3 parameters: A, mu, sigma
    for i in range(n_gaussians):
        A = params[i * 3]
        mu = params[i * 3 + 1]
        sigma = params[i * 3 + 2]
        result += gaussian(x, A, mu, sigma)
    return result

def get_peak_groups(mz, intensities, peaks, resolved_threshold=0.05):
    #Group peaks where the intensity does not drop below a relative threshold from the nearest peak between them
    peak_groups = []
    #Order the peaks by decreasing m/z
    peaks = sorted(peaks, key=lambda x: x[0])

    current_group = []
    for p in peaks:
        if len(current_group) == 0:
            current_group.append(p)
        else:
            current_idx = dt.fastnearest(mz, current_group[-1][0])
            next_idx = dt.fastnearest(mz, p[0])

            min_ratio = 1
            #Check every intensity between the current idx and the next idx
            for i in range(current_idx, next_idx + 1):
                ratio = intensities[i] / intensities[current_idx]
                if ratio < min_ratio:
                    min_ratio = ratio
            if min_ratio >= resolved_threshold:
                current_group.append(p)
            else:
                peak_groups.append(current_group)
                current_group = [p]

    return peak_groups

def get_peak_group_window(mz, intensities, group, window_threshold=0.005):
    min_peak_mz = group[0][0]
    min_peak_idx = dt.fastnearest(mz, min_peak_mz)

    max_peak_mz = group[-1][0]
    max_peak_idx = dt.fastnearest(mz, max_peak_mz)

    #Iterate out from the min and max peaks to find the window.
    #Cut off the window if the relative intensity (relative to the max intensity peak) drops below the threshold
    max_int = max(group, key=lambda x: x[1])[1]

    left_idx = min_peak_idx
    while left_idx > 0 and intensities[left_idx] / max_int >= window_threshold:
        left_idx -= 1

    right_idx = max_peak_idx
    while right_idx < len(intensities) - 1 and intensities[right_idx] / max_int >= window_threshold:
        right_idx += 1

    mz_window = mz[left_idx:right_idx + 1]
    intensities_window = intensities[left_idx:right_idx + 1]

    return mz_window, intensities_window

def model_distribution(x, sigma, *params):
    model = np.zeros_like(x)
    params = np.array(params).flatten()
    num_gaussians = len(params) // 2
    heights = params[:num_gaussians]
    positions = params[num_gaussians:]
    for i in range(len(heights)):
        model += gaussian(x, heights[i], positions[i], sigma)
    return model

def conv_fit(mz, intensities, group):
    sigma_guess = (mz[-1] - mz[0]) / len(group)
    heights = group[:, 1]
    positions = group[:, 0]


    p0 = [sigma_guess]
    for i in range(len(heights)):
        p0.append(heights[i])
    for i in range(len(positions)):
        p0.append(positions[i])


    position_variation = 0.05*sigma_guess
    sigma_variation = 2*sigma_guess

    sigma_min = sigma_guess - sigma_variation
    if sigma_min < 0:
        sigma_min = 0

    lower_bounds = [sigma_min] + [0] * len(group) + [pos - position_variation for pos in positions]
    upper_bounds = [sigma_guess + sigma_variation] + [np.inf] * len(group) + [pos + position_variation for pos in positions]

    try:
        opt, _ = curve_fit(model_distribution, mz, intensities, p0=p0, bounds=(lower_bounds, upper_bounds), ftol=1e-4, xtol=1e-4)
    except:
        print("Curve fitting failed.")
        return 0,[],[]

    optimized_sigma = opt[0]
    optimized_heights = opt[1:len(group) + 1]
    optimized_positions = opt[len(group) + 1:]


    modelled = model_distribution(mz, optimized_sigma, optimized_heights, optimized_positions)

    if len(group) > 10:
        plt.plot(mz, intensities, label="Original Data")
        plt.plot(mz, modelled, label="Fitted Model", linestyle='--')
        plt.scatter(group[:, 0], group[:, 1], color='red', label="Peaks")
        plt.show()

    return optimized_sigma, optimized_heights, optimized_positions

def get_centroids(positions, heights):
    centroids = []
    for i in range(len(positions)):
        centroids.append([positions[i],heights[i]])
    return centroids

spec = np.loadtxt("Z:\\Group Share\\JGP\\DataForJoe\\dat.txt")
# spec = np.loadtxt("Z:\\Group Share\\JGP\\DataForJoe\\Volker\\Protein+cov-binder_#212#919.txt")
# spec[:,1] = np.convolve(spec[:,1], np.ones(10)/10, mode='same')

# spec = ud.datachop(spec, 1350,1355)

mz = spec[:, 0]
intensities = spec[:, 1]

peaks = dt.fastpeakdetect(spec, window=3)
print("Peaks:", len(peaks))

centroids = []

max_groups = 100

groups_processed = 0

groups = get_peak_groups(mz, intensities, peaks,resolved_threshold=0.05)
print("Total groups:", len(groups))
for group in groups:
    groups_processed += 1
    if groups_processed % 10 == 0:
        print("Processed", groups_processed, "groups.")
    if len(group) > max_groups:
        print("Failed to process group with", len(group), "peaks.")
        continue
    mz_window, intensity_window = get_peak_group_window(mz, intensities, np.array(group))
    sigma, heights, positions = conv_fit(mz_window, intensity_window, np.array(group))
    curr_centroids = get_centroids(positions, heights)
    centroids.extend(curr_centroids)

plt.plot(mz, intensities, label="Original Data")
for i in range(len(centroids)):
    plt.plot([centroids[i][0], centroids[i][0]], [0, centroids[i][1]], color='red')
plt.show()




