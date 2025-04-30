import pandas as pd
from pyopenms import *
import os
import fnmatch
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from scipy.stats import zscore
import unidec.IsoDec.datatools as dt
from Scripts.MTM.Lipidomics.mass_defect_analysis import kendrick_analysis
from Scripts.LipiDec.AddFeaturestoDF import get_overall_tails, get_class
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.preprocessing import PolynomialFeatures

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

# def calibrate_library_rts():
#     existing_library = r"Z:\Group Share\JGP\Lipidomics_MassIVE_folder\dda_lipid_library.csv"
#     library_df = pd.read_csv(existing_library, delimiter=",")
#     masses = library_df["m/z"]
#     retention_times = library_df["Retention Time"]
#     names = library_df["Lipid Name"]
#     new_rts = [[] for i in range(len(retention_times))]
#
#     dda_file_directory = "Z:\Group Share\JGP\Lipidomics_MassIVE_folder\Raw\DDA\calibrated_mzML"
#     os.chdir(dda_file_directory)
#
#     dda_files = match_files(dda_file_directory, "*.mzML")
#
#     min_peak_intensity = 1000
#     ppm_tol = 10
#     rt_width = 2
#
#     for file in dda_files:
#         exp = MSExperiment()
#         MzMLFile().load(file, exp)
#         print("Creating indexed file")
#         ind_file = xic_xtractor.IndexedFile(exp)
#
#         for i in range(len(masses)):
#             mass_tol = (masses[i] / 1e6) * ppm_tol
#             min_rt = retention_times[i] - (rt_width / 2)
#             max_rt = retention_times[i] + (rt_width / 2)
#             xic = ind_file.extract_xic(masses[i], mass_tol, min_rt, max_rt)
#             peak_index = np.argmax(xic[:, 2])
#             peak_intensity = xic[peak_index, 2]
#
#             if peak_intensity > min_peak_intensity:
#                 rt_int = np.transpose(np.array([xic[:, 1] / 60, xic[:, 2]]))
#                 peaks = dt.fastpeakdetect(rt_int, window=10, threshold=0.05)
#                 peaks = [p for p in peaks if p[1] > min_peak_intensity]
#                 new_rt = -1
#                 if len(peaks) > 1:
#                     delta_rts = [abs(p[0] - retention_times[i]) for p in peaks]
#                     sorted_indices = np.argsort(delta_rts)
#                     new_rt = peaks[sorted_indices[0]][0]
#                 elif len(peaks) == 1:
#                     new_rt = peaks[0][0]
#
#                 if new_rt != -1:
#                     new_rts[i].append(new_rt)
#                     # plt.plot(xic[:,1]/60, xic[:,2])
#                     # plt.title(names[i])
#                     # plt.vlines(retention_times[i], 0, np.max(xic[:,2]), color='r', linestyles='dashed')
#                     # plt.vlines(new_rt, 0, np.max(xic[:,2]), color='b', linestyles='dashed')
#                     # plt.show()
#
#     rt_averages = [0.0 for i in range(len(new_rts))]
#     for i in range(len(new_rts)):
#         if (len(new_rts[i]) > 0):
#             rt_averages[i] = float(np.average(new_rts[i]))
#         else:
#             print(names[i], "not found in any file.")
#             rt_averages[i] = retention_times[i]
#
#     library_df["Retention Time"] = rt_averages
#     library_df.to_csv("Z:\\Group Share\\JGP\\Lipidomics_MassIVE_folder\\rtcalibrated_dda_library.csv")

def filter_msdial_results(results_df, to_exclude=None):
    metabolite_names = results_df["Metabolite name"]
    indices_to_remove = []
    for i,name in enumerate(metabolite_names):
        if any(exclude in name for exclude in to_exclude):
            indices_to_remove.append(i)

    results_df = results_df.drop(index=indices_to_remove)
    results_df.reset_index(drop=True, inplace=True)
    return results_df

def clean_metabolite_name(name):
    possible_prefixes = ["no MS2: ", "low score: "]

    for pref in possible_prefixes:
        if pref in name:
            name = name.removeprefix(pref)

    return name

def plot_kmds(results_df):
    mzs = []
    rts = []
    names = []
    classes = []

    names = results_df["Metabolite name"]
    names = [clean_metabolite_name(name) for name in names]
    rts = results_df["Average Rt(min)"]
    mzs = results_df["Average Mz"]
    classes = [get_class(name) for name in names]

    kmds = [0 for i in range(len(mzs))]
    rts = [float(i) for i in rts]

    for i in range(len(mzs)):
        mass, kmass, kmd = kendrick_analysis(float(mzs[i]), 1)
        tname, tl, tu = get_overall_tails(names[i])
        kmds[i] = kmd
        # plt.text(rts[i], kmd-0.001, names[i])
        plt.text(rts[i], kmd-0.001, tl)
        if classes[i] == "CL":
            if tu == 0:
                plt.plot(rts[i], kmd, 'ro')
            elif tu == 1:
                plt.plot(rts[i], kmd, 'rs')
            elif tu >= 2:
                plt.plot(rts[i], kmd, 'r^')

        elif classes[i] == "PG":
            if tu == 0:
                plt.plot(rts[i], kmd, 'bo')
            elif tu == 1:
                plt.plot(rts[i], kmd, 'bs')
            elif tu >= 2:
                plt.plot(rts[i], kmd, 'b^')

        elif classes[i] == "PE":
            if tu == 0:
                plt.plot(rts[i], kmd, 'go')
            elif tu == 1:
                plt.plot(rts[i], kmd, 'gs')
            elif tu >= 2:
                plt.plot(rts[i], kmd, 'g^')


    plt.xlabel("Retention Time")
    plt.ylabel("KMD")
    plt.show()

def filter_by_kmd(results_df):
    names = results_df["Metabolite name"]
    names = [clean_metabolite_name(name) for name in names]
    rts = results_df["Average Rt(min)"]
    mzs = results_df["Average Mz"]
    classes = [get_class(name) for name in names]

    unique_classes = list(set(classes))

    to_remove = []

    for lipid_class in unique_classes:
        unique_tus = []
        for i,name in enumerate(names):
            if classes[i] == lipid_class:
                tname,tl,tu = get_overall_tails(name)
                if tu not in unique_tus:
                    unique_tus.append(tu)

        for tu in unique_tus:
            kmds_within_tu = []
            indices = []
            for i,name in enumerate(names):
                if classes[i] == lipid_class:
                    tname, tl, curr_tu = get_overall_tails(name)
                    if curr_tu == tu:
                        mass, kmass, kmd = kendrick_analysis(float(mzs[i]), 1)
                        kmds_within_tu.append(kmd)
                        indices.append(i)

            z_scores = zscore(kmds_within_tu)
            for j in range(len(z_scores)):
                if abs(z_scores[j]) > 2:
                    to_remove.append(indices[j])

    results_df = results_df.drop(index=to_remove)
    results_df.reset_index(drop=True, inplace=True)

    return results_df

def filter_by_tl_rt(results_df):
    #The premise of this filter is that
    names = results_df["Metabolite name"]
    names = [clean_metabolite_name(name) for name in names]
    rts = results_df["Average Rt(min)"]
    mzs = results_df["Average Mz"]
    classes = [get_class(name) for name in names]

    class_rtrange_dict = {"CL":[13,17],
                          "PG":[3,10],
                          "PE":[3,14]}

    unique_classes = list(set(classes))
    for lipid_class in unique_classes:
        class_rts = []
        class_tls = []
        indices = []
        rt_range = class_rtrange_dict[lipid_class]
        for i,name in enumerate(names):
            if classes[i] != lipid_class:
                continue
            if rt_range[0] < rts[i] < rt_range[1]:
                tname, tl, tu = get_overall_tails(name)
                class_rts.append(rts[i])
                class_tls.append(tl)
                indices.append(i)

        class_rts = np.array(class_rts)
        class_tls = np.array(class_tls)

        x = np.array(class_tls).reshape(-1,1)
        y = np.array(class_rts)

        # Reshape for curve_fit
        x = x.flatten()
        y = y.flatten()

        # Quadratic model
        def quadratic(x, a, b, c):
            return a * x ** 2 + b * x + c

        # Parameters
        max_trials = 10000
        residual_threshold = 2
        min_inliers = 0.5  # fraction of data points
        bounds = ([0, -np.inf, -np.inf], [1, np.inf, np.inf])  # e.g., a >= 0


        best_model = None
        best_inlier_mask = None
        best_num_inliers = 0

        for i in range(max_trials):
            # Random 3-point sample (minimum for quadratic)
            sample_indices = np.random.choice(len(x), 5, replace=False)
            x_sample = x[sample_indices]
            y_sample = y[sample_indices]

            try:
                # Fit model to sample
                params, _ = curve_fit(quadratic, x_sample, y_sample, bounds=bounds)

                # Predict all points
                y_pred = quadratic(x, *params)

                # Compute residuals
                residuals = np.abs(y - y_pred)

                # Identify inliers
                inlier_mask = residuals < residual_threshold
                num_inliers = np.sum(inlier_mask)

                # Update best model if better
                if num_inliers > best_num_inliers and num_inliers >= min_inliers * len(x):
                    best_model = params
                    best_inlier_mask = inlier_mask
                    best_num_inliers = num_inliers

            except Exception:
                continue  # If fit fails (e.g., bad sample), skip

        # Final model fit on all inliers
        if best_model is not None:
            final_params, _ = curve_fit(quadratic, x[best_inlier_mask], y[best_inlier_mask], bounds=bounds)
            print("Final parameters (a, b, c):", final_params)

            # Plotting
            x_fit = np.linspace(min(x), max(x), 1000)
            y_fit = quadratic(x_fit, *final_params)

            plt.scatter(x, y, label='All Data', alpha=0.5)
            plt.scatter(x[best_inlier_mask], y[best_inlier_mask], color='green', label='Inliers')
            plt.plot(x_fit, y_fit, color='red', label='Quadratic Fit')
            plt.legend()
            plt.title("Custom RANSAC Quadratic Fit with Bounds")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.show()
        else:
            print("RANSAC failed to find a good model.")




if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    ms_dial_results_file = r"Z:\Group Share\JGP\Lipidomics_MassIVE_folder\Raw\DDA\Height_1_2025_04_29_13_55_07.csv"
    os.chdir(os.path.dirname(ms_dial_results_file))

    msdial_results_df = pd.read_csv(ms_dial_results_file)

    to_exclude = ["RIKEN", "Unknown", "no MS2"]
    filtered_search_results = filter_msdial_results(msdial_results_df, to_exclude)
    filtered_search_results = filter_by_kmd(filtered_search_results)
    filtered_search_results = filter_by_tl_rt(filtered_search_results)
    plot_kmds(filtered_search_results)









