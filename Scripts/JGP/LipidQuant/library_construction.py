import pandas as pd
from pyopenms import *
import os
import fnmatch
import matplotlib.pyplot as plt
import matplotlib as mpl
import unidec.IsoDec.datatools as dt
from Scripts.MTM.Lipidomics.mass_defect_analysis import kendrick_analysis
from Scripts.LipiDec.AddFeaturestoDF import get_overall_tails, get_class
from unidec.UniDecImporter.ImportTools import IndexedFile, IndexedScan
import re

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

def calibrate_library_rts(library_df, dda_mzml_directory):
    masses = library_df["m/z"]
    retention_times = library_df["Retention Time"]
    names = library_df["Lipid name"]
    new_rts = [[] for i in range(len(retention_times))]

    dda_file_directory = "Z:\Group Share\JGP\Lipidomics_MassIVE_folder\Raw\DDA\calibrated_mzML"
    os.chdir(dda_mzml_directory)

    dda_files = match_files(dda_mzml_directory, "*.mzML")

    min_peak_intensity = 1000
    ppm_tol = 10
    rt_width = 2

    for file in dda_files:
        exp = MSExperiment()
        MzMLFile().load(file, exp)
        print("Creating indexed file for", file)
        ind_file = IndexedFile()

        for i, s in enumerate(exp):
            if s.getMSLevel() == 1:
                peaks = np.transpose(np.vstack(s.get_peaks()))
                rt = s.getRT() / 60
                ind_scan = IndexedScan(peaks, rt, i + 1, 0, 1)
                ind_file.indexed_scans.append(ind_scan)

        for i in range(len(masses)):
            mass_tol = (masses[i] / 1e6) * ppm_tol
            min_rt = retention_times[i] - (rt_width / 2)
            max_rt = retention_times[i] + (rt_width / 2)
            xic = ind_file.extract_xic(masses[i], mass_tol, [min_rt, max_rt])
            peak_index = np.argmax(xic[:, 1])
            peak_intensity = xic[peak_index, 1]

            if peak_intensity > min_peak_intensity:
                rt_int = np.transpose(np.array([xic[:, 0], xic[:, 1]]))
                peaks = dt.fastpeakdetect(rt_int, window=10, threshold=0.05)
                peaks = [p for p in peaks if p[1] > min_peak_intensity]
                new_rt = -1
                if len(peaks) > 1:
                    delta_rts = [abs(p[0] - retention_times[i]) for p in peaks]
                    sorted_indices = np.argsort(delta_rts)
                    new_rt = peaks[sorted_indices[0]][0]
                elif len(peaks) == 1:
                    new_rt = peaks[0][0]

                if new_rt != -1:
                    new_rts[i].append(new_rt)

        del ind_file


    rt_averages = [0.0 for i in range(len(new_rts))]
    for i in range(len(new_rts)):
        if (len(new_rts[i]) > 0):
            rt_averages[i] = float(np.average(new_rts[i]))
        else:
            print(names[i], "not found in any file.")
            rt_averages[i] = retention_times[i]

    library_df["Retention Time"] = rt_averages
    return library_df

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

def get_charge_from_adduct(adduct):
    split = adduct.split("]")
    s = split[1]
    if s in ("","0"):
        return 0
    match = re.fullmatch(r"(\d*)([+-])", s)
    if not match:
        raise ValueError(f"Invalid charge:{adduct}")
    num,sign = match.groups()
    num = int(num) if num else 1
    polarity = "POSITIVE" if sign == "+" else "NEGATIVE"
    return num, polarity

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
        # plt.text(rts[i], kmd-0.001, tl)
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
    adducts = results_df["Adduct type"]
    charges = [get_charge_from_adduct(adduct)[0] for adduct in adducts]
    classes = [get_class(name) for name in names]

    unique_classes = list(set(classes))

    to_remove = []

    for lipid_class in unique_classes:
        class_kmds = []
        class_rts = []
        class_tus = []
        indices = []
        for i,name in enumerate(names):
            if classes[i] == lipid_class:
                tname,tl,tu = get_overall_tails(name)
                mass, kmass, kmd = kendrick_analysis(float(mzs[i]), charges[i])
                class_kmds.append(kmd)
                class_rts.append(rts[i])
                class_tus.append(tu)
                indices.append(i)


        accepted_rts = np.array(class_rts)
        accepted_kmds = np.array(class_kmds)
        accepted_tus = np.array(class_tus)
        accepted_indices = np.array(indices)

        rejected_rts = np.array([])
        rejected_kmds = np.array([])
        rejected_tus = np.array([])
        rejected_indices = np.array([])

        accepted = np.column_stack((accepted_rts, accepted_kmds, accepted_tus, accepted_indices))
        rejected = np.column_stack((rejected_rts, rejected_kmds, rejected_tus, rejected_indices))

        fig, ax = plt.subplots()
        scatter_a = ax.scatter(accepted[:, 0], accepted[:, 1], c='b', label="accepted")
        for i in range(len(accepted)):
            ax.text(accepted[i,0], accepted[i,1], str(accepted[i,2]))
        scatter_b = ax.scatter(rejected[:, 0], rejected[:, 1], c='r', label="rejected")
        for i in range(len(rejected)):
            ax.text(rejected[i,0], rejected[i,1], str(rejected[i,2]))
        ax.legend()
        ax.grid(True)
        ax.set_title(lipid_class)
        fig.canvas.draw()

        def update_plot():
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            ax.clear()
            ax.scatter(accepted[:, 0], accepted[:, 1], c='b', label="accepted")
            for i in range(len(accepted)):
                ax.text(accepted[i, 0], accepted[i, 1], str(accepted[i, 2]))
            for i in range(len(rejected)):
                ax.text(rejected[i, 0], rejected[i, 1], str(rejected[i, 2]))
            ax.scatter(rejected[:, 0], rejected[:, 1], c='r', label="rejected")
            ax.legend()
            ax.set_title(lipid_class)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.grid(True)
            fig.canvas.draw()

        def on_key(event):
            if event.key == 'r':
                ax.relim()  # Recalculate limits based on data
                ax.autoscale_view()  # Autoscale the view to the new limits
                fig.canvas.draw()

        def on_click(event):
            nonlocal accepted, rejected
            if not event.inaxes:
                return

            click_point = np.array([event.xdata, event.ydata])
            threshold = 0.5

            accepted_distances = (
                [np.linalg.norm(click_point - [pt[0], pt[1]]) for pt in accepted]
                if len(accepted) > 0 else []
            )
            rejected_distances = (
                [np.linalg.norm(click_point - [pt[0], pt[1]]) for pt in rejected]
                if len(rejected) > 0 else []
            )

            # Initialize defaults
            min_accepted_dist = float('inf')
            min_rejected_dist = float('inf')
            min_accepted = None
            min_rejected = None

            if accepted_distances:
                min_accepted = accepted_distances.index(min(accepted_distances))
                min_accepted_dist = accepted_distances[min_accepted]

            if rejected_distances:
                min_rejected = rejected_distances.index(min(rejected_distances))
                min_rejected_dist = rejected_distances[min_rejected]

            if min_accepted is not None and (
                    min_accepted_dist < min_rejected_dist and min_accepted_dist < threshold
            ):
                pt = accepted[min_accepted]
                rejected = np.vstack([rejected, pt]) if len(rejected) > 0 else np.array([pt])
                accepted = np.delete(accepted, min_accepted, axis=0)
                update_plot()
                return

            elif min_rejected is not None and min_rejected_dist < threshold:
                pt = rejected[min_rejected]
                accepted = np.vstack([accepted, pt]) if len(accepted) > 0 else np.array([pt])
                rejected = np.delete(rejected, min_rejected, axis=0)
                update_plot()
                return

            else:
                print("Click closer to a point to move it to the other group.")

        fig.canvas.mpl_connect('key_press_event', on_key)
        fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        rejected_indices = rejected[:, 3]
        to_remove.extend(rejected_indices.astype(int))


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

    indices_to_remove = []

    unique_classes = list(set(classes))
    for lipid_class in unique_classes:
        class_rts = []
        class_tls = []
        indices = []
        for i,name in enumerate(names):
            if classes[i] != lipid_class:
                continue

            tname, tl, tu = get_overall_tails(name)
            class_rts.append(rts[i])
            class_tls.append(tl)
            indices.append(i)

        class_rts = np.array(class_rts)
        class_tls = np.array(class_tls)

        accepted_rts = np.array(class_rts)
        accepted_tls = np.array(class_tls)
        accepted_indices = np.array(indices)

        rejected_rts = np.array([])
        rejected_tls = np.array([])
        rejected_indices = np.array([])

        accepted = np.column_stack((accepted_rts, accepted_tls, accepted_indices))
        rejected = np.column_stack((rejected_rts, rejected_tls, rejected_indices))

        fig,ax = plt.subplots()
        scatter_a = ax.scatter(accepted[:,0], accepted[:,1], c='b', label="accepted")
        scatter_b = ax.scatter(rejected[:,0], rejected[:,1], c='r', label="rejected")
        ax.legend()
        ax.set_title(lipid_class)
        fig.canvas.draw()

        def update_plot():
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            ax.clear()
            ax.scatter(accepted[:, 0], accepted[:, 1], c='b', label="accepted")
            ax.scatter(rejected[:, 0], rejected[:, 1], c='r', label="rejected")
            ax.legend()
            ax.set_title(lipid_class)

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            fig.canvas.draw()

        def on_click(event):
            nonlocal accepted, rejected
            if not event.inaxes:
                return

            click_point = np.array([event.xdata, event.ydata])
            threshold = 1

            accepted_distances = (
                [np.linalg.norm(click_point - [pt[0], pt[1]]) for pt in accepted]
                if len(accepted) > 0 else []
            )
            rejected_distances = (
                [np.linalg.norm(click_point - [pt[0], pt[1]]) for pt in rejected]
                if len(rejected) > 0 else []
            )

            # Initialize defaults
            min_accepted_dist = float('inf')
            min_rejected_dist = float('inf')
            min_accepted = None
            min_rejected = None

            if accepted_distances:
                min_accepted = accepted_distances.index(min(accepted_distances))
                min_accepted_dist = accepted_distances[min_accepted]

            if rejected_distances:
                min_rejected = rejected_distances.index(min(rejected_distances))
                min_rejected_dist = rejected_distances[min_rejected]

            if min_accepted is not None and (
                    min_accepted_dist < min_rejected_dist and min_accepted_dist < threshold
            ):
                pt = accepted[min_accepted]
                rejected = np.vstack([rejected, pt]) if len(rejected) > 0 else np.array([pt])
                accepted = np.delete(accepted, min_accepted, axis=0)
                update_plot()
                return

            elif min_rejected is not None and min_rejected_dist < threshold:
                pt = rejected[min_rejected]
                accepted = np.vstack([accepted, pt]) if len(accepted) > 0 else np.array([pt])
                rejected = np.delete(rejected, min_rejected, axis=0)
                update_plot()
                return

            else:
                print("Click closer to a point to move it to the other group.")

        def on_key(event):
            if event.key == 'r':
                ax.relim()  # Recalculate limits based on data
                ax.autoscale_view()  # Autoscale the view to the new limits
                fig.canvas.draw()

        fig.canvas.mpl_connect('key_press_event', on_key)
        fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        rejected_indices = rejected[:,2]
        indices_to_remove.extend(rejected_indices.astype(int))

    results_df = results_df.drop(index=indices_to_remove)
    results_df.reset_index(drop=True, inplace=True)
    return results_df

def search_results_to_library(results_df):
    names = results_df["Metabolite name"]
    names = [clean_metabolite_name(name) for name in names]
    rts = results_df["Average Rt(min)"]
    adducts = results_df["Adduct type"]
    mzs = results_df["Average Mz"]
    formulas = results_df["Formula"]
    classes = [get_class(name) for name in names]
    standards = [False] * len(names)

    library_df = pd.DataFrame()
    library_df["Lipid name"] = names
    library_df["Lipid class"] = classes
    library_df["Standard"] = standards
    library_df["m/z"] = mzs
    library_df["Retention Time"] = rts
    library_df["Adduct"] = adducts
    library_df["Chemical Formula"] = formulas

    library_df = calibrate_library_rts(library_df, "Z:\\Group Share\\JGP\\Lipidomics_MassIVE_folder\\Raw\\DDA\\calibrated_mzML")

    library_df.to_csv("rt_calibrated_filtered_library.csv", index=False)
    return


if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    ms_dial_results_file = r"Z:\Group Share\JGP\Lipidomics_MassIVE_folder\Raw\DDA\DDA_re-search\Height_11_2025_05_01_11_53_45.csv"
    os.chdir(os.path.dirname(ms_dial_results_file))
    msdial_results_df = pd.read_csv(ms_dial_results_file)


    to_exclude = ["RIKEN", "Unknown", "no MS2"]

    filtered_search_results = filter_msdial_results(msdial_results_df, to_exclude)
    filtered_search_results = filter_by_tl_rt(filtered_search_results)
    filtered_search_results = filter_by_kmd(filtered_search_results)
    search_results_to_library(filtered_search_results)
    filtered_search_results.to_csv("filtered_search_results.csv", index=False)
    plot_kmds(filtered_search_results)









