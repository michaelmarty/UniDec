import myqueue as mq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats as st
# from dunnetts import comparison as cp


def folder_handler(folder_name: str, und=8):
    """
    Function that iterates through a folder of .txt (or .dat) files from UniDecChromCD and returns an array of dictionaries.
    Assumes that the file name has the injection number after 8 sets of underscores... will fix eventually
    :param: folder_name: str; file path of folder as a string
    :param: und: int; number of underscores until you reach the injection number
    :return: array_dict: dict{key: pd.DataFrame}; a dictionary of arrays pooled from all of the files in the folder
    """
    file_queue = mq.Queue()
    counter = 0
    array_dict = {}
    try:
        any(file_queue.enqueue(file) for file in os.listdir(folder_name))
        while not file_queue.is_empty:
            file_path = os.path.join(folder_name, str(file_queue.dequeue()))
            number = file_path.split("_")[und]
            temp = pd.read_csv(file_path, delimiter=" ", header=0, usecols=["Mass", "Area", "MassFWHM"])
            key = int(number)
            array_dict[key] = temp
            del temp
            counter += 1
    except FileNotFoundError:
        print('Please Ensure the Folder Exists')
    except UnicodeDecodeError:
        print('Please make sure every file is in the correct format')
    except UnboundLocalError:
        print("Please use a supported file format")
    except IndexError:
        print("fix me")
    return array_dict


def comparison(array_dict, fname=""):
    """
    Does comparisons of each ion statistic and writes results for each into a file/group of files depending on
    :param array_dict:
    :param var: str;
    :param vers: str; the version of comparison test required. Two options: sliding and standard. Sliding compares local variation whereas standard compares each value to the highest number of injections.
    :param fname: str; name of the dataset
    :return: f: file; writes a file containing the t-scores and p-values with which all three comparison variables converge (WIP)
    """

    def p_value(a):
        if all(a["p-value (Mass)"] > 0.05 / 20) and all(a["p-value (FWHM)"] > 0.05 / 20) and all(
                a["p-value (Area)"] > 0.05 / 20):
            print("Full Convergence Found")
            return True
        else:
            if all(a["p-value (Mass)"] > 0.05 / 19):
                if all(a["p-value (Area)"] > 0.05 / 19):
                    print("Failed by FWHM Comparison")
                else:
                    print("Failed by Area Comparison")
            else:
                print("Failed by Mass Comparison")
            print("More Injections Needed")
            return False

    def sigma(a):
        return a / (2 * np.sqrt(2 * np.log(2)))

    part_path = "Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Analysis\\"
    rest_path = fname + "results.txt"
    final_path = os.path.join(part_path, rest_path)
    keys_list = list(array_dict.keys())
    with open(final_path, "w", encoding="utf-8") as f:
        standard_mass = array_dict[keys_list[-1]]["Mass"]
        standard_error = array_dict[keys_list[-1]]["SEM Mass"] / np.sqrt(keys_list[-1])
        standard_fwhm = sigma(array_dict[keys_list[-1]]["MassFWHM"])
        for i in keys_list:
            if i == 21:
                break
            dfarea = i + i + 1 - 2
            tarea = pd.DataFrame((array_dict[i]["Area"] - array_dict[i + 1]["Area"]) / np.sqrt(
                (array_dict[i]["SEM Area"] ** 2 / i + array_dict[i + 1]["SEM Area"] ** 2 / (i + 1))),
                                 columns=["t-score (Area)"])
            parea = pd.DataFrame(1 - st.t.cdf(np.abs(tarea), dfarea), columns=["p-value (Area)"])
            tmass = pd.DataFrame(
                (array_dict[i]["Mass"] - standard_mass) / np.sqrt(
                    (array_dict[i]["SEM Mass"] ** 2) / i + (standard_error ** 2)), columns=["t-score (Mass)"])
            dft = i + keys_list[-1] - 2
            if i != 1:
                dfn = i - 1
            else:
                dfn = 1
            dfd = keys_list[-1] - 1
            fstat = pd.DataFrame(
                (sigma(array_dict[i]["MassFWHM"]) ** 2 / standard_fwhm ** 2)
            ).rename(columns={"MassFWHM": "F-statistic"})
            pf = pd.DataFrame(st.f.cdf(fstat, dfn, dfd), columns=["p-value (FWHM)"])
            p = pd.DataFrame(1 - st.t.cdf(np.abs(tmass), dft), columns=["p-value (Mass)"]
                             )
            final = tmass.merge(p.merge(tarea.merge(parea.merge(fstat.merge(pf,
                                                                            left_index=True, right_index=True,
                                                                            how="outer"),
                                                                left_index=True, right_index=True, how="outer"),
                                                    left_index=True, right_index=True, how="outer"),
                                        left_index=True, right_index=True, how="outer"),
                                left_index=True, right_index=True, how="outer")
            print(final)
            if p_value(final):
                f.write(f"Full Convergence Found at {i} injections for this species using an"
                        f" independent student's t-test with a Bonferroni corrected alpha of {0.05 / (keys_list[-1] - 1)} compared to {keys_list[-1]} injections"
                        f" and an F-test for the ratio of variance\n")
                f.write(final.to_string())
                return
        f.write(
            f"All comparisons yielded p-values below {0.05 / 20}, more injections may be needed; or use a different reference value for indexing\n"
        )
        f.write(
            final.to_string()
        )


def tolerant_index(array_dict, ref_mass=0, threshold=0):
    lower = ref_mass - threshold
    upper = ref_mass + threshold
    first_form = {}
    key_list = sorted(list(array_dict.keys()))
    for key in key_list:
        for i in range(len(array_dict[key]["Mass"])):
            if lower < array_dict[key]['Mass'][i] < upper:
                first_form[key] = i
                break
    del key_list, upper, lower, ref_mass
    return first_form


def merge_replicates(replicates: list, n=3, select=None, length=0):
    """
    Averages all of the injections for each replicate in a set and returns a single dictionary of averaged arrays with standard deviation
    :param replicates: list; list of array dictionaries to be merged
    :param n: int; number of replicates
    :param select: list; (usually) a list of dictionaries containing the index of a reference mass within a set tolerance
    :param length: number of masses after the selection index to include in the final dataframe (number of peaks)
    :return: array_dict: dict{key: pd.DataFrame}; key is the injection number and pd.DataFrame is the averaged array with standard deviation for each variable
    """
    array_dict = {}
    key_list = sorted(list(replicates[0].keys()))
    if select:
        for key in key_list:
            for i in range(len(select)):
                start = select[i][key]
                end = select[i][key] + length
                replicates[i][key] = replicates[i][key][start:end]
    for key in key_list:
        averaged_vars_frame = pd.DataFrame(
            np.mean([replicates[0][key], replicates[1][key], replicates[2][key]], axis=0),
            columns=replicates[0][1].columns)
        standard_dev_frame = pd.DataFrame(
            np.std([replicates[0][key], replicates[1][key], replicates[2][key]], axis=0, ddof=1) / np.sqrt(3),
            columns=["SEM Mass", "SEM Area", "SEM MassFWHM"])
        avg_frame = averaged_vars_frame.merge(standard_dev_frame, left_index=True, right_index=True)
        array_dict[key] = avg_frame
    return array_dict

def build_single_peak_dict(avg_dict: dict, peak_idx: int):
    """
    From a merged average dictionary, build a new dictionary where each key maps
    to a single-row DataFrame corresponding to peak_idx.
    """
    single = {}
    for key in sorted(avg_dict.keys()):
        df = avg_dict[key]
        if peak_idx < 0 or peak_idx >= len(df):
            raise IndexError(f"Peak index {peak_idx} out of range for injection {key}")
        row = df.iloc[[peak_idx]].copy()
        # ensure SEM column names exist; if not, add zeros
        for c in ["SEM Mass", "SEM Area", "SEM MassFWHM"]:
            if c not in row.columns:
                row[c] = 0.0
        single[key] = row.reset_index(drop=True)
    return single


def converge_plot_single_peak(array_dict, peak_idx, ref_mass=None, **kwargs):
    """
    Plot a single peak (peak_idx) across injections. If ref_mass is provided,
    the single-row datasets are compared to that fixed ref.
    """
    # Build a single-row dict
    single = {}
    for k in array_dict:
        single[k] = array_dict[k].iloc[[peak_idx]].reset_index(drop=True)
    if ref_mass is not None:
        converge_plot(single, ref_masses=[ref_mass], **kwargs)
    else:
        converge_plot(single, **kwargs)

def ppm(ref, actual):
    """
    Calculates the ppm difference between two values
    :param ref: int; reference value
    :param actual: int; current value for comparison
    :return: ppm difference; int
    """
    return ((ref - actual) / ref) * (10 ** 6)


def converge_plot(array_dict: dict, data_name=None, colors=None,
                  x_axis='injections', save_pdf=True, output_filename='convergence_plot.pdf',
                  ref_masses=None, ref_tol=50.0):
    """
    Plot convergence data.

    New args:
      - ref_masses: None or list/tuple of reference masses. Behavior:
          * If len(ref_masses) == nrows -> each row i compares to ref_masses[i].
          * If len(ref_masses) == 1 -> same ref used for all rows.
          * Else -> auto-match each ref to the nearest peak row in the final injection
            and assign ref to that row if within ref_tol. Unmatched rows use the
            original behavior (difference to final injection mass).
      - ref_tol: tolerance (Da) used when auto-matching refs to measured peak rows.
    """
    folder_path = r"Z:\Group Share\BHT\Temp\Paper\Stats\Figures"
    final_path = folder_path + "\\" + data_name + "_" + x_axis.lower() + "_" + output_filename

    keys = list(array_dict.keys())
    if not keys:
        print("No data in array_dict to plot.")
        return
    nrows = len(array_dict[keys[0]]["Mass"])

    # Build a reference mass for each row, or None if no ref for that row
    row_refs = [None] * nrows
    if ref_masses is not None:
        refs = list(ref_masses)
        if len(refs) == nrows:
            row_refs = refs[:]  # direct mapping
        elif len(refs) == 1:
            row_refs = [refs[0]] * nrows
        else:
            # auto-match: for each ref find nearest peak row in the final injection
            final_masses = array_dict[keys[-1]]["Mass"].values
            matched = {}
            for ref in refs:
                diffs = np.abs(final_masses - ref)
                idx = int(np.argmin(diffs))
                if diffs[idx] <= ref_tol:
                    # if multiple refs map to same idx, last one wins; typically this won't happen
                    matched[idx] = ref
            for idx, ref in matched.items():
                if 0 <= idx < nrows:
                    row_refs[idx] = ref

    fig, ax = plt.subplots(nrows, 1, sharey=True, figsize=(16, 10))
    if nrows == 1:
        ax = [ax]

    for i in range(nrows):
        data = {}
        error = {}
        x_values = {}

        for key in keys:
            # If a reference mass is provided for this row, compute absolute diff to that ref.
            if row_refs[i] is not None:
                y = np.abs(row_refs[i] - array_dict[key]["Mass"][i])
            else:
                # fallback to the previous behaviour (difference to final injection mass)
                y = np.abs((array_dict[keys[-1]]["Mass"][i] - array_dict[key]["Mass"][i]))
            # errorbars still use SEM Mass (difference to final SEM Mass)
            yerr = np.abs(array_dict[keys[-1]]["SEM Mass"][i] - array_dict[key]["SEM Mass"][i])

            data[key] = y
            error[key] = yerr

            if x_axis == 'Area':
                x_values[key] = array_dict[key]["Area"][i]

        if x_axis == 'injections':
            x_plot = list(data.keys())
            x_label = 'Number of Injections'
        else:
            x_plot = list(x_values.values())
            x_label = 'Area'

        y_plot = list(data.values())
        yerr_plot = list(error.values())

        color = colors[i] if colors and i < len(colors) else None

        ax[i].scatter(x_plot, y_plot, color=color, zorder=3)
        ax[i].plot(x_plot, y_plot, color=color, alpha=0.7)
        ax[i].errorbar(x=x_plot, y=y_plot, yerr=yerr_plot,
                       fmt='o', ecolor='black', elinewidth=2,
                       capsize=5, capthick=2, color=color)

        if x_axis == "injections":
            ax[i].set_xticks(range(1, len(keys) + 1))

        # Title shows either the reference mass (if assigned) or the final injection mass
        title_mass = row_refs[i] if row_refs[i] is not None else array_dict[keys[-1]]["Mass"][i]
        ax[i].set_title(f'Mass reference: {title_mass:.2f}', fontsize=11)

        ax[i].spines["top"].set_visible(False)
        ax[i].spines["right"].set_visible(False)

    ax[-1].set_xlabel(x_label, fontsize=12)
    ax[nrows // 2].set_ylabel('Absolute Difference (Mass Da)', fontsize=10)
    title = f'Convergence Plot for {data_name}'
    fig.suptitle(title, fontsize=14, fontweight='bold')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # leave room for suptitle

    if save_pdf:
        try:
            plt.savefig(final_path, format='pdf', dpi=300, bbox_inches='tight')
            print(f"Figure saved as {final_path} with 300 DPI")
        except Exception as e:
            print(f"Failed to save figure to {final_path}: {e}")


def main():
    pmr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 1")
    pmr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 2")
    pmr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 3")
    replicates = [pmr1, pmr2, pmr3]
    # cp(replicates, fname="try_1_protein_mix")
    pmavg = merge_replicates(replicates)

    # Fixed reference masses to compare against (your requested targets)
    refs = [12360.08, 16951.48, 29024.63]

    # 1) Plot convergence where rows matched to the fixed refs (auto-match with tolerance)
    converge_plot(pmavg, data_name="Protein Mix (refs)", ref_masses=refs, ref_tol=50.0)
    converge_plot(pmavg, data_name="Protein Mix (refs)", x_axis="Area", ref_masses=refs, ref_tol=50.0)

    # 2) For each reference, find the closest peak in the final injection and run the statistical comparison
    try:
        final_masses = pmavg[max(pmavg.keys())]['Mass'].values
        for ref in refs:
            idx = int(np.argmin(np.abs(final_masses - ref)))
            # Optionally check tolerance before running comparison
            if abs(final_masses[idx] - ref) > 50.0:
                print(
                    f"Reference {ref} not within 50 Da of any peak (closest {final_masses[idx]:.2f}); skipping comparison.")
                continue
            single_peak = build_single_peak_dict(pmavg, idx)
            comparison(single_peak, fname=f"protein_mix_{int(ref)}")
    except Exception as e:
        print(f"Error running single-peak comparisons for protein mix: {e}")

        del pmr1, pmr2, pmr3, pmavg
    # Since there is an increase in inherent heterogeneity- will add a tolerance threshold for the 5 most abundant glycoforms and apply statistical testing to those.
    # For starters, will comb through until the first glycoform is found within 50 Da and the following 4 rows will be included for comparisons until Rabbit mAbs are used
    wmabr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 1", und=9)
    wmabr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 2", und=9)
    wmabr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 3", und=9)
    select_r1 = tolerant_index(wmabr1, 148220, 20)
    select_r2 = tolerant_index(wmabr2, 148220, 20)
    select_r3 = tolerant_index(wmabr3, 148220, 20)
    select = [select_r1, select_r2, select_r3]
    wmabavg = merge_replicates([wmabr1, wmabr2, wmabr3], select=select, length=5)
    converge_plot(wmabavg, "Waters mAb")
    converge_plot(wmabavg, "Waters mAb", x_axis="Area")
    comparison(wmabavg, "waters_mAb_")
    del wmabavg, wmabr3, wmabr2, wmabr1
    ahmab1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Human IgG Rep 1", und=10)
    ahmab2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Human IgG Rep 2", und=10)
    ahmab3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Human IgG Rep 3", und=10)
    select_r1 = tolerant_index(ahmab1, ref_mass=147997, threshold=25)
    select_r2 = tolerant_index(ahmab2, ref_mass=147997, threshold=25)
    select_r3 = tolerant_index(ahmab3, ref_mass=147997, threshold=25)
    select = [select_r1, select_r2, select_r3]
    ahmabavg = merge_replicates([ahmab1, ahmab2, ahmab3], select=select, length=4)
    comparison(ahmabavg, fname="amgen_human_mab_")
    converge_plot(ahmabavg, "Amgen Human mAb")
    converge_plot(ahmabavg, "Amgen Human mAb", x_axis="Area")
    del ahmabavg, ahmab1, ahmab2, ahmab3
    armabr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Rabbit IgG Rep 1", 10)
    armabr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Rabbit IgG Rep 2", 10)
    armabr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Rabbit IgG Rep 3", 10)
    select_r1 = tolerant_index(armabr1, ref_mass=143584, threshold=25)
    select_r2 = tolerant_index(armabr2, ref_mass=143584, threshold=25)
    select_r3 = tolerant_index(armabr3, ref_mass=143584, threshold=25)
    select = [select_r1, select_r2, select_r3]
    armabavg = merge_replicates([armabr1, armabr2, armabr3], select=select, length=9)
    comparison(armabavg, fname="amgen_rabbit_mab_")
    converge_plot(armabavg, "Amgen Rabbit mAb")
    converge_plot(armabavg, "Amgen Rabbit mAb", x_axis="Area")
    del armabavg, armabr1, armabr2, armabr3


if __name__ == "__main__":
    main()
