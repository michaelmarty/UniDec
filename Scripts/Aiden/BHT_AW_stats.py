import warnings
import myqueue as mq
import pandas as pd
import numpy as np
import os
from scipy import stats as st
from scipy.stats import dunnett

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
            temp = pd.read_csv(file_path, delimiter=" ", header=0, usecols=["Mass", "Area", "MassFWHM", "Height", "MassCentroid"])
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
    Does comparisons of each ion statistic using Dunnett's Method for multiple comparisons.
    The key 21 is treated as the control group, and all others as treatments.

    :param array_dict: list of dicts; List where each dict represents one experiment with keys
                       being injection counts and values being dicts with arrays of measurements
                       e.g., [{1: {"Mass": [val1, val2, ...], "Area": [...], "MassFWHM": [...]},
                               3: {"Mass": [...], "Area": [...], "MassFWHM": [...]},
                               ...
                               21: {"Mass": [...], "Area": [...], "MassFWHM": [...]}},
                              {...}, {...}]  # Multiple experiments (e.g., triplicate)
    :param fname: str; name of the dataset
    :return: Writes a file containing the Dunnett's test statistics and p-values
    """

    def sigma(a):
        """Convert FWHM to standard deviation"""
        return a / (2 * np.sqrt(2 * np.log(2)))

    def p_value_check(results_df, alpha=0.05):
        """
        Check if all comparisons show convergence (non-significant differences)

        :param results_df: DataFrame with p-values
        :param alpha: significance level
        :return: tuple (bool, str) - convergence status and message
        """
        mass_converged = all(results_df["p-value (Mass)"] > alpha)
        fwhm_converged = all(results_df["p-value (FWHM)"] > alpha)
        area_converged = all(results_df["p-value (Area)"] > alpha)
        height_converged = all(results_df["p-value (Height)"] > alpha)
        centroid_converged = all(results_df["p-value (Centroid Mass)"] > alpha)

        if mass_converged and fwhm_converged and area_converged and height_converged and centroid_converged:
            return True, "Full Convergence Found"
        else:
            failed = []
            if not mass_converged:
                failed.append("Mass")
            if not fwhm_converged:
                failed.append("FWHM")
            if not area_converged:
                failed.append("Area")
            if not height_converged:
                failed.append("Height")
            if not centroid_converged:
                failed.append("Centroid Mass")
            return False, f"Failed by {', '.join(failed)} Comparison(s)"

    def run_anova(groups):
        """
        One-way ANOVA across all injection-level groups for a single species/ion.
        Serves as the test that gates the Dunnett's post-hoc.

        Warnings emitted by scipy (e.g. precision RuntimeWarnings) are captured
        but not treated as failures: the returned F/p are kept as long as p is
        finite and in [0, 1]. Only genuine computational failures (a raised
        exception, or an out-of-range/NaN p-value) fall back to (0.0, 1.0).

        :param groups: list of 1D arrays, one per injection level (control + treatments)
        :return: tuple (F-statistic, p-value).
        """
        pooled = np.concatenate(groups)
        if np.var(pooled, ddof=1) < 1e-10:
            return 0.0, 1.0
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            try:
                result = st.f_oneway(*groups)
            except (ValueError, FloatingPointError, ZeroDivisionError):
                return 0.0, 1.0
        pv = float(result.pvalue)
        stat = float(result.statistic)
        if not np.isfinite(pv) or pv < 0.0 or pv > 1.0:
            return 0.0, 1.0
        return stat, pv

    def run_dunnett_multiarm(treatments, control, ion_idx,
                             stat_mat, pval_mat, precision_list, warning_list):
        """
        Run a SINGLE multi-arm Dunnett's test for one ion: every treatment
        injection level is compared against the shared control in one call.

        :param treatments: list of 1D arrays, ordered to match treatment_keys
        :param control: 1D array (the 21-injection control)
        :param ion_idx: row of the result matrices to fill
        :param stat_mat, pval_mat: 2D arrays, shape (num_ions, num_treatments).
                                   Row ion_idx is filled with the per-treatment result.
        :param precision_list: ions where the call failed / returned unusable values
        :param warning_list: (ion_idx, messages) where a warning was emitted but the
                             result was still usable and is retained.

        Fallback behavior:
          - A raised exception sets the whole row to stat=0.0 / p=1.0 and records
            the ion in precision_list.
          - Individual non-finite or out-of-range p-values are clamped to
            stat=0.0 / p=1.0 and the ion is recorded in precision_list.
          - Otherwise results are kept; any captured warning is logged.
        """
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                result = dunnett(*treatments, control=control, alternative='two-sided')
            except (ValueError, FloatingPointError, ZeroDivisionError):
                stat_mat[ion_idx, :] = 0.0
                pval_mat[ion_idx, :] = 1.0
                precision_list.append(ion_idx)
                return
        stats = np.atleast_1d(np.asarray(result.statistic, dtype=float)).copy()
        pvals = np.atleast_1d(np.asarray(result.pvalue, dtype=float)).copy()
        # Clamp any individually bad p-values without discarding the good ones
        bad = ~np.isfinite(pvals) | (pvals < 0.0) | (pvals > 1.0)
        if bad.any():
            stats[bad] = 0.0
            pvals[bad] = 1.0
            precision_list.append(ion_idx)
        stat_mat[ion_idx, :] = stats
        pval_mat[ion_idx, :] = pvals
        if caught:
            msgs = sorted({f"{w.category.__name__}: {w.message}" for w in caught})
            warning_list.append((ion_idx, msgs))

    # Setup file paths
    part_path = "Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Analysis\\"
    rest_path = fname + "results.txt"
    final_path = os.path.join(part_path, rest_path)

    # Get all unique injection counts from the first experiment dict
    keys_list = sorted(list(array_dict[0].keys()))
    control_key = 21  # Control group
    treatment_keys = [k for k in keys_list if k != control_key]  # All others are treatments

    num_experiments = len(array_dict)

    # Extract control group data (21 injections) from all experiments
    control_mass_list = []
    control_area_list = []
    control_fwhm_list = []
    control_height_list = []
    control_centroid_list = []
    for exp in array_dict:
        control_mass_list.append(np.array(exp[control_key]["Mass"]))
        control_area_list.append(np.array(exp[control_key]["Area"]))
        control_fwhm_list.append(np.array(exp[control_key]["MassFWHM"]))
        control_height_list.append(np.array(exp[control_key]["Height"]))
        control_centroid_list.append(np.array(exp[control_key]["MassCentroid"]))

    # Stack arrays: shape will be (num_experiments, num_ions)
    control_mass = np.vstack(control_mass_list)
    control_area = np.vstack(control_area_list)
    control_fwhm = np.vstack(control_fwhm_list)
    control_height = np.vstack(control_height_list)
    control_centroid = np.vstack(control_centroid_list)

    num_ions = control_mass.shape[1]

    # Number of treatment groups for Dunnett's adjustment
    k_treatments = len(treatment_keys)
    alpha = 0.05

    # Pre-extract all treatment data by injection count (needed for the
    # ANOVA, which spans every injection level, and reused in the post-hoc loop)
    treat_mass_by_key = {}
    treat_area_by_key = {}
    treat_fwhm_by_key = {}
    treat_height_by_key = {}
    treat_centroid_by_key = {}
    for treat_key in treatment_keys:
        tm, ta, tf, th, tc = [], [], [], [], []
        for exp in array_dict:
            tm.append(np.array(exp[treat_key]["Mass"]))
            ta.append(np.array(exp[treat_key]["Area"]))
            tf.append(np.array(exp[treat_key]["MassFWHM"]))
            th.append(np.array(exp[treat_key]["Height"]))
            tc.append(np.array(exp[treat_key]["MassCentroid"]))
        treat_mass_by_key[treat_key] = np.vstack(tm)
        treat_area_by_key[treat_key] = np.vstack(ta)
        treat_fwhm_by_key[treat_key] = np.vstack(tf)
        treat_height_by_key[treat_key] = np.vstack(th)
        treat_centroid_by_key[treat_key] = np.vstack(tc)

    #  one-way ANOVA per species/ion, across control + all treatment levels.
    # A significant result is the precondition for running Dunnett's on that ion.
    mass_anova_stat = np.zeros(num_ions)
    mass_anova_pvalues = np.zeros(num_ions)
    area_anova_stat = np.zeros(num_ions)
    area_anova_pvalues = np.zeros(num_ions)
    height_anova_stat = np.zeros(num_ions)
    height_anova_pvalues = np.zeros(num_ions)
    centroid_anova_stat = np.zeros(num_ions)
    centroid_anova_pvalues = np.zeros(num_ions)
    for ion_idx in range(num_ions):
        mass_groups = [control_mass[:, ion_idx]] + \
                      [treat_mass_by_key[k][:, ion_idx] for k in treatment_keys]
        area_groups = [control_area[:, ion_idx]] + \
                      [treat_area_by_key[k][:, ion_idx] for k in treatment_keys]
        height_groups = [control_height[:, ion_idx]] + \
                [treat_height_by_key[k][:, ion_idx] for k in treatment_keys]
        centroid_groups = [control_centroid[:, ion_idx]] + \
            [treat_centroid_by_key[k][:, ion_idx] for k in treatment_keys]
        mass_anova_stat[ion_idx], mass_anova_pvalues[ion_idx] = run_anova(mass_groups)
        area_anova_stat[ion_idx], area_anova_pvalues[ion_idx] = run_anova(area_groups)
        height_anova_stat[ion_idx], height_anova_pvalues[ion_idx] = run_anova(height_groups)
        centroid_anova_stat[ion_idx], centroid_anova_pvalues[ion_idx] = run_anova(centroid_groups)
    with open(final_path, "w", encoding="utf-8") as f:
        f.write(f"Dunnett's Method for Multiple Comparisons\n")
        f.write(f"Control Group: {control_key} injections ({num_experiments} experiments)\n")
        f.write(f"Treatment Groups: {treatment_keys}\n")
        f.write(f"Number of Experiments: {num_experiments}\n")
        f.write(f"Number of Ions/Species: {num_ions}\n")
        f.write(f"Significance Level: α = {alpha}\n")
        f.write(f"Number of Comparisons per Ion: {k_treatments}\n")

        # Degrees of freedom for each test (assuming all replicates present at
        # every injection level for every species)
        num_groups = k_treatments + 1                 # treatments + control
        n_total = num_experiments * num_groups        # observations per species
        df_between = num_groups - 1
        df_within = n_total - num_groups
        df_dunnett = n_total - num_groups
        df_fwhm_n = num_experiments - 1
        df_fwhm_d = num_experiments - 1
        f.write("\nDegrees of freedom (assuming all replicates present at every level):\n")
        f.write(f"  One-way ANOVA (per species): between = {df_between}, "
                f"within = {df_within}  ->  F({df_between}, {df_within})\n")
        f.write(f"  Dunnett's post-hoc (per species, multi-arm): pooled error df = {df_dunnett} "
                f"(= N - groups); comparisons per family = {k_treatments}\n")
        f.write(f"  FWHM two-tailed F-test (per pair): F({df_fwhm_n}, {df_fwhm_d})\n")
        f.write("=" * 80 + "\n\n")

        #  ANOVA summary (per species/ion) — gates the Dunnett's post-hoc
        anova_df = pd.DataFrame({
            "F (Mass)": mass_anova_stat,
            "p (Mass)": mass_anova_pvalues,
            "Sig (Mass)": mass_anova_pvalues < alpha,
            "F (Mass Centroid)": centroid_anova_stat,
            "p (Mass Centroid)": centroid_anova_pvalues,
            "Sig (Mass Centroid)": centroid_anova_pvalues < alpha,
            "F (Area)": area_anova_stat,
            "p (Area)": area_anova_pvalues,
            "Sig (Area)": area_anova_pvalues < alpha,
            "F (Height)": height_anova_stat,
            "p (Height)": height_anova_pvalues,
            "Sig (Height)": height_anova_pvalues < alpha,
        })
        f.write(" One-Way ANOVA across all injection levels (per species)\n")
        f.write("Dunnett's post-hoc is performed only for species with a significant ANOVA.\n")
        f.write("(FWHM reproducibility is assessed separately by F-test and is not gated by ANOVA.)\n\n")
        f.write(anova_df.to_string())
        f.write("\n\n" + "=" * 80 + "\n")

        # ------------------------------------------------------------------
        # Post-hoc computation phase: ONE multi-arm Dunnett's call per ion
        # ------------------------------------------------------------------
        # For each species/ion, all treatment injection levels are compared to
        # the shared 21-injection control in a single Dunnett's call. Results
        # come back as vectors aligned with treatment_keys; we store them in
        # (num_ions x num_treatments) matrices and slice per injection when
        # writing. A gate/precision/warning applies to every injection level of
        # that ion, so it is recorded once (per ion) rather than per comparison.
        n_treat_groups = len(treatment_keys)
        shape = (num_ions, n_treat_groups)

        mass_stat = np.zeros(shape); mass_pval = np.ones(shape)
        centroid_stat = np.zeros(shape); centroid_pval = np.ones(shape)
        area_stat = np.zeros(shape); area_pval = np.ones(shape)
        height_stat = np.zeros(shape); height_pval = np.ones(shape)
        fwhm_stat = np.ones(shape); fwhm_pval = np.ones(shape)

        anova_gated = {"Mass": [], "MassCentroid": [], "Area": [], "Height": []}
        precision_warnings = {"Mass": [], "MassCentroid": [], "Area": [], "Height": [], "FWHM": []}
        warning_records = {"Mass": [], "MassCentroid": [], "Area": [], "Height": []}

        # (metric name, ANOVA p-values, treatment data by key, control array,
        #  stat matrix, p-value matrix) for the four Dunnett-tested metrics
        dunnett_specs = [
            ("Mass", mass_anova_pvalues, treat_mass_by_key, control_mass, mass_stat, mass_pval),
            ("MassCentroid", centroid_anova_pvalues, treat_centroid_by_key, control_centroid, centroid_stat, centroid_pval),
            ("Area", area_anova_pvalues, treat_area_by_key, control_area, area_stat, area_pval),
            ("Height", height_anova_pvalues, treat_height_by_key, control_height, height_stat, height_pval),
        ]

        for ion_idx in range(num_ions):
            for name, anova_p, treat_by_key, control_arr, stat_mat, pval_mat in dunnett_specs:
                if anova_p[ion_idx] >= alpha:
                    # ANOVA non-significant -> no post-hoc; row stays 0/1
                    anova_gated[name].append(ion_idx)
                else:
                    treatments = [treat_by_key[k][:, ion_idx] for k in treatment_keys]
                    run_dunnett_multiarm(treatments, control_arr[:, ion_idx], ion_idx,
                                         stat_mat, pval_mat,
                                         precision_warnings[name], warning_records[name])

            # FWHM: variance F-test, still done per treatment/control pair
            control_sigma = sigma(control_fwhm[:, ion_idx])
            control_var = np.var(control_sigma, ddof=1)
            for j, k in enumerate(treatment_keys):
                treat_sigma = sigma(treat_fwhm_by_key[k][:, ion_idx])
                treat_var = np.var(treat_sigma, ddof=1)
                if treat_var < 1e-10 or control_var < 1e-10:
                    fwhm_stat[ion_idx, j] = 1.0
                    fwhm_pval[ion_idx, j] = 1.0
                    if ion_idx not in precision_warnings["FWHM"]:
                        precision_warnings["FWHM"].append(ion_idx)
                else:
                    dfn = max(len(treat_sigma) - 1, 1)
                    dfd = len(control_sigma) - 1
                    f_stat = treat_var / control_var
                    p_lower = st.f.cdf(f_stat, dfn, dfd)
                    p_upper = 1 - p_lower
                    fwhm_stat[ion_idx, j] = f_stat
                    fwhm_pval[ion_idx, j] = 2 * min(p_lower, p_upper)

        # ------------------------------------------------------------------
        # Diagnostics: reported once (each applies to all injection levels)
        # ------------------------------------------------------------------
        f.write("\nPost-hoc diagnostics (apply to all injection-level comparisons of that species)\n")
        f.write("-" * 80 + "\n")
        wrote_diag = False
        for metric, ion_list in anova_gated.items():
            if ion_list:
                f.write(f"NOTE: {metric} - Ions {ion_list} skipped for post-hoc "
                        f"(One-way ANOVA not significant, p >= {alpha}); p-value set to 1.0.\n")
                wrote_diag = True
        for metric, ion_list in precision_warnings.items():
            if ion_list:
                f.write(f"NOTE: {metric} - Ions {ion_list} had near-zero variance / "
                        f"unusable statistic; p-value set to 1.0.\n")
                wrote_diag = True
        for metric, records in warning_records.items():
            if records:
                f.write(f"NOTE: {metric} - Dunnett's returned a usable result but emitted "
                        f"warnings for these ions (result retained):\n")
                for ion_idx, msgs in records:
                    for msg in msgs:
                        f.write(f"      Ion {ion_idx}: {msg}\n")
                wrote_diag = True
        if not wrote_diag:
            f.write("None.\n")
        f.write("=" * 80 + "\n")

        # ------------------------------------------------------------------
        # Per-injection output tables (column j of each matrix = one level)
        # ------------------------------------------------------------------
        all_comparison_results = {}
        for j, treat_key in enumerate(treatment_keys):
            comparison_df = pd.DataFrame({
                "t-score (Mass)": mass_stat[:, j],
                "p-value (Mass)": mass_pval[:, j],
                "t-score (Centroid Mass)": centroid_stat[:, j],
                "p-value (Centroid Mass)": centroid_pval[:, j],
                "t-score (Area)": area_stat[:, j],
                "p-value (Area)": area_pval[:, j],
                "t-score (Height)": height_stat[:, j],
                "p-value (Height)": height_pval[:, j],
                "F-statistic": fwhm_stat[:, j],
                "p-value (FWHM)": fwhm_pval[:, j]
            })

            all_comparison_results[treat_key] = comparison_df

            f.write(f"\n{'=' * 80}\n")
            f.write(f"Comparison: {treat_key} injections vs {control_key} injections (Control)\n")
            f.write(f"{'=' * 80}\n\n")
            f.write(comparison_df.to_string())
            f.write("\n\n")

            converged, message = p_value_check(comparison_df, alpha)
            f.write(f"Convergence Status: {message}\n")
            if converged:
                f.write(f"✓ Full convergence achieved at {treat_key} injections!\n")
                f.write(f"  All p-values > {alpha}, indicating no significant difference from control.\n")
            else:
                f.write(f"✗ {message}\n")
                f.write(f"  Some p-values ≤ {alpha}, indicating significant differences remain.\n")

        # Summary section
        f.write(f"\n\n{'=' * 80}\n")
        f.write(f"SUMMARY: All Comparisons vs Control ({control_key} injections)\n")
        f.write(f"{'=' * 80}\n\n")

        for treat_key in treatment_keys:
            converged, message = p_value_check(all_comparison_results[treat_key], alpha)
            status = "✓ CONVERGED" if converged else "✗ NOT CONVERGED"
            f.write(f"{treat_key} injections: {status} - {message}\n")

        f.write(f"\n{'=' * 80}\n")
        f.write("Analysis complete.\n")

    print(f"Results written to: {final_path}")
    return all_comparison_results

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


def select_indices(replicates: list, select: list, length: int):
    """
    Trims each injection DataFrame in every replicate to a fixed window of rows
    starting at the reference-mass index found by tolerant_index.

    Mutates replicates in-place so that each replicates[i][key] contains only
    the ``length`` rows beginning at select[i][key].

    :param replicates: list of array dicts to trim (modified in-place)
    :param select: list of dicts {injection_key: start_row_index}, one per replicate,
                   as returned by tolerant_index
    :param length: number of rows to keep after (and including) the start index
    """
    key_list = sorted(list(replicates[0].keys()))
    for key in key_list:
        for i in range(len(replicates)):
            start = select[i][key]
            end = start + length
            replicates[i][key] = replicates[i][key].iloc[start:end].reset_index(drop=True)

def normalize_height(replicates: list):
    for rep in replicates:
        keys_list = sorted(list(rep.keys()))
        for i in keys_list:
            rep[i]["Height"] = rep[i]["Height"] / rep[i]["Height"].sum()
            rep[i]["Area"] = rep[i]["Area"] / rep[i]["Area"].sum()


def main():
    pmr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 1")
    pmr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 2")
    pmr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Protein Mix Rep 3")
    replicates = [pmr1, pmr2, pmr3]
    normalize_height(replicates)
    comparison(replicates, "protein_mix_")

    del pmr1, pmr2, pmr3, replicates

    # Since there is an increase in inherent heterogeneity- will add a tolerance threshold for the 5 most abundant glycoforms and apply statistical testing to those.
    # For starters, will comb through until the first glycoform is found within 50 Da and the following 4 rows will be included for comparisons until Rabbit mAbs are used
    wmabr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 1", und=9)
    wmabr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 2", und=9)
    wmabr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Waters mAb Rep 3", und=9)
    select_r1 = tolerant_index(wmabr1, 148220, 20)
    select_r2 = tolerant_index(wmabr2, 148220, 20)
    select_r3 = tolerant_index(wmabr3, 148220, 20)
    select = [select_r1, select_r2, select_r3]
    wmab_replicates = [wmabr1, wmabr2, wmabr3]
    select_indices(wmab_replicates, select, length=5)
    normalize_height(wmab_replicates)
    comparison(wmab_replicates, fname="waters_mab_")

    del wmabr3, wmabr2, wmabr1, wmab_replicates

    ahmab1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Human IgG Rep 1", und=10)
    ahmab2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Human IgG Rep 2", und=10)
    ahmab3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Human IgG Rep 3", und=10)
    select_r1 = tolerant_index(ahmab1, ref_mass=147997, threshold=25)
    select_r2 = tolerant_index(ahmab2, ref_mass=147997, threshold=25)
    select_r3 = tolerant_index(ahmab3, ref_mass=147997, threshold=25)
    select = [select_r1, select_r2, select_r3]
    ahmab_replicates = [ahmab1, ahmab2, ahmab3]
    select_indices(ahmab_replicates, select, length=4)
    normalize_height(ahmab_replicates)
    comparison(ahmab_replicates, fname="amgen_human_mab_")

    del ahmab1, ahmab2, ahmab3, ahmab_replicates

    armabr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Rabbit IgG Rep 1", 10)
    armabr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Rabbit IgG Rep 2", 10)
    armabr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen Rabbit IgG Rep 3", 10)
    select_r1 = tolerant_index(armabr1, ref_mass=143584, threshold=25)
    select_r2 = tolerant_index(armabr2, ref_mass=143584, threshold=25)
    select_r3 = tolerant_index(armabr3, ref_mass=143584, threshold=25)
    select = [select_r1, select_r2, select_r3]
    armab_replicates = [armabr1, armabr2, armabr3]
    select_indices(armab_replicates, select, length=9)
    normalize_height(armab_replicates)
    comparison(armab_replicates, fname="amgen_rabbit_mab_")

    del armabr1, armabr2, armabr3, armab_replicates

    aigmlf_r1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen IgM-like Fusion Rep 1", und=8)
    aigmlf_r2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen IgM-like Fusion Rep 2", und=8)
    aigmlf_r3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\Amgen IgM-like Fusion Rep 3", und=8)
    aigmlf_replicates = [aigmlf_r1, aigmlf_r2, aigmlf_r3]
    # normalize_height(aigmlf_replicates)
    comparison(aigmlf_replicates, fname="amgen_IGM-like_fusion_")

    del aigmlf_replicates, aigmlf_r1, aigmlf_r2, aigmlf_r3

    # High Throughput Section
    wmabr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Waters mAb Rep 1", und=8)
    wmabr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Waters mAb Rep 2", und=8)
    wmabr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Waters mAb Rep 3", und=8)
    select_r1 = tolerant_index(wmabr1, 148220, 40)
    select_r2 = tolerant_index(wmabr2, 148220, 40)
    select_r3 = tolerant_index(wmabr3, 148220, 40)
    select = [select_r1, select_r2, select_r3]
    wmab_replicates = [wmabr1, wmabr2, wmabr3]
    select_indices(wmab_replicates, select, length=5)
    normalize_height(wmab_replicates)
    comparison(wmab_replicates, fname="waters_mab_high_throughput_")

    del wmabr3, wmabr2, wmabr1, wmab_replicates

    ahmab1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Amgen Human IgG Rep 1", und=9)
    ahmab2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Amgen Human IgG Rep 2", und=9)
    ahmab3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Amgen Human IgG Rep 3", und=9)
    select_r1 = tolerant_index(ahmab1, ref_mass=147997, threshold=40)
    select_r2 = tolerant_index(ahmab2, ref_mass=147997, threshold=40)
    select_r3 = tolerant_index(ahmab3, ref_mass=147997, threshold=40)
    select = [select_r1, select_r2, select_r3]
    ahmab_replicates = [ahmab1, ahmab2, ahmab3]
    select_indices(ahmab_replicates, select, length=4)
    normalize_height(ahmab_replicates)
    comparison(ahmab_replicates, fname="amgen_human_mab_high_throughput_")

    del  ahmab1, ahmab2, ahmab3,ahmab_replicates

    armabr1 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Amgen Rabbit IgG Rep 1", 9)
    armabr2 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Amgen Rabbit IgG Rep 2", 9)
    armabr3 = folder_handler("Z:\\Group Share\\BHT\\Temp\\Paper\\Stats\\High Throughput\\Amgen Rabbit IgG Rep 3", 9)
    select_r1 = tolerant_index(armabr1, ref_mass=143584, threshold=25)
    select_r2 = tolerant_index(armabr2, ref_mass=143584, threshold=25)
    select_r3 = tolerant_index(armabr3, ref_mass=143584, threshold=25)
    select = [select_r1, select_r2, select_r3]
    armab_replicates = [armabr1, armabr2, armabr3]
    select_indices(armab_replicates, select, length=9)
    normalize_height(armab_replicates)
    comparison(armab_replicates, fname="amgen_rabbit_mab_high_throughput_")

    del armabr1, armabr2, armabr3, armab_replicates


if __name__ == "__main__":
    main()
