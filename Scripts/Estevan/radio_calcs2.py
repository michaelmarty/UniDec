import numpy as np
import math

def cpm_to_dpm(cpm, efficiency=0.74):
    """
    Convert counts per minute (CPM) to disintegrations per minute (DPM).
    :param cpm: Counts per minute
    :param efficiency: Detector efficiency (default is 0.74)
    :return: Disintegrations per minute
    """
    if efficiency <= 0 or efficiency > 1:
        raise ValueError("Efficiency must be between 0 and 1.")
    return cpm / efficiency

def dpm_to_ci(dpm):
    """
    Convert disintegrations per minute (DPM) to Curie.
    :param dpm: Disintegrations per minute
    :return: Ci
    """
    # Convert DPM to Ci (1 Ci = 2.22e12 DPM)
    ci = dpm / 2.22e12
    return ci

def cpm_to_ci(cpm, efficiency=0.74):
    """
    Convert counts per minute (CPM) to Curie.
    :param cpm: Counts per minute
    :param efficiency: Detector efficiency (default is 0.74)
    :return: Ci
    """
    dpm = cpm_to_dpm(cpm, efficiency)
    ci = dpm_to_ci(dpm)
    return ci

def ci_to_mol(ci, specific_activity):
    """
    Convert Curie to moles.
    :param ci: Curie
    :param specific_activity: Specific activity in Ci/mmol
    :return: Moles
    """
    if specific_activity <= 0:
        raise ValueError("Specific activity must be positive.")
    mol = ci / specific_activity
    return mol

def mol_to_conc(mol, volume):
    """
    Convert moles to concentration (M).
    :param mol: Moles
    :param volume: Volume in liters
    :return: Concentration in M
    """
    if volume <= 0:
        raise ValueError("Volume must be positive.")
    conc = mol / volume
    return conc

def predict_cpm(specific_activity, concentration, volume):
    """
    Predict total CPM given the specific activity and final concentration of hot ligand
    :param specific_activity: Specific activity in Ci/mmol
    :param concentration: Concentration in nM
    :param volume: Volume in microliters
    """
    if specific_activity <= 0:
        raise ValueError("Specific activity must be positive.")
    if volume <= 0:
        raise ValueError("Volume must be positive.")
    mol = (concentration/1e9)*(volume*1e-6)
    mmol = mol*1000
    ci = specific_activity * mmol
    dpm = ci*2.22e12
    cpm = dpm*0.74
    return cpm

def main():
    hot_sa = 75  # Ci/mol (average of the reported range)
    hot_conc = 0.25  # nM
    cold_conc = 2000  # nM
    total_conc = hot_conc + cold_conc
    final_sa = hot_sa * (hot_conc / total_conc)  # Ci/mol
    print("Calculated Activity and Receptor values using [H3] ACh")
    print('Final specific activity: {:f} Ci/mol'.format(final_sa))

    assay_vol = 200 * 1e-6  # L

    predict_cpm(hot_sa, total_conc, assay_vol)

    background = 66
    total = predict_cpm(final_sa, total_conc, assay_vol*1e6) + background
    ft = 3328/1000 # using the values from DHA which is around 1000 times more active than the acetylcholine ligand
    wash = 160 /1000
    elution = 641 /1000

    total_ci = cpm_to_ci(total - background)  # Ci
    measured_sa = total_ci / (total_conc * 1e-9 * assay_vol)  # Ci/mol
    print("Measured Specific Activity: {:f} Ci/mol".format(measured_sa))

    all = [background, total, ft, wash, elution]

    ci_all = np.array([cpm_to_ci(x) for x in all])

    mol_all = np.array([ci_to_mol(x, measured_sa) for x in ci_all])

    conc_all = np.array([mol_to_conc(x, assay_vol) for x in mol_all]) * 1e6

    print("CI", ci_all)
    print("Moles", mol_all)
    print("Conc", conc_all)

    binding_vol = 90e-6  # L
    netabs = 0.3317 - 0.0667
    epsilon = 233400 # /M/cm
    conc_r = netabs / epsilon
    mol_r = conc_r * binding_vol
    print("Moles of Receptor: ", mol_r * 1e9, "nmol")

    active_receptor = mol_all[-1] / mol_r

    print("Percent Active Receptor: {:f}".format(active_receptor * 100))
    print("\nCalulated Values using [H3] NMS at the same concentration")
    hot_sa = 75e3  # Ci/mol (average of the reported range)
    hot_conc = 0.25  # nM
    cold_conc = 2000  # nM
    total_conc = hot_conc + cold_conc
    final_sa = hot_sa * (hot_conc / total_conc)  # Ci/mol
    print('Final specific activity: {:f} Ci/mol'.format(final_sa))

    assay_vol = 200 * 1e-6  # L

    predict_cpm(hot_sa, total_conc, assay_vol)

    background = 66
    total = predict_cpm(hot_sa, total_conc, assay_vol*1e6) + background
    ft = 3328   # using the values from DHA which is comparable to the NMS ligand
    wash = 160
    elution = 641

    total_ci = cpm_to_ci(total - background)  # Ci
    measured_sa = total_ci / (total_conc * 1e-9 * assay_vol)  # Ci/mol
    print("Measured Specific Activity: {:f} Ci/mol".format(measured_sa))

    all = [background, total, ft, wash, elution]

    ci_all = np.array([cpm_to_ci(x) for x in all])

    mol_all = np.array([ci_to_mol(x, measured_sa) for x in ci_all])

    conc_all = np.array([mol_to_conc(x, assay_vol) for x in mol_all]) * 1e6

    print("Ci", ci_all)
    print("Moles", mol_all)
    print("Conc", conc_all)

    binding_vol = 90e-6  # L
    netabs = 0.3317 - 0.0667
    epsilon = 233400  # /M/cm
    conc_r = netabs / epsilon
    mol_r = conc_r * binding_vol
    print("Moles of Receptor: ", mol_r * 1e9, "nmol")

    active_receptor = mol_all[-1] / mol_r

    print("Percent Active Receptor: {:f}".format(active_receptor * 100))
if __name__ == '__main__':
    main()
