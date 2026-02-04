import pandas as pd
import os
import numpy as np
from unidec.tools import center_of_mass, weighted_std

def fix_inhomogeneous_rt(rtdata):
    lengths = [len(rt[0]) for rt in rtdata]
    # If all lengths are the same, return as is
    if len(set(lengths)) == 1:
        return np.array(rtdata)
    # Otherwise, interpolate to min to max range
    min_time = min([min(rt[0]) for rt in rtdata])
    max_time = max([max(rt[0]) for rt in rtdata])
    spacing = rtdata[0][0][1] - rtdata[0][0][0]
    new_times = np.arange(min_time, max_time, spacing)
    new_rtdata = []
    for times, intensities in rtdata:
        new_intensities = np.interp(new_times, times, intensities, left=0, right=0)
        new_rtdata.append((new_times, new_intensities))
    return np.array(new_rtdata)

def get_rt(rtrows, tcol="Times", icol="Intensities", mode="apex"):
    rtdata = []
    for index, row in rtrows.iterrows():
        times = [float(x) for x in row[tcol].split(",")]
        intensities = [float(x) for x in row[icol].split(",")]
        rtdata.append((times, intensities))
    rtdata = fix_inhomogeneous_rt(rtdata)
    times = rtdata[0,0]
    summed_eic = np.sum(rtdata[:,1], axis=0)
    rtsummed = np.transpose([times, summed_eic])

    com_rt = center_of_mass(rtsummed, relative_cutoff=0.1)
    if mode == "apex":
        max_index = np.argmax(rtsummed[:,1])
        max_pos = rtsummed[max_index,0]
        com_rt = (max_pos, com_rt[1])

    return float(com_rt[0]), float(com_rt[1])

def update_rt(df, cdf, min_window=0.5, max_window=5.0, mode="apex", nstds=8):
    # Find all unique molecule names in transition list
    unique_names = df["Molecule Name"].unique()
    outdf = df.copy()
    # Set Explicity Retention Time and Window as floats to avoid dtype issues
    outdf["Explicit Retention Time"] = outdf["Explicit Retention Time"].astype(float)
    outdf["Explicit Retention Time Window"] = outdf["Explicit Retention Time Window"].astype(float)

    for name in unique_names:
        # name = row["Molecule Name"]
        rt_row = cdf[cdf["PeptideModifiedSequence"] == name]
        # adduct = row["Precursor Adduct"][:-1]  # remove + or - at end
        # rt_row = rt_row[rt_row["PrecursorCharge"] == adduct]
        if not rt_row.empty:
            # df.at[index, "Explicit Retention Time"] = rt_row["Retention Time (min)"].values[0]
            # df.at[index, "Explicit Retention Time Window"] = rt_row["RT Window (min)"].values[0]
            rt, window = get_rt(rt_row, mode=mode)
            window = nstds * window
            window = max(min_window, min(window, max_window))
            outdf.loc[outdf["Molecule Name"] == name, "Explicit Retention Time"] = rt
            outdf.loc[outdf["Molecule Name"] == name, "Explicit Retention Time Window"] = window
            print(f"Updated RT for: {name} to {rt:.2f} +/- {window:.2f}")
        else:
            print(f"No RT found for: {name}")
    return outdf

topdir = r"Z:\Group Share\Shivam\Reducing_FP\Final_results"
os.chdir(topdir)

transition_list_path = r"Combined_TransitionList_Clean.csv"
rtfile = r"HEK.tsv"
output_path = r"Combined_TransitionList_Clean_RT.csv"


df = pd.read_csv(transition_list_path)

cdf = pd.read_csv(rtfile, sep="\t")

df = update_rt(df, cdf, mode="apex", min_window=0.5, max_window=5.0, nstds=4)

# Write to output
df.to_csv(output_path, index=False)


