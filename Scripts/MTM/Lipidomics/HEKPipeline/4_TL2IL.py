import pandas as pd
import os

def tl_to_il(df):
    il_df = pd.DataFrame()
    il_df["Compound"] = df["Molecule Name"]
    il_df["m/z"] = df["Precursor Mz"]
    il_df["z"] = tldf["Precursor Charge"]

    rts = df["Explicit Retention Time"].to_numpy()
    rtwindows = df["Explicit Retention Time Window"].to_numpy()
    il_df["t start (min)"] = rts - (rtwindows / 2)
    il_df["t stop (min)"] = rts + (rtwindows / 2)
    # il_df["Retention Time (min)"] = rts
    # il_df["RT Window (min)"] = rtwindows

    il_df["HCD Collision Energy/Energies (%)"] = 30  # Default value, can be adjusted

    # Drop duplicates
    il_df = il_df.drop_duplicates().reset_index(drop=True)
    return il_df

def derplicate_il(df, mztol=0.8):
    good_rows = []
    df = df.sort_values(by=["m/z"]).reset_index(drop=True)
    # Group similar m/z and charge
    groups = []
    used = set()
    for i, row in df.iterrows():
        if i in used:
            continue
        group_indices = [i]
        for j in range(i + 1, len(df)):
            if j in used:
                continue
            other = df.loc[j]
            mz_close = abs(row["m/z"] - other["m/z"]) <= mztol
            same_charge = row["z"] == other["z"]
            if mz_close and same_charge:
                group_indices.append(j)
        used.update(group_indices)
        group = df.loc[group_indices]
        groups.append(group)

    for group in groups:
        if len(group) == 1:
            good_rows.append(group.iloc[0].to_dict())
        else:
            # print(group.to_string())
            tstarts = group["t start (min)"].to_numpy()
            tstops = group["t stop (min)"].to_numpy()
            # Check if there are any overlaps
            overlap = False
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    if not (tstops[i] < tstarts[j] or tstops[j] < tstarts[i]):
                        overlap = True
                        break
                if overlap:
                    break
            if overlap:
                startt = min(tstarts)
                stopt = max(tstops)
                combined_row = {
                    "Compound": " & ".join(sorted(group["Compound"].astype(str).unique())),
                    "m/z": group["m/z"].mean(),
                    "z": group["z"].iloc[0],
                    "t start (min)": startt,
                    "t stop (min)": stopt,
                    "HCD Collision Energy/Energies (%)": group["HCD Collision Energy/Energies (%)"].iloc[0]
                }
                good_rows.append(combined_row)
                print("Overlap found, new combined row:")
                print(pd.Series(combined_row).to_string())
            else:
                for _, g in group.iterrows():
                    good_rows.append(g.to_dict())
                pass

    outdata = pd.DataFrame(good_rows).reset_index(drop=True)
    outdata["Retention Time (min)"] = (outdata["t start (min)"] + outdata["t stop (min)"]) / 2
    outdata["RT Window (min)"] = outdata["t stop (min)"] - outdata["t start (min)"]

    # Drop tstart and tstop columns
    outdata = outdata.drop(columns=["t start (min)", "t stop (min)"])

    return outdata


topdir = r"Z:\Group Share\Shivam\Reducing_FP\Final_results"
os.chdir(topdir)

tlpath = r"Combined_TransitionList_Clean_RT.csv"
outpath = r"Combined_IsolationList.csv"

tldf = pd.read_csv(tlpath)

ildf = tl_to_il(tldf)
print(len(ildf))
ildf = derplicate_il(ildf)
print(len(ildf))

ildf.to_csv(outpath, index=False)

