import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import os
import pickle
from PublicScripts.Lipids.FixPGinMSP import fix_pg_in_msp
import matplotlib as mpl

mpl.use("WxAgg")


def extend_df(df, rtcol="Average Rt(min)", refcol="Reference RT"):
    # Append a line beyond this to force fit to go through these points
    # List of rows to add
    template_row = {'Metabolite name': 'End1', 'Average Rt(min)': 0, 'Reference RT': 0, "S/N average": 1,
                    "Total score": 100}
    # For 0 to min rt, add x=y
    min_rt = df[refcol].min()
    new_rows = []
    for rt in np.arange(0, min_rt, 0.1):
        row = template_row.copy()
        row[rtcol] = rt
        row[refcol] = rt
        new_rows.append(row)

    max_rt = df[refcol].max()
    for rt in np.arange(max_rt, max_rt + 10, 0.1):
        row = template_row.copy()
        row[rtcol] = rt
        row[refcol] = rt
        new_rows.append(row)

    # Convert to DataFrame and concatenate
    df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
    return df


def rt_correction(df, rtcol="Average Rt(min)", refcol="Reference RT", weight_col="Total score",
                  limit=5, smooth=1, window=0.5, plot=True):
    # Create a function to correct the Reference RT to match the Average using a 4th degree polynomial
    fit = np.polyfit(df[refcol], df[rtcol], 2)
    p = np.poly1d(fit)
    # Sort by Reference RT
    df = df.sort_values(by=[refcol], ascending=True)
    x = np.array(df[refcol])
    y = np.array(df[rtcol])
    if weight_col is not None:
        w = np.array(df[weight_col])
    else:
        w = np.ones_like(y)

    minx = np.min(x)
    maxx = np.max(x)
    # Create linear spacing of window min between min and max
    xnew = np.arange(minx, maxx, window)
    # Average each y for each xnew
    ynew = []
    wnew = []
    for i in range(len(xnew)):
        xi = xnew[i]
        yvals = y[(x >= xi - window * 0.5) & (x < xi + window * 0.5)]
        wvals = w[(x >= xi - window * 0.5) & (x < xi + window * 0.5)]
        if len(yvals) > 0:
            yavg = np.average(yvals, weights=wvals)
        else:
            yavg = np.nan
        ynew.append(yavg)
        wnew.append(np.sum(wvals))
    ynew = np.array(ynew)
    # Remove NaN
    mask = ~np.isnan(ynew)
    xnew = xnew[mask]
    ynew = ynew[mask]

    # Smooth xnew and ynew using a gaussian filter with sigma of 2
    ynew = gaussian_filter1d(ynew, sigma=smooth)

    # Create an interpolation of xnew and ynew
    p = interp1d(xnew, ynew, kind='quadratic', fill_value='extrapolate')

    # Remove x and y points outside +/0- limit min of the smoothed line
    mask = (y > p(x) - limit) & (y < p(x) + limit)
    xfilt = x[mask]
    yfilt = y[mask]
    wfilt = w[mask]

    # Refit
    # ynew2 = [np.mean(y[(x >= xi-0.05) & (x < xi + 0.05)]) for xi in xnew]
    ynew2 = []
    for i in range(len(xnew)):
        xi = xnew[i]
        yvals = yfilt[(xfilt >= xi - window * 0.5) & (xfilt < xi + window * 0.5)]
        wvals = wfilt[(xfilt >= xi - window * 0.5) & (xfilt < xi + window * 0.5)]
        if len(yvals) > 0:
            yavg = np.average(yvals, weights=wvals)
        else:
            yavg = np.nan
        ynew2.append(yavg)
    ynew2 = np.array(ynew2)
    # Remove NaN
    mask = ~np.isnan(ynew2)
    xnew2 = xnew[mask]
    ynew2 = ynew2[mask]

    # Smooth xnew and ynew using a gaussian filter with sigma of 2
    ynew2 = gaussian_filter1d(ynew2, sigma=smooth)
    # Create an interpolation of xnew and ynew
    p2 = interp1d(xnew2, ynew2, kind='quadratic', fill_value='extrapolate')

    if plot:
        # Plot the data and the fit
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, 'o', label='Data', alpha=0.25)
        plt.plot(xnew, ynew, 'o', label='Averaged Data', alpha=0.5)
        plt.plot(xnew2, ynew2, 'o', label='Averaged Data', alpha=0.5)
        # plt.plot(xnew, ynew, '-', label='Averaged Data')
        plt.plot(xnew, p(xnew), '--', label='Corrected Fit')
        # plt.plot(xnew2, ynew2, '-', label='Averaged Data 2')
        plt.plot(xnew2, p2(xnew2), '--', label='Corrected Fit 2')
        plt.plot(xnew, xnew, ':', label='y=x Line')
        plt.xlabel('Reference RT')
        plt.ylabel('Average RT')
        plt.legend()
        plt.show()

    return p2

def adjust_msp_file(poslibfile, neglibfile, interpolation_file, fileoutstr="_corrected", force_overwrite=False):
    # Check if interpolation file exists
    if not os.path.exists(interpolation_file):
        print(f"Interpolation file {interpolation_file} not found. Skipping RT correction.")
        return None, None
    with open(interpolation_file, 'rb') as f:
        p = pickle.load(f)

    # Check if poslibfile and neglibfile exist
    if not os.path.exists(poslibfile):
        print(f"Error: Positive library file {poslibfile} not found. Skipping RT correction for positive library.")
        posoutfile = None
    else:
        print("Adjusting positive library file: ", poslibfile)
        # Rename outfile to replace
        posoutfile = poslibfile.replace(".msp", fileoutstr + ".msp")

        # If outfile already exists and force_overwrite is False, skip
        if os.path.exists(posoutfile) and not force_overwrite:
            print(f"Output file {posoutfile} already exists. Skipping RT correction for positive library.")
            posoutfile = None

        else:
            # Read library file
            with open(poslibfile, 'r') as f:
                lines = f.readlines()

            # Write new library file correcting RTs
            with open(posoutfile, 'w') as f:
                for i, line in enumerate(lines):
                    if line.startswith("RETENTIONTIME:"):
                        rt = float(line.split(":")[1].strip())
                        newrt = p(rt)
                        f.write(f"RETENTIONTIME: {newrt:.4f}\n")
                        # print(f"Corrected RT: {rt} -> {newrt:.4f}")
                    else:
                        f.write(line)

    if not os.path.exists(neglibfile):
        print(f"Error: Negative library file {neglibfile} not found. Skipping RT correction for negative library.")
        negoutfile = None
    else:
        print("Adjusting negative library file: ", neglibfile)
        negoutfile = neglibfile.replace(".msp", fileoutstr + ".msp")
        # If outfile already exists and force_overwrite is False, skip
        if os.path.exists(negoutfile) and not force_overwrite:
            print(f"Output file {negoutfile} already exists. Skipping RT correction for negative library.")
            negoutfile = None
        else:
            with open(neglibfile, 'r') as f:
                lines = f.readlines()
            # Write new library file correcting RTs
            with open(negoutfile, 'w') as f:
                for i, line in enumerate(lines):
                    if line.startswith("RETENTIONTIME:"):
                        rt = float(line.split(":")[1].strip())
                        newrt = p(rt)
                        f.write(f"RETENTIONTIME: {newrt:.4f}\n")
                        # print(f"Corrected RT: {rt} -> {newrt:.4f}")
                    else:
                        f.write(line)
    return posoutfile, negoutfile

def find_lib_file(df, col="Comment"):
    posdf = df[df["Polarity"] == "Positive"]
    comment = posdf[col].iloc[0]
    # Extract the library name from the comment
    poslibname = comment.split("Annotation method: ")[1].split(";")[0]

    negdf = df[df["Polarity"] == "Negative"]
    comment = negdf[col].iloc[0]
    # Extract the library name from the comment
    neglibname = comment.split("Annotation method: ")[1].split(";")[0]

    # Remove _2 that is getting added to the library name for some reason
    poslibname = poslibname.replace("_2", "")
    neglibname = neglibname.replace("_2", "")

    # Add .msp to both library names
    poslibname += ".msp"
    neglibname += ".msp"
    # print(poslibname, neglibname)
    return poslibname, neglibname


def correct_dataset_rts(df, rtcol, ref_dataset="Stellar", namecol="Metabolite name", adductcol="Adduct type", plot=True):
    df = df.copy()
    refdf = df[df["Dataset"] == ref_dataset]

    if plot:
        # Plot 1-1 RT line for comparison
        plt.figure(figsize=(8, 8))
        plt.subplot(1,2,1)
        plt.plot([0, 20], [0, 20], "k--", label="1:1 line")
        plt.subplot(1,2,2)
        plt.plot([0, 20], [0, 20], "k--", label="1:1 line")
    # Take only unique names in the reference dataset, and their RTs
    for dataset in df["Dataset"].unique():
        if dataset == ref_dataset:
            continue
        # Get subset of data for this dataset
        subset = df[df["Dataset"] == dataset]

        xmatch = []
        ymatch = []
        for y in subset.iterrows():
            match = refdf[(refdf[namecol] == y[1][namecol]) & (refdf[adductcol] == y[1][adductcol])]
            if len(match) > 0:
                xmatch.append(match.iloc[0][rtcol])
                ymatch.append(y[1][rtcol])

        print(len(xmatch), len(ymatch))

        # Sort the matches by xmatch
        xmatch = np.array(xmatch)
        ymatch = np.array(ymatch)
        sort_idx = np.argsort(xmatch)
        xmatch = xmatch[sort_idx]
        ymatch = ymatch[sort_idx]
        # Drop any duplicates in xmatch and ymatch
        _, unique_idx = np.unique(xmatch, return_index=True)
        xmatch = xmatch[unique_idx]
        ymatch = ymatch[unique_idx]

        # Remove outliers that are more than 1 min away from the 1-1 line
        median_shift = np.median(ymatch - xmatch)
        print(f"Median RT shift for {dataset}: {median_shift:.2f} minutes")
        # Remove points that are more than 1 min away from the median shift
        mask = np.abs((ymatch - xmatch) - median_shift) < 1
        xmatch = xmatch[mask]
        ymatch = ymatch[mask]

        # cs = CubicSpline(xmatch, ymatch, bc_type='natural')
        fit = np.polyfit(ymatch, xmatch, deg=3)
        correction = np.poly1d(fit)
        # Define a function to apply the correction
        newrts = correction(df.loc[df["Dataset"] == dataset, rtcol])
        deltart = newrts - df.loc[df["Dataset"] == dataset, rtcol]
        print(f"Mean RT shift for {dataset}: {np.mean(deltart):.2f} minutes")
        df.loc[df["Dataset"] == dataset, rtcol] = newrts
        df.loc[df["Dataset"] == dataset, "RT Shift"] = deltart

        if plot:
            plt.subplot(1,2,1)
            # Plot these RTs to check for linearity
            plt.scatter(xmatch, ymatch, label=dataset)
            plt.xlabel("Reference RT")
            plt.ylabel("Dataset RT")
            print("SSE:", np.sum((xmatch - ymatch) ** 2))

            # Plot corrected values
            plt.subplot(1,2,2)
            plt.scatter(xmatch, correction(ymatch))
            print("SSE after correction:", np.sum((xmatch - correction(ymatch)) ** 2))
            plt.xlabel("Reference RT")
            plt.ylabel("Corrected Dataset RT")

    if plot:
        plt.subplot(1,2,1)
        plt.legend()
        plt.show()
    return df

if __name__ == "__main__":

    topdir = r"Z:\Group Share\Shivam\DDA"
    os.chdir(topdir)
    filepath = r"Combined_DDA_Full.xlsx"

    # Read in the combined DDA data with different sheets for each dataset
    dfs = pd.read_excel(filepath, sheet_name=None)

    # Individual datasets are keys with Total in them
    datasets = [key for key in dfs.keys() if "Total" in key]

    for dataset in datasets:
        print(f"Processing dataset: {dataset}")
        df = dfs[dataset]
        p2 = rt_correction(df, weight_col="Simple dot product", plot=False)

        # Save fit to p2 function to a pickle file
        with open(dataset + '_RTcorr.pkl', 'wb') as f:
            pickle.dump(p2, f)

    for dataset in datasets:
        if "Stellar" in dataset:
            print("Stellar dataset detected, skipping RT correction for library files.")
            print(find_lib_file(dfs[dataset]))
            continue
        interpolation_file = dataset + '_RTcorr.pkl'
        poslibfile, neglibfile = find_lib_file(dfs[dataset])
        posout, negout = adjust_msp_file(poslibfile, neglibfile, interpolation_file, fileoutstr="_" + dataset + "_DDA")
        if posout is None or negout is None:
            print(f"Skipping PG fix for {dataset} due to missing library files.")
            continue
        # If fixed not in name, also fix PG in the new library file
        if "fixed" not in negout.lower():
            fix_pg_in_msp(negout)


