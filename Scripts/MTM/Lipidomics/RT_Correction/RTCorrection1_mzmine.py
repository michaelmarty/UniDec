import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import os
import pickle

def cleanup_df(df):
    # Remove rows with null in MS/MS Spectrum
    df = df[df['MSMS spectrum'].notnull()]
    # Remove name if Unknown
    df = df[~df['Name'].str.contains("Unknown", na=False)]
    df = df[~df['Name'].str.contains("low score", na=False)]
    return df

def import_data(posfile, negfile=None):
    posdf = pd.read_csv(posfile)
    # posdf = cleanup_df(posdf)
    # Grab relevant columns: Metabolite name, Average Rt(min), Reference RT, S/N average, Total score
    posrtdf = posdf[['Metabolite name', 'Average Rt(min)', 'Reference RT']].reset_index(drop=True)

    if negfile is not None:
        negdf = pd.read_csv(negfile)
        # negdf = cleanup_df(negdf)
        negrtdf = negdf[['Metabolite name', 'Average Rt(min)', 'Reference RT']].reset_index(drop=True)

        # Contatenate dataframes
        df = pd.concat([negrtdf, posrtdf])
    else:
        df = posrtdf

    # Append a line beyond this to force fit to go through these points
    # List of rows to add
    template_row = {'Metabolite name': 'End1', 'Average Rt(min)': 0, 'Reference RT': 0,}
    # For 0 to min rt, add x=y
    min_rt = df['Reference RT'].min()
    new_rows = []
    for rt in np.arange(0, min_rt, 0.1):
        row = template_row.copy()
        row['Average Rt(min)'] = rt
        row['Reference RT'] = rt
        new_rows.append(row)

    max_rt = df['Reference RT'].max()
    for rt in np.arange(max_rt, max_rt + 10, 0.1):
        row = template_row.copy()
        row['Average Rt(min)'] = rt
        row['Reference RT'] = rt
        new_rows.append(row)

    # Convert to DataFrame and concatenate
    df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

    return df


def rt_fit(df, smooth=1, window=0.5, plot=True):
    # Create a function to correct the Reference RT to match the Average using a 4th degree polynomial
    fit = np.polyfit(df['Reference RT'], df['Average Rt(min)'], 2)
    p = np.poly1d(fit)
    # Sort by Reference RT
    df = df.sort_values(by=['Reference RT'], ascending=True)
    x = np.array(df['Reference RT'])
    y = np.array(df['Average Rt(min)'])
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

    # Save fit to p2 function to a pickle file
    with open('RTInterpolation.pkl', 'wb') as f:
        pickle.dump(p2, f)

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
        plt.legend()
        plt.show()

def rt_correction(poslibfile, neglibfile=None, corr_fuction=None, outputext="-Apex"):
    if corr_fuction is None:
        interpolationfunction = r"RTInterpolation.pkl"
        with open(interpolationfunction, 'rb') as f:
            p = pickle.load(f)
    else:
        p = corr_fuction

    # Rename outfile to replace C1 with Stellar 1
    posoutfile = poslibfile.replace(".msp", outputext+ ".msp")

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

    if neglibfile is not None:

        with open(neglibfile, 'r') as f:
            lines = f.readlines()

        negoutfile = neglibfile.replace(".msp", outputext + ".msp")

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

    print(f"Corrected library files saved as {posoutfile} and {negoutfile if neglibfile is not None else 'N/A'}")


os.chdir(r"Z:\Group Share\Annika\DDA_Compare\mzmine")

posfile = r"ITIT_Pos_features_data_msdial.csv"
negfile = None

poslibfile = r"..\MSDIAL-TandemMassSpectralAtlas-VS69-Pos-StellarRTP2-IQX_AqX.msp"
poslibfile = None
# poslibfile = r"..\MSDIAL-TandemMassSpectralAtlas-VS69-Neg-RTP2_Fixed_no_hcoo.msp"
neglibfile = None

extadd = "-IQX_AqX"

limit = 5
smooth= 1
window = 0.5

df = import_data(posfile, negfile)
p2 = rt_fit(df, smooth=smooth, window=window, plot=True)


if poslibfile is not None:
    do_corr = input("Do you want to apply this correction to the library files? (y/n): ")
    if do_corr.lower() == 'y':
        rt_correction(poslibfile, neglibfile, corr_fuction=p2, outputext=extadd)


exit()


