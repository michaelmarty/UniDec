import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import os

def cleanup_df(df):
    # Remove rows with null in MS/MS Spectrum
    df = df[df['MSMS spectrum'].notnull()]
    # Remove name if Unknown
    df = df[~df['Name'].str.contains("Unknown", na=False)]
    df = df[~df['Name'].str.contains("low score", na=False)]
    return df

os.chdir(r"Z:\Group Share\Annika\Stellar\HEK")

posfile = r"20251023_AS_HEK_PosDDA_MSDial.csv"
negfile = r"20251023_AS_HEK_NegDDA_MSDial.csv"

posfile = "20250929_AS_HEK_DDA_20_Pos.txt"
negfile = "20250929_AS_HEK_DDA_20_Neg.txt"

limit = 5
smooth= 1
window = 0.5

posdf = pd.read_table(posfile)
posdf = cleanup_df(posdf)

posrtdf = posdf[['Name', 'RT (min)', 'Reference RT', "Area", "Total score"]].rename(columns={'RT (min)': 'Average Rt(min)'}).reset_index(drop=True)

negdf = pd.read_table(negfile)
negdf = cleanup_df(negdf)
negrtdf = negdf[['Name', 'RT (min)', 'Reference RT', "Area", "Total score"]].rename(columns={'RT (min)': 'Average Rt(min)'}).reset_index(drop=True)

# Contatenate dataframes
df = pd.concat([negrtdf, posrtdf])
# df = negrtdf

# Append a line beyond this to force fit to go through these points
# List of rows to add
new_rows = [
    {'Name': 'End1', 'Average Rt(min)': 30, 'Reference RT': 30, "Area":1, "Total score":100},
    {'Name': 'End2', 'Average Rt(min)': 20, 'Reference RT': 20, "Area":1, "Total score":100},
     {'Name': 'End3', 'Average Rt(min)': 25, 'Reference RT': 25, "Area":1, "Total score":100},
{'Name': 'End3', 'Average Rt(min)': 22, 'Reference RT': 22, "Area":1, "Total score":100},
{'Name': 'End3', 'Average Rt(min)': 27, 'Reference RT': 27, "Area":1, "Total score":100},
]


# Convert to DataFrame and concatenate
df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)


# Create a function to correct the Reference RT to match the Average using a 4th degree polynomial
fit = np.polyfit(df['Reference RT'], df['Average Rt(min)'], 2)
p = np.poly1d(fit)
# Sort by Reference RT
df = df.sort_values(by=['Reference RT'], ascending=True)
x = np.array(df['Reference RT'])
y = np.array(df['Average Rt(min)'])
w = np.array(df['Total score'])

minx = np.min(x)
maxx = np.max(x)
print(minx, maxx)
# Create linear spacing of window min between min and max
xnew = np.arange(minx, maxx, window)
# Average each y for each xnew
ynew = []
wnew = []
for i in range(len(xnew)):
    xi = xnew[i]
    yvals = y[(x >= xi-window*0.5) & (x < xi + window*0.5)]
    wvals = w[(x >= xi-window*0.5) & (x < xi + window*0.5)]
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
print(np.amin(xnew), np.amax(xnew))


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
    yvals = yfilt[(xfilt >= xi-window*0.5) & (xfilt < xi + window*0.5)]
    wvals = wfilt[(xfilt >= xi-window*0.5) & (xfilt < xi + window*0.5)]
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
import pickle
with open('RTInterpolation.pkl', 'wb') as f:
    pickle.dump(p2, f)

# Plot the data and the fit
plt.figure(figsize=(10, 6))


plt.plot(x, y, 'o', label='Data', alpha=0.25)
plt.plot(xnew , ynew, 'o', label='Averaged Data', alpha=0.5)
plt.plot(xnew2 , ynew2, 'o', label='Averaged Data', alpha=0.5)
# plt.plot(xnew, ynew, '-', label='Averaged Data')
plt.plot(xnew, p(xnew), '--', label='Corrected Fit')
# plt.plot(xnew2, ynew2, '-', label='Averaged Data 2')
plt.plot(xnew2, p2(xnew2), '--', label='Corrected Fit 2')
plt.plot(xnew, xnew, ':', label='y=x Line')
plt.legend()
plt.show()
exit()


