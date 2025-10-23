import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import os


os.chdir(r"Z:\Group Share\Annika\Stellar\HEK")

posfile = r"20251023_AS_HEK_PosDDA_MSDial.csv"
negfile = r"20251023_AS_HEK_NegDDA_MSDial.csv"

limit = 2
smooth= 2

posdf = pd.read_csv(posfile)
posrtdf = posdf[['Name', 'RT (min)', 'Reference RT']].rename(columns={'RT (min)': 'Average Rt(min)'}).reset_index(drop=True)

negdf = pd.read_csv(negfile)
negrtdf = negdf[['Name', 'RT (min)', 'Reference RT']].rename(columns={'RT (min)': 'Average Rt(min)'}).reset_index(drop=True)

# Contatenate dataframes
df = pd.concat([negrtdf, posrtdf])

# Append a line beyond this to force fit to go through these points
# List of rows to add
new_rows = [
    {'Name': 'End1', 'Average Rt(min)': 30, 'Reference RT': 30},
    {'Name': 'End2', 'Average Rt(min)': 20, 'Reference RT': 20},
     {'Name': 'End3', 'Average Rt(min)': 25, 'Reference RT': 25},
{'Name': 'End3', 'Average Rt(min)': 22, 'Reference RT': 22},
{'Name': 'End3', 'Average Rt(min)': 27, 'Reference RT': 27}
]


# Convert to DataFrame and concatenate
df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)


# Create a function to correct the Reference RT to match the Average using a 4th degree polynomial
fit = np.polyfit(df['Reference RT'], df['Average Rt(min)'], 5)
p = np.poly1d(fit)
# Sort by Reference RT
df = df.sort_values(by=['Reference RT'], ascending=True)
x = np.array(df['Reference RT'])
y = np.array(df['Average Rt(min)'])

minx = np.min(x)
maxx = np.max(x)
# Create linear spacing of 0.1 min between min and max
xnew = np.arange(minx, maxx, 0.1)
# Average each y for each xnew
ynew = [np.mean(y[(x >= xi-0.05) & (x < xi + 0.05)]) for xi in xnew]
ynew = np.array(ynew)
# Remove NaN
mask = ~np.isnan(ynew)
xnew = xnew[mask]
ynew = ynew[mask]

# Smooth xnew and ynew using a gaussian filter with sigma of 2
ynew = gaussian_filter1d(ynew, sigma=smooth)

# Create an interpolation of xnew and ynew
p = interp1d(xnew, ynew, kind='linear', fill_value='extrapolate')

# Remove x and y points outside +/0- limit min of the smoothed line
mask = (y > p(x) - limit) & (y < p(x) + limit)
x = x[mask]
y = y[mask]

# Refit
ynew2 = [np.mean(y[(x >= xi-0.05) & (x < xi + 0.05)]) for xi in xnew]
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
plt.plot(xnew, ynew, '-', label='Averaged Data')
# plt.plot(xnew, p(xnew), '--', label='Corrected Fit')
plt.plot(xnew2, ynew2, '-', label='Averaged Data 2')
# plt.plot(xnew2, p2(xnew2), '--', label='Corrected Fit 2')
plt.plot(xnew, xnew, ':', label='y=x Line')
plt.legend()
plt.show()
exit()


exit()
data = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\NegIDs.csv"
os.chdir("Z:\\Group Share\\MTM\\IQ-X\\MartyLab_AX")

df = pd.read_csv(data)

libfile = r"C:\Users\mm96978\Downloads\MSDIAL-TandemMassSpectralAtlas-VS69-Neg.msp"

# Read library file
with open(libfile, 'r') as f:
    lines = f.readlines()

# Loop through df
rts = []
for index, row in df.iterrows():
    name = row['Metabolite name']
    rt = row['Average Rt(min)']
    refrt = row['Reference RT']
    librt = 0
    try:
        n2 = name.split("|")[1]
    except:
        n2 = name
    try:
        # Replace O2 and O3 with 2O and 3O
        n2 = n2.replace("O2", "2O").replace("O3", "3O")
    except:
        pass
    print(n2)
    # Find name in library
    for i, line in enumerate(lines):
        if line.startswith("NAME: "+n2):
            print(line)
            # Look for RT line
            for j in range(i, min(i+20, len(lines))):
                if lines[j].startswith("RETENTIONTIME:"):
                    librt = float(lines[j].split(":")[1].strip())
                    print("Found RT:", librt)
                    break
            break
    rts.append([name, rt, refrt, librt])

# Convert to DataFrame
rtdf = pd.DataFrame(rts, columns=['Name', 'Average Rt(min)', 'Reference RT', 'Library RT'])
print(rtdf)

# Save to CSV
rtdf.to_csv("RTComparison.csv", index=False)