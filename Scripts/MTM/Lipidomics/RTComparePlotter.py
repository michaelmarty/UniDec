import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import os

data = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\RTComparison.csv"
posdata = r"Z:\Group Share\MTM\IQ-X\MartyLab_AX\RTComparisonPos.csv"

os.chdir("Z:\\Group Share\\MTM\\IQ-X\\MartyLab_AX")

df = pd.read_csv(data)
# Drop where Library RT is 0
df = df[df['Library RT'] > 0]

dfpos = pd.read_csv(posdata)
# Drop where Library RT is 0
dfpos = dfpos[dfpos['Library RT'] > 0]

# Contatenate dataframes
df = pd.concat([df, dfpos])

# Append a line beyond this to force fit to go through these points
# List of rows to add
new_rows = [
    {'Name': 'End1', 'Average Rt(min)': 30, 'Reference RT': 30, 'Library RT': 30},
    {'Name': 'End2', 'Average Rt(min)': 20, 'Reference RT': 20, 'Library RT': 20},
     {'Name': 'End3', 'Average Rt(min)': 25, 'Reference RT': 25, 'Library RT': 25},
{'Name': 'End3', 'Average Rt(min)': 22, 'Reference RT': 22, 'Library RT': 22},
{'Name': 'End3', 'Average Rt(min)': 27, 'Reference RT': 27, 'Library RT': 27}
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
ynew = gaussian_filter1d(ynew, sigma=1)

# Create an interpolation of xnew and ynew
p = interp1d(xnew, ynew, kind='linear', fill_value='extrapolate')

# Remove x and y points outside +/0- 3 min of the smoothed line
mask = (y > p(x) - 3) & (y < p(x) + 3)
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
ynew2 = gaussian_filter1d(ynew2, sigma=1)
# Create an interpolation of xnew and ynew
p2 = interp1d(xnew2, ynew2, kind='quadratic', fill_value='extrapolate')

# Save fit to p2 function to a pickle file
import pickle
with open('RTInterpolation.pkl', 'wb') as f:
    pickle.dump(p2, f)



# Plot the data and the fit
plt.figure(figsize=(10, 6))


plt.plot(x, y, 'o', label='Data')
plt.plot(xnew, ynew, '-', label='Averaged Data')
plt.plot(xnew, p(xnew), '--', label='Corrected Fit')
plt.plot(xnew2, ynew2, '-', label='Averaged Data 2')
plt.plot(xnew2, p2(xnew2), '--', label='Corrected Fit 2')
plt.plot(xnew, xnew, ':', label='y=x Line')
plt.show()
exit()


df['Corrected RT'] = p(df['Reference RT'])
print(fit)



# Plot Average Rt(min) vs Reference RT
plt.figure(figsize=(10, 6))
plt.subplot(131)
plt.scatter(df['Reference RT'], df['Average Rt(min)'], color='blue', label='Data Points')
plt.plot([df['Reference RT'].min(), df['Reference RT'].max()],
         [df['Reference RT'].min(), df['Reference RT'].max()],
         color='red', linestyle='--', label='y=x Line')
plt.scatter(df['Reference RT'], df['Corrected RT'], color='orange', linestyle='-', label='4th Degree Fit')
plt.plot(np.arange(5,15), p(np.arange(5,15)), color='green', linestyle=':', label='Corrected Line')
plt.xlabel('Reference RT (min)')
plt.ylabel('Average Rt (min)')
plt.title('Average Rt vs Reference RT')
plt.legend()
plt.grid()

# Plot Average Rt(min) vs Library RT
plt.subplot(132)
plt.scatter(df['Library RT'], df['Average Rt(min)'], color='green', label='Data Points')
plt.plot([df['Library RT'].min(), df['Library RT'].max()],
         [df['Library RT'].min(), df['Library RT'].max()],
         color='red', linestyle='--', label='y=x Line')
plt.xlabel('Library RT (min)')
plt.ylabel('Average Rt (min)')
plt.title('Average Rt vs Library RT')
plt.legend()
plt.grid()

# Plot Reference RT vs Library RT
plt.subplot(133)
plt.scatter(df['Reference RT'], df['Library RT'], color='purple', label='Data Points')
plt.plot([df['Reference RT'].min(), df['Reference RT'].max()],
         [df['Reference RT'].min(), df['Reference RT'].max()],
         color='red', linestyle='--', label='y=x Line')
plt.xlabel('Reference RT (min)')
plt.ylabel('Library RT (min)')
plt.title('Reference RT vs Library RT')
plt.legend()
plt.grid()

plt.show()