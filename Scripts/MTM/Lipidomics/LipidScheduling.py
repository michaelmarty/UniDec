import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

topdir = r"Y:\Data\Lipidomics\Libraries"

os.chdir(topdir)

df = pd.read_csv("neglib.csv")

print(len(df))
dropclasses = ["PIM", "EBDS", "AHex"]
# Drop rows where the class contains any of the dropclasses as substrings
df = df[~df['compoundclass'].str.contains('|'.join(dropclasses), case=False, na=False)]
print(len(df))
# Print unique compoundclass values
print(df['compoundclass'].unique())

# Create a new column with rounded values of m/z
df['mz_rounded'] = df['precursor_mz'].round()

# create a colum for rounded retention_time values to the first decimal
df['rt_rounded'] = df['retention_time'].round(1)

# Create bins for m/z and retention time
mzbins = np.arange(df['mz_rounded'].min(), df['mz_rounded'].max() + 1, 1)
rtbins = np.arange(df['rt_rounded'].min(), df['rt_rounded'].max() + 0.1, 0.1)

# Create a 2D histogram of mz_rounded and rt_rounded
hist, xedges, yedges = np.histogram2d(df['mz_rounded'], df['rt_rounded'], bins=[mzbins, rtbins])

# At each time point, count the number of non-zero values in the histogram
# This will give us the number of unique m/z values at each retention time
unique_counts = np.count_nonzero(hist, axis=0)
plt.plot(yedges[:-1], unique_counts, marker='o', linestyle='-', color='blue')
plt.show()

exit()
# Plot the 2D histogram
plt.figure(figsize=(10, 6))
plt.imshow(np.log(hist.T), origin='lower', aspect='auto', cmap='viridis',
           extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='Counts')
plt.xlabel('m/z (rounded)')
plt.ylabel('Retention Time (rounded)')
plt.title('2D Histogram of m/z and Retention Time')
plt.show()

print(df)