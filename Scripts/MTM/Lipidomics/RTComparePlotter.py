import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_smoothing_spline
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

# Create a function to correct the Reference RT to match the Average using a 4th degree polynomial
fit = np.polyfit(df['Reference RT'], df['Average Rt(min)'], 8)
p = np.poly1d(fit)
# Sort by Reference RT
df = df.sort_values(by=['Reference RT'], ascending=True)
x = np.array(df['Reference RT'])
y = np.array(df['Average Rt(min)'])
p = make_smoothing_spline(x, y, w=None, lam=0.5)
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