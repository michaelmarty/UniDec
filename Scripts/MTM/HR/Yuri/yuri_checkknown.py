import os
import numpy as np
import pandas as pd

import unidec.IsoDec.compare_tools
import unidec.IsoDec.match as match
import matplotlib.pyplot as plt
import unidec.IsoDec.plots as isoplots
import unidec.tools as ud

dir = "C:\\Data\\Yuri\\"
os.chdir(dir)

infile = "known_masses.csv"
df = pd.read_csv(infile)
# Only keep lines where "Ion" inclucdes "c"
df = df[df["Ion"].str.contains("c", na=False)]

datafile = "CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_400k_Charge Reduction-yes_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL1000_768ms.csv"
datafile = "CA2_Sim_a_y_c_z_34+_LTQ-FT 21T_400k_Charge Reduction-yes_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL1000_768ms_NT3p5.csv"
datafile = "CA2_Sim_c_34+_LTQ-FT 21T_Charge Reduction-no_adduct_N-term_MP2_pi200_NoP_IntGaus_high_range_NL0.csv"


# data = np.loadtxt(datafile, delimiter=",", skiprows=1)

cdatafile = os.path.splitext(datafile)[0] + "_1_centroided.txt"
cdata = np.loadtxt(cdatafile)

# Convert to MatchedCollection
mc = unidec.IsoDec.compare_tools.df_to_matchedcollection(df, monoiso="Mass", peakmz="m/z", z="Charge", scan=None, intensity=None, peakmass=None, ion=None)

colormap = plt.get_cmap("rainbow_r")


cssdat = []
for pk in mc.peaks:
    isodist = pk.isodist
    # Match this against the data
    cutdata = ud.datachop(cdata, pk.mz - 1.5, pk.mz + 2.5)
    css = match.calc_css_from_data(cutdata, isodist)
    cssdat.append(css)
    pk.css = css
    # Set color based on css value
    color = colormap(css)
    pk.color = color

# plt.hist(cssdat, bins=50, color='blue', alpha=0.7)
# plt.xlabel('CSS Value')
# plt.ylabel('Count')
# plt.show()
# exit()

print("Number of peaks with CSS above 0.8:", len([pk for pk in mc.peaks if pk.css > 0.8]))
print("Out of total peaks:", len(mc.peaks))
print("Percent: ", len([pk for pk in mc.peaks if pk.css > 0.8]) / len(mc.peaks) * 100)

isoplots.plot_pks(mc, centroids=cdata, plotmass=False, ccolor="k")
plt.show()