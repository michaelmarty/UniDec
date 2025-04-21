from LCImporter import *
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.sans-serif'] = "Arial"


p1 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240816_Estevan's sample test 2024-08-16 12-05-27\\20240816_sf GFP.D"

p2 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240816_Estevan's sample test 2024-08-16 12-26-58\\20240816_sf GFP 10 ul.D"

p3 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240816_Estevan's sample test 2024-08-16 12-44-04\\20240816_sf GFP 20 ul.D"

p4 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240816_Estevan's sample test 2024-08-16 13-02-33\\20240816_sf GFP 30 ul.D"

p5 = "C:\\Users\\Public\\Documents\\ChemStation\\1\\Data\\OF_20240816_Estevan's sample test 2024-08-16 13-21-38\\20240816_sf GFP 40 ul.D"

d1 = extract_data(p1, index = 0)
d2 = extract_data(p2, index = 0)
d3 = extract_data(p3, index = 0)
d4 = extract_data(p4, index = 0)
d5 = extract_data(p5, index = 0)

plt.plot(d1[:,0], d1[:,1], label="5 ul sf GFP", color="black")



#plt.subplot(3, 1, 2)
plt.plot(d2[:,0], d2[:,1], label="10 ul sf GFP", color="blue")



#plt.subplot(3, 1, 3)
plt.plot(d3[:,0], d3[:,1], label="20 ul sf GFP", color="red")


#plt.subplot(3, 1, 3)
plt.plot(d4[:,0], d4[:,1], label="30 ul sf GFP", color="orchid")


#plt.subplot(3, 1, 3)
plt.plot(d5[:,0], d5[:,1], label="40 ul sf GFP", color="green")

plt.ylim(0.19, 3)
plt.xlim(0, 7)

# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel('Retention time (min)')
plt.ylabel('FL (LU)')
plt.title('Minimum volume required in sample vial')
plt.legend()

# Adjust layout
plt.tight_layout()

# Save figure
plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240816_minimum volume in vial.png", dpi=1000)


plt.show()