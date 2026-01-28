import pandas as pd
import os
from Scripts.LipiDec.AddFeaturestoDF import *
import matplotlib.pyplot as plt
import numpy as np

topdir = r"Y:\Data\Lipidomics\Libraries\LCpredefined"

# Find all CSV files in this directory and its subdirectories
files = []
for root, dirs, filenames in os.walk(topdir):
    for filename in filenames:
        if filename.endswith(".csv"):
            files.append(os.path.join(root, filename))

taillist = []
classlist = []

for f in files:
    print(f)
    df = pd.read_csv(f, header=None, names=["Metabolite name"])
    df = set_class_name(df)
    df = set_tails(df, vectors=True)
    tl = get_tail_list(df)
    for t in tl:
        taillist.append(t)

    cl = get_class_list(df)
    for c in cl:
        classlist.append(c)

taillist = np.array(taillist)
classlist = np.array(classlist)
#
# lib = r"Y:\Data\Lipidomics\Libraries\LCpredefined\Mouse\Heart.csv"
#
# # Read list and set the column name as "Metabolite name"
# df = pd.read_csv(lib, header=None, names=["Metabolite name"])
# df = set_tails(df, vectors=True)
# print(df)
# taillist = get_tail_list(df)

tlrange = [np.amin(taillist[:,0]), np.amax(taillist[:,0])]
tursange = [np.amin(taillist[:,1]), np.amax(taillist[:,1])]

# Get histogram of tail lengths and unsaturations
tailhist, xedges, yedges = np.histogram2d(taillist[:,0], taillist[:,1], bins=[np.arange(tlrange[0], tlrange[1]+2)-0.5, np.arange(tursange[0], tursange[1]+2)-0.5])
print(tailhist)

tailhits = []
# For all hist values that aren't zero, add them to a list
for i in range(tailhist.shape[0]):
    for j in range(tailhist.shape[1]):
        if tailhist[i,j] > 0:
            tailhits.append((xedges[i]+0.5, yedges[j]+0.5, tailhist[i,j]))
tailhits = np.array(tailhits)
print(tailhits, len(tailhits))

ucl = np.unique(classlist)
print(ucl)

# Create a long string from the tailhits for pasting into excel
tailhitstr = ""
for th in tailhits:
    tailhitstr += f"{int(th[0])}:{int(th[1])}\t{int(th[2])}\n"
print(tailhitstr)
exit()

# create a 2D histogram of the tail lengths and unsaturations
plt.figure(figsize=(10, 6))
plt.hist2d(taillist[:,0], taillist[:,1], bins=[np.arange(tlrange[0], tlrange[1]+2)-0.5, np.arange(tursange[0], tursange[1]+2)-0.5], cmap='viridis')
plt.colorbar(label='Counts')
plt.xlabel('Tail Length')
plt.ylabel('Tail Unsaturations')
plt.title('2D Histogram of Tail Lengths and Unsaturations')
plt.xticks(np.arange(tlrange[0], tlrange[1]+1))
plt.yticks(np.arange(tursange[0], tursange[1]+1))
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()



