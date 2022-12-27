import os
import unidec.modules.unidectools as ud
from unidec.modules.waters_importer.WatersImporter import WatersDataImporter as WDI
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import pandas as pd

def correct_mass(data, exact, window=1):
    """
    Function to internally calibrate the data to a known m/z
    :param data: m/z in column 0, intensity in column 1, numpy array
    :param exact: float, m/z value of the exact mass
    :param window: float, m/z window for considering the correction. The exact mass must be within this window.
    """
    d2 = ud.datachop(data, exact - window, exact + window) # Chop the data down to the +/- window around exact
    maxval = d2[np.argmax(d2[:, 1]), 0] # Find the maximum m/z in this window
    data[:, 0] = data[:, 0] + exact - maxval #correct the m/z values so that the peak matches exact
    return data # return the data


def get_vars(path):
    """
    Function to help decode the file names into lipid types and ratios
    :param path: the file path
    """
    lipids = ["DMPC", "DMPG", "Ecoli", "DPPC","POPC", "DLPC", "Free",  "FDAP", "fdap", "fdapb", "Only"] # Potential lipid string matches
    ratios = ["3to", "9to", "18to"] # Potential ratio matches
    freenames = ["free", "fdap", "fdapb", "only"]
    # Initialize
    lipid = ""
    ratio = ""
    # Loop over the lipid strings and see if any fit
    for l in lipids:
        if l.lower() in path.lower():
            lipid = l
            if l.lower() in freenames:
                lipid = "Free"
    # Loop over the ratios and see if any fit
    for r in ratios:
        if r.lower() in path.lower():
            ratio = r[:-2]
    print(lipid, ratio)
    return lipid, ratio


matplotlib.use('WxAgg')

# Directory Access
#topdir = "Z:\\Group Share\\Dash\\Data\\Synapt Data\\Daptomycin\\DaptoScan4Oct"
topdir = "Z:\\Group Share\\Deseree\\FPOP Paper 2\\For Public Repository\\FPOP Data"
#files = glob.glob(f'{topdir}/*.raw', recursive=True)
files = ud.match_dirs_recursive(topdir, ".raw")
print(files)

# Set constants
monomass = 1619.710
carbmass = 12.000
fluormass = 18.998403
iodomass = 126.904477
z = 2

# Set potential m/z values to look for
# NOTE: These are meant to be flexible and be able to search for other potential mods for future studies
mz = (monomass + z * ud.protmass) / z
mzox = mz + (ud.oxmass / z) #(+16)
mzox1 = mz - (ud.oxmass / z) #(-16)
mzox2 = mz + ((ud.oxmass * 2) / z) #(+32)
mzox3 = mz + ((ud.oxmass * 3) / z) #(+48)
mzox4 = mz - ((ud.oxmass * 2) / z) #(-32)
mzox5 = mz - ((ud.oxmass + 6 * ud.protmass) / z) #(-22)
mzox6 = mz - ((ud.oxmass - 6 * ud.protmass) / z) #(-10)
mzox7 = mz - (5 * ud.protmass / z) #(+5)
mzox8 = mz + ((ud.oxmass - 2 * ud.protmass) / z) #(+14)
mzox9 = mz - ((ud.oxmass * 3 - 5 * ud.protmass) / z) #(-43)
mzox10 = mz - ((2 * ud.oxmass - 2 * ud.protmass)/ z) #(-30)
mzox11 = mz - (2 * ud.protmass / z) #(-2)
mzox12 = mz - ((ud.oxmass + 7 * ud.protmass) / z) #(-23)
mzox13 = mz + (4 * ud.protmass / z) #(+4)

mzfluorox = mz + ((carbmass + 3 * fluormass) / z)
mziodox = mz + (iodomass / z)

oxlist = [mzox, mzox1, mzox2, mzox3, mzox4, mzox5, mzox6, mzox7, mzox8, mzox9, mzox10, mzox11, mzox12, mzox13, mzfluorox, mziodox]
names = ["+Ox", "-Ox", "+2Ox", "+3Ox", "-2Ox", "-Ox-6H", "-Ox+6H", "-5H", "+Ox-2H", "-3Ox+5H", "-Ox-2H", "-2H", "-Ox7H", "+4H", "+C+3F", "+I"]

# Set windows to consider
window = 0.1
window1 = 0.1
window2 = 0.1
time_range = [9, 12]

# Set outputs
output = []
outfile = "spectrum.npz"
# Loop through files
for i, path in enumerate(files[:]):
    print(path)
    os.chdir(path)

    # Check if there is a spectrum file already produced
    # If there is load it
    # If not, average the time range to produce a spectrum file
    if os.path.isfile(outfile):
        data = np.load(outfile, allow_pickle=True)['data']
    else:
        try:
            d = WDI(path)
        except Exception as e:
            print("File Failed", path)
            print(e)

        data = d.get_data(time_range=time_range)
        np.savez_compressed(outfile, data=data)

    # Correct the mass
    data = correct_mass(data, mz, window=1)
    data = ud.datachop(data, 805, 830)

    # Integrate the main peak
    mzint = ud.integrate(data, mz - window1, mz + window1)[0]
    # Integrate the other peaks
    oxints = np.array([ud.integrate(data, ox-window1, ox+window1)[0] for ox in oxlist])
    reloxs = (oxints/mzint) * 100

    # Some sanity check to make sure that there is some reasonable peak area detected
    if mzint < 10000:
        print("ERROR: ", path)

    # Parse file name to get the correct lipid and ratio from the file names
    lipid, ratio = get_vars(path)

    # Collect the data
    out1 = np.array([path, lipid, ratio])
    out2 = reloxs.astype(str)
    out3 = np.concatenate((out1, out2))
    output.append(out3)
    print(out3)

    # Plot the restuls
    plt.plot(data[:, 0], data[:, 1] / np.amax(data[:, 1]) - i * 0.05)
    plt.vlines(mz - window2, -7, 1)
    plt.vlines(mz + window1, -7, 1)
    plt.vlines(mzox - window2, -7, 1)
    plt.vlines(mzox + window1, -7, 1)
    # plt.xlim(805, 830)

# Save the outputs
os.chdir(topdir)
# Create the column heads
head = ["File Name", "Lipid Type", "Ratio"]
header = head + names
print(header)
df = pd.DataFrame(output, dtype=str, columns=header)
df.to_excel("RawExtractedPercents.xlsx")

# Show the plots
plt.show()
