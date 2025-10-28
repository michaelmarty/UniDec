# --- Import necessary libraries ---
import rainbow as rb
from matplotlib import pyplot as plt
import re
import os

# --- Set matplotlib font properties ---
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# --- Prompt user for main HPLC data folder and create Figures directory ---
hplc_dir = input("Enter the HPLC data folder path: ").strip().replace('"', '')
figures_dir = os.path.join(hplc_dir, "Figures")
os.makedirs(figures_dir, exist_ok=True)

# --- Ask user how many .D folders to overlay and collect their paths ---
num_folders = int(input("How many .D folders would you like to overlay? "))
d_folder_paths = []
for i in range(num_folders):
    folder_path = input(f"Enter the full path for .D folder #{i+1}: ").strip().replace('"', '')
    d_folder_paths.append(folder_path)

# --- Create a new figure for overlaying absorbance plots ---
plt.figure(figsize=(16, 9))

# --- Process each .D folder and plot its absorbance data ---
for d_folder_path in d_folder_paths:
    print(f"Processing: {d_folder_path}")
    datadir = rb.read(d_folder_path)
    print(datadir)
    for datafile in datadir.datafiles:
        print(datafile.name, datafile.detector)
    print(datadir.get_info())

    # --- Plot only if there is at least one datafile ---
    if datadir.datafiles:
        datafile = datadir.datafiles[0]
        times = datafile.xlabels
        wavelengths = datafile.ylabels
        absorbances = datafile.data

        # --- Determine y-axis label based on wavelength information ---
        if wavelengths:
            match = re.search(r"\d+\.?\d*", str(wavelengths[0]))
            if match:
                wavelength_nm = float(match.group())
                ylabel = f"Absorbance at {wavelength_nm:.0f} nm (mAu)"
            else:
                ylabel = "Absorbance (mAu)"
        else:
            ylabel = "Absorbance (mAu)"

        # --- Plot the absorbance data ---
        plt.plot(times, absorbances, label=os.path.basename(d_folder_path))

# --- Finalize and save the figure ---
plt.xlabel("Time (minutes)")
plt.ylabel(ylabel)
plt.legend()
figure_filename = "Absorbance Overlay Plot.png"
figure_path = os.path.join(figures_dir, figure_filename)
plt.savefig(figure_path, dpi=300, bbox_inches="tight")
plt.close()