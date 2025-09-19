# --- Imports ---
import rainbow as rb
from matplotlib import pyplot as plt
import re
import os

# --- Plotting Configuration ---
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28

# --- TUI for Filepaths ---
hplc_dir = input("Enter the HPLC data folder path: ").strip().replace('"', '')
figures_dir = os.path.join(hplc_dir, "Figures")
os.makedirs(figures_dir, exist_ok=True)  # Ensure Figures folder exists

# --- Find all .d folders ---
d_folders = [f for f in os.listdir(hplc_dir) if f.endswith('.D') and os.path.isdir(os.path.join(hplc_dir, f))]

for d_folder_name in d_folders:
    hplc_data_path = os.path.join(hplc_dir, d_folder_name)
    print(f"Processing: {hplc_data_path}")

    # --- Load HPLC Data ---
    datadir = rb.read(hplc_data_path)
    print(datadir)  # Inspect datadir object

    # --- Inspect Datafiles ---
    for datafile in datadir.datafiles:
        print(datafile.name, datafile.detector)
    print(datadir.get_info())  # Print additional info

    # --- Extract and Plot Data ---
    if datadir.datafiles:
        datafile = datadir.datafiles[0]
        times = datafile.xlabels         # Time points
        wavelengths = datafile.ylabels   # Wavelengths
        absorbances = datafile.data      # Absorbance values

        # Determine y-axis label from wavelength
        if wavelengths:
            match = re.search(r"\d+\.?\d*", str(wavelengths[0]))
            if match:
                wavelength_nm = float(match.group())
                ylabel = f"Absorbance at {wavelength_nm:.0f} nm (mAu)"
            else:
                ylabel = "Absorbance (mAu)"
        else:
            ylabel = "Absorbance (mAu)"

        # Plot absorbance vs. time
        plt.figure(figsize=(16, 9))
        plt.plot(times, absorbances)
        plt.xlabel("Time (minutes)")
        plt.ylabel(ylabel)

        # --- Save Figure ---
        figure_filename = f"Absorbance Plot {d_folder_name}.png"
        figure_path = os.path.join(figures_dir, figure_filename)
        plt.savefig(figure_path, dpi=300, bbox_inches="tight")
        plt.close()