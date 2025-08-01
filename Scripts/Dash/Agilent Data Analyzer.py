import os
import numpy as np
import rainbow as rb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import colormaps
from matplotlib import rcParams
import pandas as pd
import colorsys
import glob

topdir = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\20250729 TD UIB_FLDTest_MultiplexBoost 2025-07-30 11-12-02"
files = glob.glob(os.path.join(topdir, "*.D"))
if not files:
    raise FileNotFoundError("No .D files found in the specified directory.")

channel_indices = [3, 4] # Define the channel indices for multi-channel extraction only. Numbering starts from 0, so 3 and 4 correspond to the 4th and 5th channels in the data files.

def get_filtered_colors(cmap_name='cool', num_colors=256, luminosity_factor=0.8, saturation_factor=2, hue_shift=0):
    #A value of 0.5 would reduce the luminosity/sat to half, while 1.0 would keep it unchanged,  A value greater than 1.0 will increase saturation, while 1.0 will keep it unchanged.

    viridis = colormaps.get_cmap(cmap_name)
    colors = []

    for i in range(num_colors):
        color = viridis(i / num_colors)  # Sample uniformly
        r, g, b, _ = color  # Unpack the color and ignore the alpha value
        h, l, s = colorsys.rgb_to_hls(r, g, b)


        new_hue = (h + hue_shift) % 1.0
        new_light = max(0, l * luminosity_factor)  # Ensure it doesn't go below 0
        new_sat = min(1, s * saturation_factor)  # Ensure it doesn't go below 0
        new_color = colorsys.hls_to_rgb(new_hue, new_light, new_sat)

        colors.append((new_color[0], new_color[1], new_color[2], color[3]))  # Retain the original alpha

    return colors


def get_UIB1B_files(files):
    for path in files:
        print(f"Processing: {path}")
        datadir = rb.read(path)
        #print(datadir)

        datafile = datadir.get_file("UIB1B.ch")
        time = np.array(datafile.xlabels)
        #print(time)
        voltages = np.array(datafile.data)
        #print(voltages)

        print("time shape:", time.shape)
        print("voltages shape:", voltages.shape)

        outdir = os.path.basename(path).replace(".D", "")
        os.makedirs(outdir, exist_ok=True)
        np.savetxt(os.path.join(outdir, "UIB1B_data.txt"), np.column_stack((time, voltages)), header="Time (s)\tVoltage (V)", delimiter="\t")
        np.savez_compressed(os.path.join(outdir, "UIB1B_data.npz"), time, voltages)

        plt.figure(figsize=(16, 9))
        plt.plot(time, voltages, marker='o', markersize=2, linestyle='-', linewidth=0.5, color='b', label='UIB1B Data')
        plt.savefig(os.path.join(path, "UIB1B_data_raw.png"))
        plt.close()

def extract_all_files(files):
    for path in files:
        print(f"Processing: {path}")
        datadir = rb.read(path)
        print(datadir)

        ch_files = [f for f in os.listdir(path) if f.endswith('.ch')]
        for ch_file in ch_files:
            print(f"Processing channel file: {ch_file}")
            datafile = datadir.get_file(ch_file)
            time = np.array(datafile.xlabels)
            voltages = np.array(datafile.data)

            outdir = os.path.join(path, "processed")
            os.makedirs(outdir, exist_ok=True)
            np.savetxt(os.path.join(outdir, f"{ch_file.replace('.ch', '_data.txt')}"),
                       np.column_stack((time, voltages)), header="Time (s)\tVoltage (V)", delimiter="\t")
            np.savez_compressed(os.path.join(outdir, f"{ch_file.replace('.ch', '_data.npz')}"), time, voltages)

def extract_multi_channel_data(files, channel_indices):
    for path in files:
        print(path)
        os.chdir(path)
        data_dict = {}
        for idx in channel_indices:
            try:
                result = extract_data(path, index=idx)
                # Normalize y-axis data to np.amax
                y_max = np.amax(result[:, 1])
                if y_max != 0:
                    result[:, 1] = result[:, 1]  # / y_max
                data_dict[f'channel_{idx}'] = result
                # Optionally, save each channel as a separate txt
                output_txt = f"extracted_data_channel_{idx}.txt"
                with open(output_txt, 'w') as f:
                    np.savetxt(f, result, delimiter=',', header='x,y', comments='')
            except Exception as e:
                print(f"Failed to process {path} for channel {idx}: {e}")
        # Save all channels for this file in a single npz
        multioutfile = os.path.join(path, "Agilent" + str(channel_indices) + ".txt")
        np.savez_compressed(multioutfile, **data_dict)









