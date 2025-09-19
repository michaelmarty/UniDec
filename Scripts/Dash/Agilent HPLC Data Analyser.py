import os
import glob
import numpy as np
from LCImporter import extract_data
from LCImporter import parse_uv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import colormaps
from matplotlib import rcParams
from scipy.signal import find_peaks
import pandas as pd
import colorsys

# RC Parameters
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 12
rcParams['font.weight'] = 500
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.titleweight'] = 500
rcParams['axes.titlesize'] = 22
rcParams['axes.linewidth'] = 3
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1


def get_filtered_colors(cmap_name='cool', num_colors=256, luminosity_factor=0.8, saturation_factor=2, hue_shift=0): #A value of 0.5 would reduce the luminosity/sat to half, while 1.0 would keep it unchanged,  A value greater than 1.0 will increase saturation, while 1.0 will keep it unchanged.
    """Returns a list of colors from the colormap with lowered luminosity."""
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

# Define the top directory and file pattern
topdir = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\C4 Tests_Prelim 2025-08-13 15-40-41\Steep_Fast"

files = glob.glob(f'{topdir}/*.D', recursive=False)

# Define the channel index to extract
# FLD=0
# UV=1
channel_index = 2

# Define the channel indices for multi-channel extraction only
channel_indices = [2, 5]

# outfile setup
outfile = "spectra.npz"
multioutfile = "spectra" + str(channel_indices) + ".npz"


# Function to extract data and write to text files
def extract_and_write_data(files):
    for path in files:
        print(path)
        os.chdir(path)
        if True:
            try:
                result = extract_data(path, index=channel_index)
                x, y = result[:, 0], result[:, 1]

                # Write x and y data to a text file in the original directory
                output_txt = "extracted_data.txt"
                with open(output_txt, 'w') as f:
                    np.savetxt(f, np.column_stack((x, y)), delimiter=',', header='x,y', comments='')
                np.savez_compressed(outfile, result)

            except Exception as e:
                print(f"Failed to process {path}: {e}")


def extract_and_write_data_multi_channel(files, channel_indices, multioutfile):
    for path in files:
        print(path)
        os.chdir(path)
        data_dict = {}
        for idx in channel_indices:
            try:
                result = extract_data(path, index=idx)
                # Optionally normalize y-axis data
                # y_max = np.amax(result[:, 1])
                # if y_max != 0:
                #     result[:, 1] = result[:, 1] / y_max
                data_dict[f'channel_{idx}'] = result
                # Save each channel as a separate txt
                output_txt = f"extracted_data_channel_{idx}.txt"
                with open(output_txt, 'w') as f:
                    np.savetxt(f, result, delimiter=',', header='x,y', comments='')
            except Exception as e:
                print(f"Failed to process {path} for channel {idx}: {e}")
        # Save all channels for this file in a single npz
        np.savez_compressed(multioutfile, **data_dict)

def extract_spectra(files, channel_index=0):
    """
    Extracts spectra data from the given files and saves it to a .npz file.
    :param files: List of file paths to extract data from.
    :param channel_index: Index of the channel to extract (0 for FLD, 1 for UV).
    """
    for path in files:
        print(f"Processing {path}...")
        os.chdir(path)
        try:
            result = parse_uv(path, index=channel_index)
            x, y = result[:, 0], result[:, 1]

            # Save the extracted data to a .npz file
            output_npz = os.path.join(path, "Test.npz")
            np.savez_compressed(output_npz, result)

        except Exception as e:
            print(f"Failed to process {path}: {e}")

# Function to read data from npz files
def read_data(files):
    data = []
    file_names = []
    for path in files:
        try:
            output_npz = os.path.join(path, "spectra.npz")
            if os.path.exists(output_npz):
                result = np.load(output_npz)['arr_0']
                data.append(result)
                file_names.append(os.path.basename(path))
        except Exception as e:
            print(f"Failed to read {path}: {e}")
    return data, file_names


# Function to create waterfall plot
def create_waterfall_plot(data, file_names):
    # Sort data and file_names based on file_names
    sorted_indices = sorted(range(len(file_names)), key=lambda i: file_names[i])
    sorted_data = [data[i] for i in sorted_indices]
    sorted_file_names = [file_names[i] for i in sorted_indices]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = get_filtered_colors(num_colors=len(sorted_data))

    for i, d in enumerate(sorted_data):
        x, y = d[:, 0], d[:, 1]
        label = sorted_file_names[i]
        color = colors[i % len(colors)]  # Cycle through the filtered colors
        ax.plot(x, y, zs=i, zdir='y', label=label, color=color)

        # Create vertices for the filled area
        verts = [(x[j], i, y[j]) for j in range(len(x))] + [(x[-1], i, 0), (x[0], i, 0)]
        poly = Poly3DCollection([verts], color=color, alpha=1)
        ax.add_collection3d(poly)

    ax.set_xlabel('Retention Time (min)')
    ax.set_ylabel('Replicate')
    ax.set_zlabel('Intensity')
    ax.set_title("Waterfall Plot of HPLC Data")

    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    plt.savefig(os.path.join(topdir, "WF" + str(channel_index) + ".pdf"), dpi=300)
    plt.savefig(os.path.join(topdir, "WF" + str(channel_index) + ".png"), dpi=300)
    plt.show()


# Function to calculate areas and retention times
def calculate_areas_and_retention_times(data):
    results = []

    for dataset in data:
        # Assuming each dataset is an array with two columns: time (x) and intensity (y)
        x = dataset[:, 0]  # First column for time
        y = dataset[:, 1]  # Second column for intensity
        # print(x)
        # print(y)

        # Find peaks in the y data
        peak_indices, _ = find_peaks(y)
        # print(peak_indices)

        # Define the time threshold for flowthrough and elution peaks
        threshold_time = 4.3  # in minutes

        # Initialize result dictionary for this dataset
        result = {
            'flowthrough_peak': {},
            'elution_peak': {}
        }

        # Identify flowthrough peak (before 4 minutes)
        flowthrough_idx = peak_indices[x[peak_indices] < threshold_time]
        if len(flowthrough_idx) > 0:
            # flowthrough_peak_idx = flowthrough_idx[0]  # Take the first peak before 4 mins
            flowthrough_peak_idx = flowthrough_idx[np.argmax(y[flowthrough_idx])]
            # print(flowthrough_idx)
            # print(flowthrough_peak_idx)
            flowthrough_array = np.array([[x[flowthrough_peak_idx]], [y[flowthrough_peak_idx]]])
            print(flowthrough_array)

            x_lower_bound = (x >= 1) & (x < flowthrough_array[0])
            selected_y = y[x_lower_bound]
            selected_x = x[x_lower_bound]
            lower_indices = np.argmin(selected_y)
            lower_bound_time = selected_x[lower_indices]
            print("Lower Bound Time: ", lower_bound_time)
            # lower_bound_intensity=selected_y[lower_indices]
            # print(lowe_bound_time)
            # print(lowe_bound_intensity)

            x_upper_bound = (x >= flowthrough_array[0]) & (x < threshold_time)
            selected_up_y = y[x_upper_bound]
            selected_up_x = x[x_upper_bound]
            upper_indices = np.argmin(selected_up_y)
            upper_bound_time = selected_up_x[upper_indices]
            print("Upper Bound Time: ", upper_bound_time)
            # upper_bound_intensity=selected_up_y[upper_indices]

            x_sub = (x >= lower_bound_time) & (x <= upper_bound_time)
            # y_sub=(y>=lowe_bound_intensity)&(y<=upper_bound_intensity)
            x_data = x[x_sub]
            y_data = y[x_sub]
            # print(y_data)
            # print(x_data.shape)
            # print(y_data.shape)

            # Calculate the area under the curve for the flowthrough peak
            area_flowthrough = np.trapz(y_data, x_data)
            print(area_flowthrough)

            # Store results
            result['flowthrough_peak'] = {
                'data': flowthrough_array,
                'flowthrough_retention_time': x[flowthrough_peak_idx],
                'lower_bound_time': lower_bound_time,
                'upper_bound_time': upper_bound_time,
                'flowthrough_area': area_flowthrough
            }

        # Identify elution peaks (after 4 minutes)
        elution_idx = peak_indices[x[peak_indices] > threshold_time]
        if len(elution_idx) > 0:
            elution_peak_idx = elution_idx[np.argmax(y[elution_idx])]
            # print(elution_idx)
            # print(elution_peak_idx)
            elution_array = np.array([[x[elution_peak_idx]], [y[elution_peak_idx]]])
            print(elution_array)

            x_el_lower_bound = (x >= threshold_time) & (x < elution_array[0])
            selected_el_y = y[x_el_lower_bound]
            selected_el_x = x[x_el_lower_bound]
            lower_el_indices = np.argmin(selected_el_y)
            lower_bound_el_time = selected_el_x[lower_el_indices]
            print("Lower Bound Elution Time: ", lower_bound_el_time)
            # lower_bound_intensity=selected_y[lower_indices]
            # print(lowe_bound_time)
            # print(lowe_bound_intensity)

            x_el_upper_bound = (x >= elution_array[0]) & (x < 7)
            selected_el_up_y = y[x_el_upper_bound]
            selected_el_up_x = x[x_el_upper_bound]
            upper_el_indices = np.argmin(selected_el_up_y)
            upper_bound_el_time = selected_el_up_x[upper_el_indices]
            print("Upper Bound Elution Time: ", upper_bound_el_time)
            # upper_bound_intensity=selected_up_y[upper_indices]

            x_el_sub = (x >= lower_bound_el_time) & (x <= upper_bound_el_time)
            # y_sub=(y>=lowe_bound_intensity)&(y<=upper_bound_intensity)
            x_el_data = x[x_el_sub]
            y_el_data = y[x_el_sub]
            # print(y_data)
            # print(x_data.shape)
            # print(y_data.shape)

            # Calculate the area under the curve for the flowthrough peak
            area_elution = np.trapz(y_el_data, x_el_data)
            print(area_elution)

            # Store results
            result['elution_peak'] = {
                'data': elution_array,
                'elution_retention_time': x[elution_peak_idx],
                'lower_el_bound_time': lower_bound_el_time,
                'upper_el_bound_time': upper_bound_el_time,
                'max_el': y[elution_peak_idx],
                'elution_area': area_elution
            }

        # Append the result for this dataset to results list
        results.append(result)

    # export results to a csv file
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(topdir, "results.csv"), index=False)

    return results


def create_annotated_ridgeline_plot(data, file_names, results, offset=100):
    # Sort data, results, and file_names based on file_names
    sorted_indices = sorted(range(len(file_names)), key=lambda i: file_names[i])
    sorted_data = [data[i] for i in sorted_indices]
    sorted_results = [results[i] for i in sorted_indices]
    sorted_file_names = [file_names[i] for i in sorted_indices]

    plt.figure(figsize=(12, 12))
    colors = get_filtered_colors(num_colors=len(sorted_data))


    for i, (d, result) in enumerate(zip(sorted_data, sorted_results)):
        x, y = d[:, 0], d[:, 1]
        color = colors[i % len(colors)]
        plt.fill_between(x, y + i * offset, i * offset, alpha=0.4, color=color)
        plt.plot(x, y + i * offset, color=color, linewidth=0.5)
        plt.text(0, i * offset + 2, sorted_file_names[i], fontsize=7, fontweight='bold', ha='left', color=color)

        # Annotate flowthrough if flowthrough_peak exists
        if 'flowthrough_peak' in result and 'flowthrough_retention_time' in result['flowthrough_peak'] and 'flowthrough_area' in result['flowthrough_peak']:
            plt.annotate(
                f"Flowthrough RT: {result['flowthrough_peak']['flowthrough_retention_time']:.2f}\nArea: {result['flowthrough_peak']['flowthrough_area']:.2f}",
                xy=(result['flowthrough_peak']['flowthrough_retention_time'], result['flowthrough_peak']['data'][1] + i * offset + 2),
                xytext=(result['flowthrough_peak']['flowthrough_retention_time'] + 0.5, result['flowthrough_peak']['data'][1] + i * offset + 2),
                #arrowprops=dict(facecolor= 'black', shrink=0.05, width=0.5, headwidth=4, connectionstyle="arc3,rad=0.3"),
                fontsize=7, fontweight='bold',
                color=color
            )

        # Annotate elution if elution_retention_time and elution_area exist
        if 'elution_peak' in result and 'elution_retention_time' in result['elution_peak'] and 'elution_area' in result['elution_peak']:
            plt.annotate(
                f"Elution RT: {result['elution_peak']['elution_retention_time']:.2f}\nArea: {result['elution_peak']['elution_area']:.2f}",
                xy=(result['elution_peak']['elution_retention_time'], min(d[:,1]) + i * offset + 2),
                xytext=(result['elution_peak']['elution_retention_time'] + 0.5, min(d[:,1]) + i * offset + 2),
                #arrowprops=dict(facecolor= 'black', shrink=0.05, width=0.5, headwidth=4, connectionstyle="arc3,rad=0.3"),
                fontsize=7, fontweight='bold',
                color=color
            )
        else:
            print(f"Warning: Missing elution data for result #{i}. Keys must include 'elution_retention_time' and 'elution_area'.")

    plt.xlabel("Time (min)", fontweight='bold')
    plt.ylabel("")  # Remove y-axis label
    plt.title("Annotated Ridgeline Plot of HPLC Data", fontweight='bold')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['bottom'].set_bounds(0, 7)
    plt.subplots_adjust(hspace = -0.3)

    plt.yticks([])  # Remove y-axis ticks

    plt.savefig(os.path.join(topdir, "ARP" + str(channel_index) + ".pdf"), dpi=300)
    plt.savefig(os.path.join(topdir, "ARP" + str(channel_index) + ".png"), dpi=300)
    plt.show()

#Use data from the extracted data to create a 2D plot
def create_2D_plot(data, file_names):
    # Sort data and file_names based on file_names
    sorted_indices = sorted(range(len(file_names)), key=lambda i: file_names[i])
    sorted_data = [data[i] for i in sorted_indices]
    sorted_file_names = [file_names[i] for i in sorted_indices]

    plt.figure(figsize=(16, 7))
    colors = get_filtered_colors(num_colors=len(sorted_data))

    for i, d in enumerate(sorted_data):
        x, y = d[:, 0], d[:, 1]
        y_sub = y - np.min(y)  # Subtract the minimum to avoid negative values
        label = sorted_file_names[i]
        color = colors[i % len(colors)]  # Cycle through the filtered colors
        plt.plot(x, y_sub, label=label, color=color, linewidth=1)
        #plt.fill_between(x, y, label=label, color=color, alpha=0.3)  # Fill the area under the curve

    plt.xlabel('Retention Time (min)', fontweight='bold')
    #plt.xlim(0, 10) # Set x-axis limits
    plt.ylabel('Intensity [LU]', fontweight='bold')
    #plt.ylim(0, 180)
    plt.legend()
    #plt.title("Peak Overlay", fontweight='bold')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    # Overlay timetable as a step plot on a secondary y-axis
    # t1 = np.array([0, 4])
    # p1 = np.array([0, 0])
    #
    # t2 = np.linspace(4, 4.01, 20)  # Slope segment
    # p2 = np.linspace(0, 100, 20)
    #
    # t3 = np.array([4.01, 5.49])
    # p3 = np.array([100, 100])
    #
    # t4 = np.linspace(5.49, 5.5, 20)  # Slope down
    # p4 = np.linspace(100, 0, 20)
    #
    # t5 = np.array([5.5, 8.10])
    # p5 = np.array([0, 0])
    #
    # t6 = np.linspace(8.10, 8.11, 10)  # Slope up
    # p6 = np.linspace(0, 100, 10)
    #
    # t7 = np.array([8.11, 9.59])
    # p7 = np.array([100, 100])
    #
    # t8 = np.linspace(9.59, 9.6, 10)  # Slope down
    # p8 = np.linspace(100, 0, 10)
    #
    # t9 = np.array([9.6, 10])
    # p9 = np.array([0, 0])
    #
    # # Concatenate all segments
    # timetable_times = np.concatenate([t1, t2[1:], t3[1:], t4[1:], t5[1:], t6[1:], t7[1:], t8[1:], t9[1:]])
    # timetable_percent = np.concatenate([p1, p2[1:], p3[1:], p4[1:], p5[1:], p6[1:], p7[1:], p8[1:], p9[1:]])
    #
    # # In your plotting function, replace the timetable overlay with:
    # ax2 = plt.gca().twinx()
    # ax2.plot(timetable_times, timetable_percent, color='black', linewidth=1, label='Timetable (%)', alpha=0.5)
    # ax2.set_ylabel('B (%)', fontweight='bold', color='black')
    # ax2.tick_params(axis='y', labelcolor='black')
    # ax2.set_ylim(-5, 105)
    # # x2.legend(loc='upper right')
    # plt.gca().spines['top'].set_visible(False)

    #plt.legend()
    plt.savefig(os.path.join(topdir, "2D" + str(channel_index) + ".pdf"), dpi=300)
    plt.savefig(os.path.join(topdir, "2D" + str(channel_index) + ".png"), dpi=300)
    plt.show()


def calculate_average_and_rsd(areas):
    average = np.mean(areas)
    std_dev = np.std(areas)
    rsd = (std_dev / average) * 100 if average != 0 else 0
    return average, rsd

def calculate_relative_elution_area_and_rsd(results):
    relative_elution_areas = []
    for result in results:
        if 'elution_peak' in result and 'flowthrough_peak' in result:
            elution_area = result['elution_peak']['elution_area']
            flowthrough_area = result['flowthrough_peak']['flowthrough_area']
            if elution_area + flowthrough_area != 0:
                relative_elution_area = ((elution_area / (elution_area + flowthrough_area))*100)
                relative_elution_areas.append(relative_elution_area)

    average_relative_elution_area = np.mean(relative_areas)
    std_dev_relative_elution_area = np.std(relative_areas)
    rsd_relative_elution_area = (std_dev_relative_elution_area / average_relative_elution_area) * 100 if average_relative_area != 0 else 0

    return average_relative_elution_area, rsd_relative_elution_area

#function to calculate relative flowthrough areas and RSD
def calculate_relative_flowthrough_area_and_rsd(results):
    relative_flowthrough_areas = []
    for result in results:
        if 'elution_peak' in result and 'flowthrough_peak' in result:
            elution_area = result['elution_peak']['elution_area']
            flowthrough_area = result['flowthrough_peak']['flowthrough_area']
            if elution_area + flowthrough_area != 0:
                relative_flowthrough_area = ((flowthrough_area / (elution_area + flowthrough_area))*100)
                relative_flowthrough_areas.append(relative_flowthrough_area)

    average_relative_flowthrough_area = np.mean(relative_flowthrough_areas)
    std_dev_relative_flowthrough_area = np.std(relative_flowthrough_areas)
    rsd_relative_flowthrough_area = (std_dev_relative_flowthrough_area / average_relative_flowthrough_area) * 100 if average_relative_flowthrough_area != 0 else 0

    return average_relative_flowthrough_area, rsd_relative_flowthrough_area


# Function to create flowthrough area scatter plot
def create_flowthrough_area_plot(results, file_names):
    # Extract flowthrough areas and sort file names
    flowthrough_areas = [result['flowthrough_peak']['flowthrough_area'] for result in results if
                         'flowthrough_peak' in result]
    sorted_indices = sorted(range(len(file_names)), key=lambda i: file_names[i])
    sorted_file_names = [file_names[i] for i in sorted_indices]
    sorted_flowthrough_areas = [flowthrough_areas[i] for i in sorted_indices if i < len(flowthrough_areas)]

    # Calculate average and %RSD
    average_flowthrough, rsd_flowthrough = calculate_average_and_rsd(sorted_flowthrough_areas)

    color = get_filtered_colors(num_colors=len(sorted_flowthrough_areas))

    # for i in range(len(sorted_flowthrough_areas)):
    #     color = colors[i % len (colors)]


    plt.figure(figsize=(12, 6))
    plt.scatter(sorted_file_names, sorted_flowthrough_areas, color = color, alpha=0.7, label='Flowthrough Area')
    plt.xlabel('File Names', fontweight='bold')
    plt.ylabel('Flowthrough Area', fontweight='bold')
    plt.title('Flowthrough Areas Tracking', fontweight='bold')

    # Annotate average and %RSD
    plt.annotate(f'Avg: {average_flowthrough:.2f}\n%RSD: {rsd_flowthrough:.2f}',
                 xy=(0.99, 0.99), xycoords='axes fraction',
                 fontsize=18, ha='right', fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))

    plt.xticks(rotation=90)
    plt.tight_layout()
    # plt.grid()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(os.path.join(topdir, "Flowthrough_Area_Scatter_Plot.pdf"), dpi=300)
    plt.savefig(os.path.join(topdir, "Flowthrough_Area_Scatter_Plot.png"), dpi=300)
    plt.show()


# Function to create elution area scatter plot
def create_elution_area_plot(results, file_names):
    # Extract elution areas and sort file names
    elution_areas = [result['elution_peak']['elution_area'] for result in results if 'elution_peak' in result]
    sorted_indices = sorted(range(len(file_names)), key=lambda i: file_names[i])
    sorted_file_names = [file_names[i] for i in sorted_indices]
    sorted_elution_areas = [elution_areas[i] for i in sorted_indices if i < len(elution_areas)]

    # Calculate average and %RSD
    average_elution, rsd_elution = calculate_average_and_rsd(sorted_elution_areas)

    color = get_filtered_colors(num_colors=len(sorted_elution_areas))

    # for i in range(len(sorted_elution_areas)):
    #     color = colors[i % len(colors)]

    plt.figure(figsize=(12, 6))
    plt.scatter(sorted_file_names, sorted_elution_areas, color= color, alpha=0.7, label='Elution Area')
    plt.xlabel('File Names', fontweight='bold')
    plt.ylabel('Elution Area', fontweight='bold')
    plt.title('Elution Areas Tracking', fontweight='bold')

    # Annotate average and %RSD
    plt.annotate(f'Avg: {average_elution:.2f}\n%RSD: {rsd_elution:.2f}',
                 xy=(0.99, 0.99), xycoords='axes fraction',
                 fontsize=18, ha='right', fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))

    plt.xticks(rotation=90)
    plt.tight_layout()
    #plt.grid()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(os.path.join(topdir, "Elution_Area_Scatter_Plot.pdf"), dpi=300)
    plt.savefig(os.path.join(topdir, "Elution_Area_Scatter_Plot.png"), dpi=300)
    plt.show()

#function to create relative elution area scatter plot
def create_relative_elution_area_plot(results, file_names):
    # Extract relative elution areas and sort file names
    relative_areas = []
    for result in results:
        if 'elution_peak' in result and 'flowthrough_peak' in result:
            elution_area = result['elution_peak']['elution_area']
            flowthrough_area = result['flowthrough_peak']['flowthrough_area']
            if elution_area + flowthrough_area != 0:
                relative_elution_area = (elution_area / (elution_area + flowthrough_area))*100
                relative_areas.append(relative_elution_area)

    sorted_indices = sorted(range(len(file_names)), key=lambda i: file_names[i])
    sorted_file_names = [file_names[i] for i in sorted_indices]
    sorted_relative_areas = [relative_areas[i] for i in sorted_indices if i < len(relative_areas)]

    # Calculate average and %RSD
    average_relative_area, rsd_relative_area = calculate_average_and_rsd(sorted_relative_areas)

    color = get_filtered_colors(num_colors=len(sorted_relative_areas))

    # for i in range(len(sorted_relative_areas)):
    #     color = colors[i % len(colors)]

    plt.figure(figsize=(12, 6))
    plt.scatter(sorted_file_names, sorted_relative_areas, color=color, alpha=0.7, label='Relative Elution Area')
    plt.xlabel('File Names', fontweight='bold')
    plt.ylabel('Relative Elution Area', fontweight='bold')
    plt.title('Relative Elution Areas Tracking', fontweight='bold')

    # Annotate average and %RSD
    plt.annotate(f'Avg: {average_relative_area:.2f}\n%RSD: {rsd_relative_area:.2f}',
                 xy=(0.99, 0.99), xycoords='axes fraction',
                 fontsize=18, ha='right', fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))

    plt.xticks(rotation=90)
    plt.tight_layout()
    # plt.grid()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(os.path.join(topdir, "Relative_Elution_Area_Scatter_Plot.pdf"), dpi=300)
    plt.savefig(os.path.join(topdir, "Relative_Elution_Area_Scatter_Plot.png"), dpi=300)
    plt.show()

# Function to create relative flowthrough area scatter plot
def create_relative_flowthrough_area_plot(results, file_names):
    # Extract relative flowthrough areas and sort file names
    relative_areas = []
    for result in results:
        if 'elution_peak' in result and 'flowthrough_peak' in result:
            elution_area = result['elution_peak']['elution_area']
            flowthrough_area = result['flowthrough_peak']['flowthrough_area']
            if elution_area + flowthrough_area != 0:
                relative_flowthrough_area = (flowthrough_area / (elution_area + flowthrough_area)*100)
                relative_areas.append(relative_flowthrough_area)

    sorted_indices = sorted(range(len(file_names)), key=lambda i: file_names[i])
    sorted_file_names = [file_names[i] for i in sorted_indices]
    sorted_relative_areas = [relative_areas[i] for i in sorted_indices if i < len(relative_areas)]

    # Calculate average and %RSD
    average_relative_area, rsd_relative_area = calculate_average_and_rsd(sorted_relative_areas)

    color = get_filtered_colors(num_colors=len(sorted_relative_areas))

    # for i in range(len(sorted_relative_areas)):
    #     color = colors[i % len(colors)]

    plt.figure(figsize=(12, 6))
    plt.scatter(sorted_file_names, sorted_relative_areas, color=color, alpha=0.7, label='Relative Flowthrough Area')
    plt.xlabel('File Names', fontweight='bold')
    plt.ylabel('Relative Flowthrough Area', fontweight='bold')
    plt.title('Relative Flowthrough Areas Tracking', fontweight='bold')

    # Annotate average and %RSD
    plt.annotate(f'Avg: {average_relative_area:.2f}\n%RSD: {rsd_relative_area:.2f}',
                 xy=(0.99, 0.99), xycoords='axes fraction',
                 fontsize=18, ha='right', fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))

    plt.xticks(rotation=90)
    plt.tight_layout()
    # plt.grid()
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig(os.path.join(topdir, "Relative_Flowthrough_Area_Scatter_Plot.pdf"), dpi=300)
    plt.savefig(os.path.join(topdir, "Relative_Flowthrough_Area_Scatter_Plot.png"), dpi=300)
    plt.show()

def get_timetable():
    t1 = np.array([0, 1])
    p1 = np.array([10,10])
    t2 = np.linspace(1, 1.01, 20)
    p2 = np.linspace(10, 100, 20)
    t3 = np.array([1.01, 5.0])
    p3 = np.array([100, 100])
    t4 = np.linspace(5.0, 6.0, 20)
    p4 = np.linspace(100, 10, 20)
    t5 = np.array([6.0, 7.0])
    p5 = np.array([10, 10])
    # t6 = np.linspace(8.10, 8.6, 10)
    # p6 = np.linspace(0, 100, 10)
    # t7 = np.array([8.6, 9.1])
    # p7 = np.array([100, 100])
    # t8 = np.linspace(9.1, 9.6, 10)
    # p8 = np.linspace(100, 0, 10)
    # t9 = np.array([9.6, 10])
    # p9 = np.array([0, 0])
    timetable_times = np.concatenate([t1, t2[1:], t3[1:], t4[1:], t5[1:]]) #t6[1:], t7[1:], t8[1:], t9[1:]])
    timetable_percent = np.concatenate([p1, p2[1:], p3[1:], p4[1:], p5[1:]]) #p6[1:], p7[1:], p8[1:], p9[1:]])
    return timetable_times, timetable_percent


def create_multi_channel_npz(topdir, channel_indices, title=""):
    """
    Reads all multi-channel .npz files in topdir and plots each channel as a 2D overlay,
    and also shows each channel as a separate subplot for each file.
    """
    # Find all .npz files in topdir that match the multi-channel pattern
    npz_files = [
        os.path.join(root, file)
        for root, _, files in os.walk(topdir)
        for file in files
        if file.startswith("spectra") and file.endswith(".npz")
    ]
    if not npz_files:
        print(f"No multi-channel .npz files found in {topdir}")
        return

    # Timetable definition (Gradient Information for %B)


    timetable_times, timetable_percent = get_timetable()
    colors = get_filtered_colors(num_colors=len(channel_indices))

    # Overlay plot for all files and channels
    plt.figure(figsize=(16, 9))
    for multioutfile in npz_files:
        data = np.load(multioutfile)
        for i, idx in enumerate(channel_indices):
            key = f'channel_{idx}'
            if key in data:
                arr = data[key]
                x, y = arr[:, 0], arr[:, 1]
                if idx == 2:
                    y = y * 25.0
                label = f'{os.path.basename(multioutfile)} - Channel {idx}'
                plt.plot(x, y, label=label, color=colors[i % len(colors)], alpha=0.7)
            else:
                print(f"Channel {idx} not found in {multioutfile}")

    plt.xlabel('Retention Time (min)', fontweight='bold')
    plt.ylabel('Intensity', fontweight='bold')
    plt.title("Peak Overlay", fontweight='bold')
    plt.legend(loc='upper right', fontsize=4)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    ax2 = plt.gca().twinx()
    ax2.plot(timetable_times, timetable_percent, color='black', linewidth=1, label='Timetable (%)', alpha=0.5)
    ax2.set_ylabel('B (%)', fontweight='bold', color='black')
    ax2.tick_params(axis='y', labelcolor='black')
    ax2.set_ylim(-5, 105)
    plt.tight_layout()
    plt.savefig(os.path.join(topdir, "Multi_Channel_Overlay.pdf"), dpi=300)
    plt.savefig(os.path.join(topdir, "Multi_Channel_Overlay.png"), dpi=300)
    plt.show()

    # Individual subplots for each file, one subplot per channel
    for multioutfile in npz_files:
        data = np.load(multioutfile)
        num_channels = len(channel_indices)
        fig, axes = plt.subplots(num_channels + 1, 1, figsize=(16, 4 * (num_channels + 1)), sharex=True)
        if num_channels == 1:
            axes = [axes]
        for i, idx in enumerate(channel_indices):
            key = f'channel_{idx}'
            ax = axes[i]
            if key in data:
                arr = data[key]
                x, y = arr[:, 0], arr[:, 1]
                if idx == 2:
                    y = y * 25.0
                ax.plot(x, y, color=colors[i % len(colors)], alpha=0.7)
                ax.set_ylabel(f'Intensity (Channel {idx})', fontweight='bold')
            else:
                ax.set_ylabel(f'Channel {idx} not found', fontweight='bold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        timetable_ax = axes[-1]
        timetable_ax.plot(timetable_times, timetable_percent, color='black', linewidth=2, label='Timetable (%)', alpha=0.7)
        timetable_ax.set_ylabel('B (%)', fontweight='bold', color='black')
        timetable_ax.set_ylim(-5, 105)
        timetable_ax.spines['top'].set_visible(False)
        timetable_ax.spines['right'].set_visible(False)
        timetable_ax.set_xlabel('Time (min)', fontweight='bold')
        #timetable_ax.set_title('Timetable', fontweight='bold')
        #fig.suptitle(f"Channels for {os.path.basename(multioutfile)}", fontweight='bold')
        plt.savefig(os.path.join(topdir, f"Multi_Channel_{os.path.basename(multioutfile)}.pdf"), dpi=300)
        plt.savefig(os.path.join(topdir, f"Multi_Channel_{os.path.basename(multioutfile)}.png"), dpi=300)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.close()

# Main function to parse files and process data
def main():
    extract_and_write_data(files)
    extract_and_write_data_multi_channel(files, channel_indices, multioutfile)
    #extract_spectra(files, channel_index=channel_index)
    data, file_names = read_data(files)
    get_timetable()
    #create_waterfall_plot(data, file_names)
    #results = calculate_areas_and_retention_times(data)
    #create_annotated_ridgeline_plot(data, file_names, results, offset=100)
    create_2D_plot(data, file_names)
    # create_flowthrough_area_plot(results, file_names)
    # create_elution_area_plot(results, file_names)
    create_multi_channel_npz(topdir, channel_indices, title=" ")
    #create_relative_elution_area_plot(results, file_names)
    #create_relative_flowthrough_area_plot(results, file_names)



if __name__ == "__main__":
    main()


