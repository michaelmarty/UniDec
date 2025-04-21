import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams
import pandas as pd

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 24
rcParams['font.weight'] = 700
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.titleweight'] = 700
rcParams['axes.titlesize'] = 16
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1


areas_df = pd.DataFrame(columns=['File', 'Area1', 'Area2', 'RelArea'])

def datachop(data, min_x, max_x):
    mask = (data[:, 0] >= min_x) & (data[:, 0] <= max_x)
    return data[mask]

def parse_and_plot(topdir, min_x1, max_x1, min_x2, max_x2):
    for root, dirs, files in os.walk(topdir):
        eic_pos_file = None
        eic_neg_file = None

        for file in files:
            if file == 'eic_pos.npz':
                eic_pos_file = os.path.join(root, file)
            elif file == 'eic_neg.npz':
                eic_neg_file = os.path.join(root, file)

        if eic_pos_file and eic_neg_file:
            eic_pos_data = np.load(eic_pos_file)
            eic_neg_data = np.load(eic_neg_file)

            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

            # Plot eic_neg data in blue and calculate area1
            for key in eic_neg_data.files:
                chopped_neg_data = datachop(eic_neg_data[key], min_x1, max_x1)
                ax1.plot(eic_neg_data[key][:, 0], eic_neg_data[key][:, 1], color='blue')
                area1 = np.trapz(chopped_neg_data[:, 1], chopped_neg_data[:, 0])
            ax1.set_title('EIC Negative Mode_DPPG', loc='right')
            ax1.fill_between(eic_neg_data[key][:, 0], eic_neg_data[key][:, 1], 0, where=(eic_neg_data[key][:, 0] > min_x1) & (eic_neg_data[key][:, 0] < max_x1), alpha=0.5, color='blue')
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax1.set_ylabel("Intensity")
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.spines['left'].set_linestyle('-')
            ax1.spines['bottom'].set_linestyle('-')
            ax1.spines['left'].set_visible(True)
            ax1.spines['bottom'].set_visible(True)
            ax1.spines['left'].set_position(('outward', 2))
            ax1.spines['bottom'].set_position(('outward', 2))
            ax1.spines['left'].set_bounds(0, ax1.get_ylim()[1])
            ax1.spines['bottom'].set_bounds(0, len(eic_neg_data[key][:, 0]) - 1)
            ax1.text(1.0, ax1.get_ylim()[1] * 0.85, "Area1: " + str(round(area1, 4)), fontsize=12, ha='right')

            # Plot eic_pos data in red and calculate area2
            for key in eic_pos_data.files:
                chopped_pos_data = datachop(eic_pos_data[key], min_x2, max_x2)
                ax2.plot(eic_pos_data[key][:, 0], eic_pos_data[key][:, 1], color='red')
                area2 = np.trapz(chopped_pos_data[:, 1], chopped_pos_data[:, 0])
            ax2.set_title('EIC Positive Mode_DPPC', loc='right')
            ax2.fill_between(eic_pos_data[key][:, 0], eic_pos_data[key][:, 1], 0, where=(eic_pos_data[key][:, 0] > min_x2) & (eic_pos_data[key][:, 0] < max_x2), alpha=0.5, color='red')
            ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax2.set_ylabel("Intensity")
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['left'].set_linestyle('-')
            ax2.spines['bottom'].set_linestyle('-')
            ax2.spines['left'].set_visible(True)
            ax2.spines['bottom'].set_visible(True)
            ax2.spines['left'].set_position(('outward', 2))
            ax2.spines['bottom'].set_position(('outward', 2))
            ax2.spines['left'].set_bounds(0, ax2.get_ylim()[1])
            ax2.spines['bottom'].set_bounds(0, len(eic_pos_data[key][:, 0]) - 1)
            ax2.text(1.0, ax2.get_ylim()[1] * 0.85, "Area2: " + str(round(area2, 4)), fontsize=12, ha='right')

            # Combined plot and calculate relative area
            for key in eic_neg_data.files:
                ax3.plot(eic_neg_data[key][:, 0], eic_neg_data[key][:, 1], color='blue')
            for key in eic_pos_data.files:
                ax3.plot(eic_pos_data[key][:, 0], eic_pos_data[key][:, 1], color='red')
            relarea = area1 / (area1 + area2) * 100
            ax3.set_title('Combined EIC', loc='right')
            ax3.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax3.set_xlabel("Retention Time (min)", fontsize=12, fontweight=700)
            ax3.set_ylabel("Intensity")
            ax3.spines['top'].set_visible(False)
            ax3.spines['right'].set_visible(False)
            ax3.spines['left'].set_linestyle('-')
            ax3.spines['bottom'].set_linestyle('-')
            ax3.spines['left'].set_visible(True)
            ax3.spines['bottom'].set_visible(True)
            ax3.spines['left'].set_position(('outward', 2))
            ax3.spines['bottom'].set_position(('outward', 2))
            ax3.spines['left'].set_bounds(0, ax3.get_ylim()[1])
            ax3.spines['bottom'].set_bounds(0, len(eic_neg_data[key][:, 0]) - 1)
            ax3.text(1.0, ax3.get_ylim()[1] * 0.85, "Rel. Area: " + str(round(relarea, 4)) + "%", fontsize=12, ha='right')

            # Set the figure title as the name of the subfolder
            subfolder_name = os.path.basename(root)
            fig.suptitle(subfolder_name, fontsize=8, fontweight=700)

            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            os.chdir(root)
            plt.savefig("EIC_Combined.png")
            plt.close()

            areas_df.loc[len(areas_df)] = [os.path.basename(root), area1, area2, relarea]


def save_areas_to_excel(output_path):
    areas_df.to_excel(output_path, index=False)

# Example usage
topdir = "Z:\\Group Share\\Dash\\Synapt Runs 20241129\\DPPG-DPPC\\Pos mode"
parse_and_plot(topdir, 1.40, 1.75, 1.40, 2)
save_areas_to_excel(os.path.join(topdir, "areas.xlsx"))