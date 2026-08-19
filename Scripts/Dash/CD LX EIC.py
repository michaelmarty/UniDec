import os
import glob
import unidec.tools as ud
from UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import pandas as pd


# RC Params
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 12
rcParams['font.weight'] = 200
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.titleweight'] = 700
rcParams['axes.titlesize'] = 16
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1

topdir = "Z:\\Group Share\\Dash\\Current Data\\Synapt 96 Wells\\Samples"
files = glob.glob(f'{topdir}/*.raw', recursive=False)

output = []
print(files)

mass1 = 426.2875
disc1 = "contaminant"
mass2 = 734.5621
disc2 = "DPPC"

ppm = 100
mtol1 = ppm * mass1 / 1e6
mtol2 = ppm * mass2 / 1e6

print(mtol1, mtol2)

mtol = 0.002

outfile = "eic.npz"

# DataFrame to store area1, area2, relarea1 and file name
areas_df = pd.DataFrame(columns=['File', 'Area1', 'Area2', 'RelArea1', 'RelArea2'])

def process_files(files):
    for path in files[:]:
        print(path)
        os.chdir(path)
        if not os.path.isfile(outfile):
            try:
                d = ImporterFactory.create_importer(path)
            except Exception as e:
                print("File Failed", path)
                print(e)
                continue
            tic1 = d.get_eic(mass1, mtol)
            tic2 = d.get_eic(mass2, mtol)

            tic1_peak_time = tic1[np.argmax(tic1[:, 1]), 0]
            tic2_peak_time = tic2[np.argmax(tic2[:, 1]), 0]

            expected_time1 = tic1_peak_time
            expected_time2 = tic2_peak_time

            if abs(tic1_peak_time - expected_time1) > 0.2 or abs(tic2_peak_time - expected_time2) > 0.2:
                print("Peak Time1:", tic1_peak_time, "Peak Time2:", tic2_peak_time)
                raise Exception("Peak Time Error")

            else:
                sd1 = 0.6
                t1_range = [tic1_peak_time - sd1, tic1_peak_time + sd1]
                sd2 = 0.6
                t2_range = [tic2_peak_time - sd2, tic2_peak_time + sd2]

            print("Peak Time1:", t1_range, "Peak Time2:", t2_range)

            area1 = ud.integrate(tic1, t1_range[0], t1_range[1])[0]
            area2 = ud.integrate(tic2, t2_range[0], t2_range[1])[0]
            if area1 > 0 and area2 > 0 and (area1 / area2 >= 100 or area2 / area1 >= 100):
                if area1 < area2:
                    area1 = 0
                else:
                 area2 = 0
            relarea1 = (area1 / (area1 + area2)) * 100
            relarea2 = (area2 / (area1 + area2)) * 100
            print("Area1:", area1, "Area2:", area2, "RelArea1:", relarea1, "RelArea2:", relarea2)

            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

            ax1.plot(tic1[:, 0], tic1[:, 1], color='blue')
            ax1.set_title('mass = ' + str(mass1), loc='right')
            ax1.fill_between(tic1[:, 0], tic1[:, 1], 0, where=(tic1[:, 0] > t1_range[0]) & (tic1[:, 0] < t1_range[1]), alpha=0.5, color='blue')
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax1.set_ylabel("Intensity")
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)
            ax1.text(1.0, ax1.get_ylim()[1] * 0.86, "Area1: " + str(round(area1, 4)), fontsize=12, ha='left', va='top')
            ax1.spines['left'].set_linestyle('-')
            ax1.spines['bottom'].set_linestyle('-')
            ax1.spines['left'].set_visible(True)
            ax1.spines['bottom'].set_visible(True)
            ax1.spines['left'].set_position(('outward', 2))
            ax1.spines['bottom'].set_position(('outward', 2))
            ax1.spines['left'].set_bounds(0, ax1.get_ylim()[1])
            ax1.spines['bottom'].set_bounds(0, len(tic1[:, 0]) - 1)

            ax2.plot(tic2[:, 0], tic2[:, 1], color='red')
            ax2.set_title('mass = ' + str(mass2), loc='right')
            ax2.fill_between(tic2[:, 0], tic2[:, 1], 0, where=(tic2[:, 0] > t2_range[0]) & (tic2[:, 0] < t2_range[1]), alpha=0.5, color='red')
            ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            ax2.set_ylabel("Intensity")
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.text(1.0, ax2.get_ylim()[1] * 0.86, "Area2: " + str(round(area2, 4)), fontsize=12, ha='left', va='top')
            ax2.spines['left'].set_linestyle('-')
            ax2.spines['bottom'].set_linestyle('-')
            ax2.spines['left'].set_visible(True)
            ax2.spines['bottom'].set_visible(True)
            ax2.spines['left'].set_position(('outward', 2))
            ax2.spines['bottom'].set_position(('outward', 2))
            ax2.spines['left'].set_bounds(0, ax2.get_ylim()[1])
            ax2.spines['bottom'].set_bounds(0, len(tic1[:, 0]) - 1)

            ax3.plot(tic1[:, 0], tic1[:, 1], color='blue')
            ax3.set_title('Combined', loc='right')
            ax3.plot(tic2[:, 0], tic2[:, 1], color='red')
            ax3.annotate("DMPC Rel. Area: " + str(round(relarea1, 2)) + "%", xy=(1.0, ax3.get_ylim()[1] * 0.86),
             xytext=(10, 10), textcoords='offset points', fontsize=12, ha='left', va='top')
            ax3.annotate("DPPC Rel. Area: " + str(round(relarea2, 2)) + "%", xy=(1.0, ax3.get_ylim()[1] * 0.83),
             xytext=(10, -10), textcoords='offset points', fontsize=12, ha='left', va='top')
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
            ax3.spines['bottom'].set_bounds(0, len(tic1[:, 0]) - 1)

            folder_name = os.path.basename(os.path.dirname(outfile))
            fig.suptitle(folder_name, fontsize=16, fontweight='bold')

            plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
            plt.savefig(os.path.join(path, ("Chrom" + str(path.split("\\")[-1]) + ".png")))
            plt.close()
            np.savez_compressed(outfile, array1=tic1, array2=tic2)

            # Append areas to DataFrame
            areas_df.loc[len(areas_df)] = [os.path.basename(path), area1, area2, relarea1, relarea2]




def read_npz_files(topdir):
    areas_df = pd.DataFrame(columns=['File', 'Area1', 'Area2', 'RelArea1', 'RelArea2'])
    for root, dirs, files in os.walk(topdir):
        for file in files:
            if file == 'eic2.npz':
                npz_file_path = os.path.join(root, file)
                data = np.load(npz_file_path)
                tic1 = data['array1']
                tic2 = data['array2']

                # Calculate peak times
                tic1_peak_time = tic1[np.argmax(tic1[:, 1]), 0]
                tic2_peak_time = tic2[np.argmax(tic2[:, 1]), 0]

                # Define time ranges
                sd1 = 0.6
                t1_range = [tic1_peak_time - sd1, tic1_peak_time + sd1]
                sd2 = 0.6
                t2_range = [tic2_peak_time - sd2, tic2_peak_time + sd2]

                # Calculate areas using ud.integrate
                area1 = ud.integrate(tic1, t1_range[0], t1_range[1])[0]
                area2 = ud.integrate(tic2, t2_range[0], t2_range[1])[0]

                # Adjust areas if necessary
                if area1 > 0 and area2 > 0 and (area1 / area2 >= 100 or area2 / area1 >= 100):
                    if area1 < area2:
                        area1 = 0
                    else:
                        area2 = 0

                # Calculate relative areas
                relarea1 = (area1 / (area1 + area2)) * 100 if (area1 + area2) != 0 else 0
                relarea2 = (area2 / (area1 + area2)) * 100 if (area1 + area2) != 0 else 0

                # Append to DataFrame
                areas_df.loc[len(areas_df)] = [os.path.basename(root), area1, area2, relarea1, relarea2]

    return areas_df

def plot_ridgeline(df, topdir):
    fig, ax = plt.subplots(figsize=(16, 12))
    offset = 5
    df = df.sort_values(by='File')

    for index, row in df.iterrows():
        npz_file_path = os.path.join(topdir, row['File'], 'eic.npz')

        if not os.path.exists(npz_file_path):
            print(f"File not found: {npz_file_path}")
            continue

        data = np.load(npz_file_path)
        tic1 = data['array1']
        tic2 = data['array2']

        ax.fill_between(tic1[:, 0], tic1[:, 1] + offset, offset, alpha=0.5, color='blue', label=disc1 if index == 0 else "")
        ax.plot(tic1[:, 0], tic1[:, 1]/np.amax(tic2[:,1]*100) + offset, color='blue', linewidth=0.5)
        ax.fill_between(tic2[:, 0], tic2[:, 1] + offset, offset, alpha=0.5, color='red', label=disc2 if index == 0 else "")
        ax.plot(tic2[:, 0], tic2[:, 1]/np.amax(tic2[:,1]*100) + offset, color='red', linewidth=0.5)

        offset += max(tic1[:, 1].max(), tic2[:, 1].max()) * 1.1

    ax.set_xlabel("Retention Time (min)", fontsize=12, fontweight=700)
    ax.set_ylabel("Intensity")
    ax.set_title("Ridgeline Plot of TIC1 and TIC2", fontsize=16, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_position(('outward', 2))
    ax.spines['bottom'].set_bounds(0, len(df) - 1)
    ax.set_yticks([])  # Remove y-axis ticks
    ax.legend()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.savefig(os.path.join(topdir, "Ridgeline_Plot.png"))
    plt.show()


def save_to_excel(df, output_path):
    df.to_excel(output_path, index=False)
    print(f"Data saved to {output_path}")

def plot_areas(df):
    plt.figure(figsize=(12, 8))
    df = df.sort_values(by='File')
    df['File'] = df['File'].apply(lambda x: x.split("_")[-2] if len(x.split("_")) > 1 else x)
    plt.plot(df['File'], df['Area1'], label=disc1, marker='o', linestyle='--', color='#008080')
    plt.plot(df['File'], df['Area2'], label=disc2, marker='o', linestyle='--', color='red')
    for i, txt in enumerate(df['Area1']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['Area1'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='#008080', fontsize=8)
    for i, txt in enumerate(df['Area2']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['Area2'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='red', fontsize=8)
    plt.xlabel('Fraction')
    plt.ylabel('Area [Raw]')
    plt.title('Area1 and Area2 vs Fraction')
    plt.legend()
    plt.xticks(rotation=45)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_linestyle('-')
    plt.gca().spines['bottom'].set_linestyle('-')
    plt.gca().spines['left'].set_visible(True)
    plt.gca().spines['bottom'].set_visible(True)
    plt.gca().spines['left'].set_position(('outward', 2))
    plt.gca().spines['bottom'].set_position(('outward', 2))
    #plt.gca().spines['left'].set_bounds(df['Area1'].min(), df['Area1'].max())
    plt.gca().spines['bottom'].set_bounds(0, len(df['File']) - 1)
    plt.tight_layout()
    plt.savefig(os.path.join(topdir, "Area vs File[raw].pdf"))
    plt.show()

def plot_areas_percent(df):
    plt.figure(figsize=(12, 8))
    df = df.sort_values(by='File')
    df['File'] = df['File'].apply(lambda x: x.split("_")[-2] if len(x.split("_")) > 1 else x)
    min_area = min(df['Area1'].min(), df['Area2'].min())
    df['Area1'] = ((df['Area1'] - min_area) / (df['Area1'].max() - min_area) * 100)
    df['Area2'] = ((df['Area2'] - min_area) / (df['Area2'].max() - min_area) * 100)
    plt.plot(df['File'], df['Area1'], label=disc1, marker='o', linestyle='--', color='#008080')
    plt.plot(df['File'], df['Area2'], label=disc2, marker='o', linestyle='--', color='red')
    for i, txt in enumerate(df['Area1']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['Area1'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='#008080', fontsize=8)
    for i, txt in enumerate(df['Area2']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['Area2'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='red', fontsize=8)
    plt.xlabel('Fraction')
    plt.ylabel('Area [Normalized]')
    plt.title('Area1 and Area2 vs Fraction Normalized')
    plt.legend()
    plt.xticks(rotation=90)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().spines['left'].set_color('black')
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_linestyle('-')
    plt.gca().spines['bottom'].set_linestyle('-')
    plt.gca().spines['left'].set_visible(True)
    plt.gca().spines['bottom'].set_visible(True)
    plt.gca().spines['left'].set_position(('outward', 2))
    plt.gca().spines['bottom'].set_position(('outward', 2))
    plt.gca().spines['left'].set_bounds(0, 100)
    plt.gca().spines['bottom'].set_bounds(0, len(df['File']) - 1)
    plt.tight_layout()
    plt.savefig(os.path.join(topdir, "Area vs File [Norm].pdf"))
    plt.show()

def plot_relareas_percent(df):
    plt.figure(figsize=(12, 8))
    df = df.sort_values(by='File')
    df['File'] = df['File'].apply(lambda x: x.split("_")[-2] if len(x.split("_")) > 1 else x)
    plt.plot(df['File'], df['RelArea1'], label=f'%{disc1}', marker='o', linestyle='--', color='#008080')
    plt.plot(df['File'], df['RelArea2'], label=f'%{disc2}', marker='o', linestyle='--', color='red')
    for i, txt in enumerate(df['RelArea1']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['RelArea1'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='#008080', fontsize=8)
    for i, txt in enumerate(df['RelArea2']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['RelArea2'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='red', fontsize=8)
    plt.xlabel('Fraction')
    plt.ylabel('Relative % [Normalized]')
    plt.title('Relative Areas vs Fraction')
    plt.legend()
    plt.xticks(rotation=90)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    # plt.gca().spines['left'].set_bounds(0, 100)
    # plt.gca().spines['bottom'].set_bounds(0, len(df['File']) - 1)
    plt.xlim(-1, len(df['File']) - 1)
    plt.ylim(0, 100)
    #plt.tight_layout()
    # plt.savefig(os.path.join(topdir, "Relative Areas vs File.pdf"))
    # plt.savefig(os.path.join(topdir, "Relative Areas vs File.png"))
    plt.show()

def main():
    process_files(files)
    areas_df = read_npz_files(topdir)
    output_path = os.path.join(topdir, "areas.xlsx")
    save_to_excel(areas_df, output_path)
    plot_ridgeline(areas_df, topdir)
    #plot_areas(areas_df)
    #plot_areas_percent(areas_df)
    #plot_relareas_percent(areas_df)

if __name__ == "__main__":
    main()