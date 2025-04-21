import os
import glob
import unidec.tools as ud
from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# RC Params
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

topdir = "Z:\\Group Share\\Dash\\Synapt Runs 20241218"
files = glob.glob(f'{topdir}/*.raw', recursive=False)

mass1 = 678.3890
mass2 = 760.4525

ppm = 100
mtol1 = ppm * mass1 / 1e6
mtol2 = ppm * mass2 / 1e6

mtol = 0.02
outfile = "eictest.npz"

# DataFrame to store area1, area2, relarea1 and file name
areas_df = pd.DataFrame(columns=['File', 'Area1', 'Area2', 'RelArea1', 'RelArea2'])

def process_file(path):
    for path in files[:]:
        print(path)
        os.chdir(path)
        if True: #not os.path.isfile(outfile):
            try:
                d = ImporterFactory.create_importer(path)
            except Exception as e:
                print("File Failed", path)
                print(e)
                return None

            tic1 = d.get_eic(mass1, mtol)
            tic2 = d.get_eic(mass2, mtol)

            tic1_peak_time = tic1[np.argmax(tic1[:, 1]), 0]
            tic2_peak_time = tic2[np.argmax(tic2[:, 1]), 0]

            expected_time1 = 1.7
            expected_time2 = 1.6

            if abs(tic1_peak_time - expected_time1) > 0.2 or abs(tic2_peak_time - expected_time2) > 0.2:
                print("Peak Time1:", tic1_peak_time, "Peak Time2:", tic2_peak_time)
                return None

            sd1 = 0.2
            t1_range = [tic1_peak_time - sd1, tic1_peak_time + sd1]
            sd2 = 0.2
            t2_range = [tic2_peak_time - sd2, tic2_peak_time + sd2]

            area1 = ud.integrate(tic1, t1_range[0], t1_range[1])[0]
            area2 = ud.integrate(tic2, t2_range[0], t2_range[1])[0]
            relarea1 = (area1 / (area1 + area2)) * 100
            relarea2 = (area2 / (area1 + area2)) * 100

            np.savez_compressed(os.path.join(path, outfile), array1=tic1, array2=tic2)

            return [os.path.basename(path), area1, area2, relarea1, relarea2]

def process_files(files):
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_file, files))
    for result in results:
        if result:
            areas_df.loc[len(areas_df)] = result

def save_to_excel(df, output_path):
    df.to_excel(output_path, index=False)
    print(f"Data saved to {output_path}")

def plot_areas(df):
    plt.figure(figsize=(12, 8))
    df['File'] = df['File'].apply(lambda x: x.split("_")[-2] if len(x.split("_")) > 1 else x)
    plt.plot(df['File'], df['Area1'], label='DMPC(-)', marker='o', linestyle='--', color='blue')
    plt.plot(df['File'], df['Area2'], label='POPC(+)', marker='o', linestyle='--', color='red')
    for i, txt in enumerate(df['Area1']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['Area1'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='blue', fontsize=8)
    for i, txt in enumerate(df['Area2']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['Area2'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='red', fontsize=8)
    plt.xlabel('Fraction')
    plt.ylabel('Area [Raw]')
    plt.title('Area1 and Area2 vs Fraction')
    plt.legend()
    plt.xticks(rotation=90)
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
    plt.gca().spines['left'].set_bounds(0, df['Area1'].max())
    plt.gca().spines['bottom'].set_bounds(0, len(df['File']) - 1)
    plt.tight_layout()
    plt.savefig(os.path.join(topdir, "Area vs File[raw]_test.pdf"))
    plt.show()

def plot_areas_percent(df):
    plt.figure(figsize=(12, 8))
    df['File'] = df['File'].apply(lambda x: x.split("_")[-2] if len(x.split("_")) > 1 else x)
    min_area = min(df['Area1'].min(), df['Area2'].min())
    df['Area1'] = ((df['Area1'] - min_area) / (df['Area1'].max() - min_area) * 100)
    df['Area2'] = ((df['Area2'] - min_area) / (df['Area2'].max() - min_area) * 100)
    plt.plot(df['File'], df['Area1'], label='DMPC(-)', marker='o', linestyle='--', color='blue')
    plt.plot(df['File'], df['Area2'], label='POPC(+)', marker='o', linestyle='--', color='red')
    for i, txt in enumerate(df['Area1']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['Area1'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='blue', fontsize=8)
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
    plt.savefig(os.path.join(topdir, "Area vs File [Norm]_test.pdf"))
    plt.show()

def plot_relareas_percent(df):
    plt.figure(figsize=(12, 8))
    df['File'] = df['File'].apply(lambda x: x.split("_")[-2] if len(x.split("_")) > 1 else x)
    plt.plot(df['File'], df['RelArea1'], label='%DMPC', marker='o', linestyle='--', color='blue')
    plt.plot(df['File'], df['RelArea2'], label='%POPC', marker='o', linestyle='--', color='red')
    for i, txt in enumerate(df['RelArea1']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['RelArea1'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='blue', fontsize=8)
    for i, txt in enumerate(df['RelArea2']):
        plt.annotate(f'{txt:.2f}', (df['File'][i], df['RelArea2'][i]), textcoords="offset points", xytext=(0,10), ha='center', color='red', fontsize=8)
    plt.xlabel('Fraction')
    plt.ylabel('Relative % [Normalized]')
    plt.title('Relative Areas vs Fraction')
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
    plt.savefig(os.path.join(topdir, "Relative Areas vs File_test.pdf"))
    plt.show()

def main():
    process_files(files)
    output_path = os.path.join(topdir, "areastest.xlsx")
    save_to_excel(areas_df, output_path)
    plot_areas(areas_df)
    plot_areas_percent(areas_df)
    plot_relareas_percent(areas_df)

if __name__ == "__main__":
    main()