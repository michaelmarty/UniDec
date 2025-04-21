# noinspection PyUnresolvedReferences
import pandas as pd
# noinspection PyUnresolvedReferences
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
from Scripts.LipiDec.AddFeaturestoDF import *

norms = ["", "N1", "N2", "N3"]

def calc_fold_changes_all(df, g1=1, g0=0, n0=3, n1=3):
    for norm in norms:
        a1 = "AvgA" + norm + "_" + str(g1)
        a0 = "AvgA" + norm + "_" + str(g0)
        s1 = "SD" + norm + "_" + str(g1)
        s0 = "SD" + norm + "_" + str(g1)
        df["FoldChange"+norm] = df[a1] / df[a0]
        df = df.sort_values("FoldChange"+norm)
        sd = df["FoldChange"+norm] * np.sqrt((df[s1] / df[a1]) ** 2 + (df[s0] / df[a0]) ** 2)
        df["SDFold"+norm] = sd

        m1 = df[a1].to_numpy()
        m0 = df[a0].to_numpy()
        sd1 = df[s1].to_numpy()
        sd0 = df[s0].to_numpy()
        pvals = [stats.ttest_ind_from_stats(m0[i], sd0[i], n0, m1[i], sd1[i], n1, equal_var=True)[1] for i in
                 range(len(m0))]
        df["Pval"+norm] = np.array(pvals).tolist()
        df["Log2Fold"+norm] = pd.Series(np.log2(df["FoldChange"]))
        df["Log10P"+norm] = pd.Series(-np.log10(df["Pval"]))
    return df


def find_file_headers(df):
    keys = df.keys()
    filenames = []
    for k in keys:
        if "Total Ion Current Area" in k:
            f = k[:-23]
            filenames.append(f)
    return np.array(filenames)


def set_colors(df):
    colors = []
    for i, row in df.iterrows():
        c = row["Molecule List Name"]
        if c in colorscheme:
            color = colorscheme[c]
        else:
            color = "grey"
        colors.append(color)
    df["colors"] = colors
    return df


def set_stds(df):
    df = df.copy(deep=True)
    stds = []
    stdclassname=find_std_class(df)
    print("Standard Class:", stdclassname)
    stddf = df[df["Molecule List Name"] == stdclassname]
    stdnames = stddf["Molecule Name"]
    for i, row in df.iterrows():
        n = row["Molecule List Name"]
        l = len(n)
        std = ''
        for j, stdname in enumerate(stdnames):
            if stdname[:l + 1] == n + " ":
                std = stdname
        if row["Molecule List Name"] == stdclassname:
            std = row["Molecule Name"]
        stds.append(std)
    df["Std"] = stds
    return df


def norm_to_std(df, col_name):
    stdnorm = []
    stddf = df[df["Molecule List Name"] == find_std_class(df)]
    for i, row in df.iterrows():
        a = float(row[col_name])
        stdname = row["Std"]
        stdtest = stddf["Molecule Name"].str.contains(stdname)
        stdrow = stddf[stdtest]
        try:
            stda = stdrow[col_name].to_numpy()[0]
            norma = a / stda
        except:
            print("Std Not Found:", stdname)
            norma = 0
        stdnorm.append(norma)
    return np.array(stdnorm)


def set_avg_std(df, areas, i, name=""):
    areas = np.array(areas)
    avgareas = np.mean(areas, axis=0)
    sd = np.std(areas, axis=0)
    df["AvgA" + name + "_" + str(i)] = avgareas
    df["SD" + name + "_" + str(i)] = sd

def int_filter(df, threshold, norm="N3", g1=1, g0=0):
    a1 = "AvgA" + norm + "_" + str(g1)
    a0 = "AvgA" + norm + "_" + str(g0)

    b1 = df[a1]>threshold
    b0 = df[a0]>threshold

    return df[b1 + b0]


def sd_filter(df, threshold):
    b1 = df["SDFold"] < threshold
    return df[b1]

def set_avg_rt_mz(df):
    keys = df.keys()
    rtcols = []
    mzcols = []
    for k in keys:
        if "Retention Time" in k:
            rtcols.append(k)
        if "Mass Error" in k:
            mzcols.append(k)
    df["Average Rt(min)"] = df[rtcols].mean(axis=1)
    df["Average Mass Error"] = df[mzcols].mean(axis=1)
    df["Average Mz"]= df["Average Mass Error"] / 1e6 * df["Precursor Mz"] + df["Precursor Mz"]
    return df

def translate(df):
    df["Ontology"] = df["Molecule List Name"]
    return df

def filter_skyline(df):
    b1 = df["Fragment Ion"] == "precursor"
    try:
        b2 = df["Fragment Ion Type"] == "precursor"
    except:
        print('Fragmention Ion Type not found. Skipping Filtering')
        b2 = b1
    b3 = np.logical_not(df["Molecule Name"].str.contains("Unsettled"))
    return df[b1 * b3 * b2]

def rel_quant(df, ngroups, nreps):
    headers = find_file_headers(df)
    grouped = headers.reshape((ngroups, nreps))
    std_class = find_std_class(df)
    not_std_bool = df["Molecule List Name"] != std_class
    std_bool = df["Molecule List Name"] == std_class
    for i, g in enumerate(grouped):
        areas = []
        mtic_norm_areas = []
        stdnorms = []
        mticstd = []
        for j, f in enumerate(g):
            print(i, j, f)
            col_name = f + " Area"
            a = df[col_name]
            areas.append(a)

            loga = np.log(a)

            mtic_norm = a / np.sum(a[not_std_bool])
            #mtic_norm = np.exp(loga - np.average(loga[not_std_bool]))
            mtic_norm_areas.append(mtic_norm)

            stdnorm = norm_to_std(df, col_name)
            stdnorms.append(stdnorm)

            mtic_norm2 = stdnorm / np.sum(stdnorm[not_std_bool])
            mticstd.append(mtic_norm2)

        set_avg_std(df, areas, i, name="")
        set_avg_std(df, mtic_norm_areas, i, name="N1")
        set_avg_std(df, stdnorms, i, name="N2")
        set_avg_std(df, mticstd, i, name="N3")

    df = calc_fold_changes_all(df, n0=nreps, n1=nreps, g0=0, g1=1)
    return df

def full_processing(df):
    # Filtering
    df = filter_skyline(df)

    df = set_stds(df)
    df = set_tails(df, "Molecule Name")
    df = parse_names(df, "Molecule Name", parse_adduct=False)
    df = set_avg_rt_mz(df)
    df = mass_defect_df(df, ref_col="Precursor Mz", z_col="Precursor Charge")
    df = translate(df)
    df = set_colors(df)

    return df

def find_std_class(df):
    classes = np.unique(df["Molecule List Name"])
    if "EquiSPLASH" in classes:
        return "EquiSPLASH"
    elif "IS" in classes:
        return "IS"
    elif "LightSPLASH" in classes:
        return "LightSPLASH"
    else:
        return None


colorscheme = {"CL": "blue",
               "Cer": "dodgerblue",
               "EqualSplash": "gray",
                "HexCer": "rosybrown",
                "LPG": "tomato",
               "LPC": "darkmagenta",
               "LPE": "violet",
               "LPEO": "plum",
               "PC": "maroon",
               "PCO": "red",
               "PE": "orange",
               "PEO": "lightcoral",
               "PEP": "goldenrod",
               "PI": "lightseagreen",
               "PG": "green",
               "PIO": "darkturquoise",
               "PS": "lightskyblue",
               "PSO": "darkblue",
               "SM": "mediumvioletred",
               "TG": "lightpink",
               "None": "grey"}

if __name__ == "__main__":
    file = "Z:\\Group Share\\Melanie Odenkirk\\20220815_Ecoli_d7_lowerCone\\Exported_TLs\\221117_Brain_D1_Temperature_fixed_Areps_filtered.csv"
    #file = "Z:\\Group Share\\Melanie Odenkirk\\20220815_Ecoli_d7_lowerCone\\Exported_TLs\\220816_DIA_Ecoli_lipid_ND_res_outputs.csv"
    #file = "Z:\\Group Share\\Melanie Odenkirk\\20220815_Ecoli_d7_lowerCone\\Exported_TLs\\220816_SONAR_Ecoli_lipid_ND_res_outputs.csv"
    #file = "Z:\\Group Share\\Melanie Odenkirk\\20220815_Ecoli_d7_lowerCone\\Exported_TLs\\220816_IMS_DIA_Ecoli_lipid_ND_res_outputs.csv"
    os.chdir(os.path.dirname(file))
    outfile = file[:-4]+"_filtered.csv"

    nreps = 3
    ngroups = 2

    df = pd.read_csv(file)
    print(df.keys())

    # Processing
    df = full_processing(df)

    # Quant
    df = rel_quant(df, ngroups, nreps)

    #df = int_filter(df, 1e-3)
    #df = sd_filter(df, 1)

    labels = df["Simp Name"]
    fold = df["FoldChange"]
    sd = df["SDFold"]

    '''
    plt.scatter(np.abs(df["Average Mass Error"]), fold)
    plt.show()
    exit()'''


    plt.bar(np.arange(len(labels)), fold, yerr=sd, color=df.colors)
    plt.hlines(1, 0, len(labels), linestyle="--")
    plt.gca().set_xticks(np.arange(len(labels)))
    plt.gca().set_xticklabels(labels, rotation="vertical")
    plt.gca().set_ylabel("Fold Change")
    plt.tight_layout()
    df.to_csv(outfile)
    plt.show()

    image_format = 'eps'  # e.g .png, .svg, etc.
    image_name = 'myimage.eps'

    plt.savefig(image_name, format=image_format, dpi=1200
                )
