import matplotlib.pyplot as plt
import os
import scipy.stats as stats
from unidec.LipiDec.AddFeaturestoDF import *


def calc_fold_changes(df, g1=1, g0=0, norm="", n0=3, n1=3):
    a1 = "AvgA" + norm + "_" + str(g1)
    a0 = "AvgA" + norm + "_" + str(g0)
    s1 = "SD" + norm + "_" + str(g1)
    s0 = "SD" + norm + "_" + str(g1)
    df["FoldChange"] = df[a1] / df[a0]
    df = df.sort_values("FoldChange")
    sd = df["FoldChange"] * np.sqrt((df[s1] / df[a1]) ** 2 + (df[s0] / df[a0]) ** 2)
    df["SDFold"] = sd

    m1 = df[a1].to_numpy()
    m0 = df[a0].to_numpy()
    sd1 = df[s1].to_numpy()
    sd0 = df[s0].to_numpy()
    pvals = [stats.ttest_ind_from_stats(m0[i], sd0[i], n0, m1[i], sd1[i], n1, equal_var=True)[1] for i in
             range(len(m0))]
    df["Pval"] = np.array(pvals).tolist()
    df["Log2Fold"] = pd.Series(np.log2(df["FoldChange"]))
    df["Log10P"] = pd.Series(-np.log10(df["Pval"]))
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

            mtic_norm = a / np.amax(a)
            mtic_norm_areas.append(mtic_norm)

            stdnorm = norm_to_std(df, col_name)
            stdnorms.append(stdnorm)

            mtic_norm2 = stdnorm / np.amax(stdnorm)
            mticstd.append(mtic_norm2)

        set_avg_std(df, areas, i, name="")
        set_avg_std(df, mtic_norm_areas, i, name="N1")
        set_avg_std(df, stdnorms, i, name="N2")
        set_avg_std(df, mticstd, i, name="N3")

    df = calc_fold_changes(df, norm="N3", n0=nreps, n1=nreps, g0=0, g1=1)


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
    if "IS" in classes:
        return "IS"
    else:
        return None

def find_std_class(df):
    classes = np.unique(df["Molecule List Name"])
    if "EquiSPLASH" in classes:
        return "EquiSPLASH"
    if "IS" in classes:
        return "IS"
    else:
        return None


colorscheme = {"CL": "slateblue", "PE": "goldenrod", "PG": "mediumseagreen", "LPE": "violet", "LPG" : "tomato", "IS": "grey", "EquiSPLASH": "grey"}

if __name__ == "__main__":
    file = "Z:\\Group Share\\Melanie Odenkirk\\20220815_Ecoli_d7_lowerCone\\Exported_TLs\\220825_Ecoli_dre_on-off_comparison.csv"

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
