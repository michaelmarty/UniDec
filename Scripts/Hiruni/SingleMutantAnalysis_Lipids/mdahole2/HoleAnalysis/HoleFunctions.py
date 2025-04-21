import MDAnalysis as mda
import numpy as np
import os
from mdahole2.analysis import HoleAnalysis
import pandas as pd
import warnings

warnings.filterwarnings("ignore")
import mdahole2.analysis as ha
from scipy.interpolate import interp1d


def mergedata(data1, data2):
    """
    Uses a 1D intepolation to merge the intensity of data2 into the axis of data1
    :param data1: Reference data (N x 2 array)
    :param data2: Data to merge (M x 2 array)
    :return: data2 interpolated into the same axis as data1 (N x 2)
    """
    f = interp1d(data2[:, 0], data2[:, 1], bounds_error=False, fill_value=0)
    intx = data1[:, 0]
    inty = f(intx)
    newdat = np.column_stack((intx, inty))
    return newdat


def process_hole_inputs(file):
    print(file)
    try:
        results = ha.utils.collect_hole(file)[0]
        rxn = np.array([r[0] for r in results])
        radius = np.array([r[1] for r in results])
        line = np.array([r[2] for r in results])
        df = pd.DataFrame({"rxn": rxn, "radius": radius, "line": line})
        return df
    except Exception as e:
        print(e)
        return pd.DataFrame({"rxn": [], "radius": [], "line": []})


def process_hole_trajectory(header, xcol="rxn", ycol="radius"):
    # find files
    files = [f for f in os.listdir() if f.startswith(header) and f.endswith(".out")]
    frames = [int(f.split(".")[0].split(header)[1]) for f in files]
    dfs = [process_hole_inputs(f) for f in files]
    # Interpolate the data
    template = np.transpose([dfs[0][xcol].values, dfs[0][ycol].values])
    outdata = []
    for i, df in enumerate(dfs):
        if len(df) == 0:
            continue
        data = np.transpose([df[xcol].values, df[ycol].values])
        merged = mergedata(template, data)
        outdata.append(merged[:, 1])
    avg = np.average(outdata, axis=0)
    std = np.std(outdata, axis=0)
    rxn = dfs[0][xcol].values
    return np.transpose([rxn, avg, std]), outdata


def process_multiholes(labels=["A", "B", "C", "D"], xcol="rxn", ycol="radius"):
    outdata = []
    eachone = []
    for i, label in enumerate(labels):
        data, alldat = process_hole_trajectory(label + "hole", xcol, ycol)
        eachone.append(data)
        if i == 0:
            template = data
            for dat in alldat:
                outdata.append(dat)
        else:
            for dat in alldat:
                merged = mergedata(template, np.transpose([data[:, 0], dat]))[:, 1]
                outdata.append(merged)
    rxn = template[:, 0]
    avg = np.average(outdata, axis=0)
    std = np.std(outdata, axis=0)
    return np.transpose([rxn, avg, std]), outdata, eachone


def run_hole(file1, file2, cps=None, stopframe=None, labels=['A', 'B', 'C', 'D']):
    u = mda.Universe(file1, file2)
    print("# of frames:", len(u.trajectory))
    if cps is None:
        cps = guess_cpoint(u)

    if not os.path.exists('./hole'):
        os.mkdir('./hole')

    if stopframe is None:
        stopframe = len(u.trajectory)

    df = pd.DataFrame()
    # add column for frames
    df['frame'] = range(stopframe)
    # add columns for each label
    for label in labels:
        df[label] = np.nan

    for i, cp in enumerate(cps):
        label = labels[i]
        print(label)
        ha = HoleAnalysis(u, executable='/mnt/c/VMD/hole2/exe/hole', prefix="./hole/" + label, cpoint=cp,
                          cvect=[0, 0, -1])
        ha.run(start=0, stop=int(stopframe))
        mr = ha.min_radius()
        print(mr)
        df[label] = mr[:, 1]
        np.savetxt("./hole/minrad" + label + ".csv", mr, delimiter=",")
        try:
            ha.create_vmd_surface(filename="./hole/hole" + label + ".vmd")
        except Exception as e:
            print(e, "Error in creating vmd surface")

    df.to_excel("./hole/minrad.xlsx", index=False)


def stitch_holes(labels=['A', 'B', 'C', 'D']):
    df = pd.DataFrame()
    for i, label in enumerate(labels):
        mr = np.loadtxt("./hole/minrad" + label + ".csv", delimiter=",")
        if i == 0:
            df['frame'] = mr[:, 0]
        df[label] = mr[:, 1]
    df.to_excel("./hole/minrad.xlsx", index=False)


def guess_cpoint(u):
    segs = ["A", "B", "C", "D"]
    avgs = []
    for seg in segs:
        segname = "segid PRO" + seg
        selection1 = "resid 189 and name HH22"
        selection2 = "resid 174 and name HD2"
        selection3 = "resid 183 and name O"
        selection4 = "resid 43 and name HE2"

        s = ""
        for i, selection in enumerate([selection1, selection2, selection3, selection4]):
            s += "(" + segname + " and " + selection
            if i != 3:
                s += ") or "
            else:
                s += ")"

        sel = u.select_atoms(s)

        centers = []
        for ts in u.trajectory:
            centers.append(sel.center_of_geometry())
        centers = np.array(centers)
        avg = np.average(centers, axis=0)
        print("PRO" + seg + ":", avg)
        avgs.append(avg)
    return avgs


if __name__ == "__main__":
    dir = "/mnt/z/Group Share/Group/VMD/AqpZMutants/WT_AA_PE_PG_CL_2/namd"
    os.chdir(dir)

    file1 = "step5_input.pdb"
    file2 = "wca_s5.dcd"

    u = mda.Universe(file1, file2)

