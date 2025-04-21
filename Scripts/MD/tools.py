import numpy as np
import MDAnalysis as md
import MDAnalysis.analysis.distances
import MDAnalysis.analysis.hbonds
import MDAnalysis.analysis.rms
import h5py
import sys

sys.path.append("/mnt/c/Python/unidec")
import unidec.tools as ud


def ncutoff(array, cutoff):
    boo1 = array < cutoff
    return np.sum(boo1)


def export_pos(u, hdf5file, startindex=-1, endindex=10000):
    amtbsel = "segid PROA or segid PROB or segid PROC"
    pgsel = "segid MEMB"
    phossel = "name P"
    phos = u.select_atoms(phossel)
    amtb = u.select_atoms(amtbsel)
    pg = u.select_atoms(pgsel)
    resids = pg.residues.resids
    protpos = []
    pos = []
    ppos = []
    for index, ts in enumerate(u.trajectory):
        if startindex <= index < endindex:
            repos = []
            for p in pg.residues:
                repos.append(np.mean(p.atoms.positions, axis=0))
            pos.append(repos)
            ppos.append(phos.atoms.positions)
            protpos.append(amtb.atoms.positions)

    protpos = np.array(protpos)
    pos = np.array(pos)
    ppos = np.array(ppos)

    hdf = h5py.File(hdf5file)
    group = hdf.require_group("positions")
    replace_dataset(group, "avg_lipid", pos)
    replace_dataset(group, "avg_phos", ppos)
    replace_dataset(group, "prot", protpos)
    replace_dataset(group, "resids", resids)
    hdf.close()


def extract_mins(u, hdf5file, protsel, lipsel, group_label, startindex=-1, endindex=10000):
    prot = u.select_atoms(protsel)
    lip = u.select_atoms(lipsel)
    resids = lip.residues.resids
    resnames = lip.residues.resnames
    d4 = []
    d = np.empty((len(prot), len(lip)))
    print("Trajectory Length:", len(u.trajectory), group_label)
    for index, ts in enumerate(u.trajectory):
        if startindex <= index < endindex:
            md.analysis.distances.distance_array(prot.positions, lip.positions, result=d)
            if ud.isempty(d):
                d4.append([])
            else:
                d2 = np.amin(d, axis=0)
                d3 = np.zeros(len(resids)) + 100
                for i, a in enumerate(lip.atoms):
                    resid = a.residue.resid
                    pos = ud.nearest(resids, resid)
                    if d2[i] < d3[pos]:
                        d3[pos] = d2[i]
                print(index, np.amin(d3))
                d4.append(d3)

    d4 = np.array(d4)

    hdf = h5py.File(hdf5file)
    group = hdf.require_group(group_label)
    replace_dataset(group, "min_dists", d4)
    replace_dataset(group, "resids", resids)
    replace_dataset_strings(group, "resnames", resnames)
    hdf.close()


def get_rmsf(u, hdf5file, sel, group_label, startindex=0, endindex=10000):
    lip = u.select_atoms(sel)
    names = np.unique(lip.names)
    resnames = lip.residues.resnames
    resids = lip.residues.resids

    rmsf = MDAnalysis.analysis.rms.RMSF(lip, start=startindex, stop=endindex)
    rmsf.run()

    d3 = np.zeros((len(resids), len(names))) + 100
    for i, a in enumerate(lip.atoms):
        resid = a.residue.resid
        name = a.name
        pos = ud.nearest(resids, resid)
        pos2 = int(np.where([x == name for x in names])[0][0])
        d3[pos, pos2] = rmsf.rmsf[i]

    hdf = h5py.File(hdf5file)
    group = hdf.require_group(group_label)
    replace_dataset(group, "rmsf", d3)
    replace_dataset(group, "resids", resids)
    replace_dataset_strings(group, "names", names)
    replace_dataset_strings(group, "resnames", resnames)
    hdf.close()

    pass


def extract_hbonds(u, hdf5file, protsel, lipsel, group_label, startindex=-1, endindex=10000):
    h = HBA_other(u, lipsel, protsel, start=startindex, end=endindex)
    h.run()
    h.generate_table()
    print(h.table)
    print(h.count_by_time())

    pass


def find_cutoff(mids, hist, window=5, threshold=0.1):
    data = np.transpose([mids, hist])
    hmax = np.amax(hist)
    peaks = ud.peakdetect(data, window=window, threshold=threshold)

    cutoffs = []
    npeaks = len(peaks)
    # print("Peaks:", peaks)
    for i, p in enumerate(peaks):
        if i + 1 < npeaks:
            isodat = ud.datachop(data, peaks[i, 0], peaks[i + 1, 0])
            minarg = np.argmin(isodat[:, 1])
            mpos = isodat[minarg, 0]
            cutoffs.append(mpos)
    # print("Mins:", cutoffs)

    return cutoffs, peaks[:, 0]


def cutoff_count(md, cutoffs, barray):
    means = []
    percents = []
    for c in cutoffs:
        # print(c)
        lens = []
        pers = []
        for i, d in enumerate(md):
            b1 = d < c
            try:
                b2 = barray[b1]
                tmean = np.average(b2)
                pers.append(tmean)
            except:
                pers.append(-1)
            l = np.sum(b1)
            lens.append(l)
        # lens = np.array(lens)
        m = np.mean(lens)
        means.append(m)
        try:
            t = np.mean(pers)
        except:
            t = -1
        percents.append(t)
    means = np.array(means)
    percents = np.array(percents)
    return means, percents

    """
    means = []
    for c in cutoffs:
        # print(c)
        lens = []
        for i, d in enumerate(md):
            l = np.sum(d < c)
            lens.append(l)
        lens = np.array(lens)
        m = np.mean(lens)
        means.append(m)
    means = np.array(means)"""


class HBA_other(MDAnalysis.analysis.hbonds.HydrogenBondAnalysis):
    DEFAULT_DONORS = {'GLYCAM06': ('OW', 'N3', 'NT', 'OH', 'N'),
                      'other': (
                          'OH', 'OG', 'NE2', 'OG1', 'NE', 'N', 'ND1', 'NZ', 'NH1', 'NH2', 'OW', 'ND2', 'SG', 'OH2',
                          'NE1'),
                      'CHARMM27': (
                          'OH', 'OG', 'NE2', 'OG1', 'NE', 'N', 'ND1', 'NZ', 'NH1', 'NH2', 'OW', 'ND2', 'SG', 'OH2',
                          'NE1')}
    DEFAULT_ACCEPTORS = {'GLYCAM06': ('SM', 'OH', 'NT', 'O', 'N', 'OY', 'OW', 'OS', 'O2'),
                         'other': (
                             'OH', 'OG', 'OD1', 'OD2', 'SD', 'O', 'OW', 'NE2', 'OE2', 'OC1', 'OG1', 'OC2', 'SG', 'OE1',
                             'OH2', 'ND1'),
                         'CHARMM27': (
                             'OH', 'OG', 'OD1', 'OD2', 'SD', 'O', 'OW', 'NE2', 'OE2', 'OC1', 'OG1', 'OC2', 'SG', 'OE1',
                             'OH2', 'ND1')}
