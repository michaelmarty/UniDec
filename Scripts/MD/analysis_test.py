from tools import *

topdir = "/mnt/c/VMD/AmtB_Nanodiscs/"
dirs = ["POPG_E3_4NH2_all_atom", "POPC_E3_4NH2_all_atom", "PGPC_E3_4NH2_all_atom"]

#dirs = ["PGPC_E3_4NH2_all_atom"]
for i, d in enumerate(dirs):
    directory = os.path.join(topdir, d)

    os.chdir(directory)

    psffile = "step5_fixed.psf"
    dcdfile = "namd/aligned.dcd"
    hdf5file = "analysis.hdf5"

    u = md.Universe(psffile, dcdfile)

    pcatsel = "(segid PROA PROB PROC) and (type NH2 NH3)"  # and (resname ARG LYS) and (name NH* NZ)"
    pansel = "(segid PROA PROB PROC) and (type OC)"  # and (resname ASP GLU) and (name OE* OD*)"
    amtbsel = "segid PROA or segid PROB or segid PROC"
    lipsel = "segid MEMB"
    lcatsel = "segid MEMB and type NTL"
    lansel = "segid MEMB and type PL"
    lhaccsel = "(segid MEMB) and (type OBL OSL O2L OSLP OHL)"
    lhdonsel = "(segid MEMB) and (type OHL)"
    phdonsel = "(segid PROA PROB PROC) and (name OH OG OG1 NE2 NE N ND1 NZ NH1 NH2 ND2 NE1 SG)"
    phaccsel = "(segid PROA PROB PROC) and (name OH OG OG1 OD1 OD2 OE1 OE2 O NE2 ND1 SG SD)"
    lipheads = "(segid MEMB) and (name C1 C2 C3 C11 C12 C13 C14 C15 C21 C31 O* N* P*)"

    lnoh = "(segid MEMB) and not (name H*)"
    pnoh = "(segid PROA PROB PROC) and not (name H*)"

    sindex = 100
    eindex = 1000

    extract_mins(u, hdf5file, phaccsel, lhdonsel, "lipid_hdon", startindex=sindex, endindex=eindex)
    extract_mins(u, hdf5file, phdonsel, lhaccsel, "lipid_hacc", startindex=sindex, endindex=eindex)
    #break
    get_rmsf(u, hdf5file, lipheads, "rmsf", sindex, eindex)
    extract_mins(u, hdf5file, pnoh, lipheads, "lipid_heads", startindex=sindex, endindex=eindex)
    extract_mins(u, hdf5file, pnoh, lnoh, "all_noh", startindex=sindex, endindex=eindex)

    extract_mins(u, hdf5file, pansel, lcatsel, "lipid_cation", startindex=sindex, endindex=eindex)
    extract_mins(u, hdf5file, pcatsel, lansel, "lipid_anion", startindex=sindex, endindex=eindex)
    extract_mins(u, hdf5file, amtbsel, lipsel, "all_atoms", startindex=sindex, endindex=eindex)
    export_pos(u, hdf5file, startindex=sindex, endindex=eindex)

