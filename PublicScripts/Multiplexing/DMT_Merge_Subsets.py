import os
import re
from unidec.UniDecImporter.ImporterFactory import ImporterFactory as ImpFac
import numpy as np

def strip_dmt_files(folder):
    for file in os.listdir(folder):
        if re.match(r'.*\.dmt$', file, re.IGNORECASE):
            newname = os.path.splitext(file)[0]
            newname = newname.rsplit('_', 1)[0] + '.dmt'
            print("Editing DMT File:", file, "to", newname)
            os.rename(file, newname)

def sort_key(filename):
    filename = os.path.splitext(filename)[0]
    parts = filename.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        return (parts[0], int(parts[1]))
    else:
        return (filename, -1)

def get_file_list(folder):
    """Return a dict mapping base filename -> sorted list of numeric suffixes.

    Only files matching <base>_<number>.raw are included. Files without a
    numeric suffix are ignored for merging.
    """
    rawfiles = []
    for file in os.listdir(folder):
        if re.match(r'.*\.raw$', file, re.IGNORECASE):
            rawfiles.append(file)
    rawfiles.sort()

    # Split each raw file by the last "_" and collect numeric suffixes
    bases = {}
    for f in rawfiles:
        base, num = sort_key(f)
        if num == -1:
            # ignore files without numeric suffix for merging
            continue
        bases.setdefault(base, []).append(num)

    # Sort the numeric lists
    for base in list(bases.keys()):
        bases[base] = sorted(bases[base])

    return bases

def merge_looper(file_dict):
    """For each base, create incremental subset merges.

    If a base has N injections (e.g. base_1.raw ... base_N.raw) then this
    creates N merges: first 1 injection, first 2 injections, ..., first N
    injections. Each merge is written into its own output folder named
    <base>_merged_<k>_unidecfiles (and a .npz saved there).
    """
    outfiles = []
    for base, nums in file_dict.items():
        if not nums:
            continue
        # For each subset size k (1..N) merge the first k injections
        # Create a parent folder for this base's merged outputs
        base_outdir = f"{base}_merged"
        os.makedirs(base_outdir, exist_ok=True)

        for k in range(1, len(nums) + 1):
            subset = sorted(nums)[:k]
            files = [f"{base}_{num}.raw" for num in subset]
            dmtfiles = [f"{base}_{num}.dmt" for num in subset]
            outfile = f"{base}_merged_{k}"

            # Save outputs directly under the base_outdir so the .npz and the
            # corresponding "<outfile>_unidecfiles" folder live together and are
            # not placed inside separate '1 injection', '2 injection', ... folders.
            parent_save_dir = base_outdir
            print(f"Merging {len(files)} files for base '{base}' into: {parent_save_dir} (outfile: {outfile})")
            try:
                # Pass the parent directory so merge_files will create the unidecfiles
                # folder inside the base_outdir and save the .npz there as well.
                merge_files(files, dmtfiles, outfile, parent_dir=parent_save_dir)
                outfiles.append(os.path.join(parent_save_dir, outfile))
            except Exception as e:
                print(f"Error merging {outfile} into {parent_save_dir}: {e}")
    return outfiles

def merge_files(rawlist, dmtlist, outfile, parent_dir=None):
    start_scan = 0
    start_time = 0
    outdata = []
    if len(rawlist) != len(dmtlist):
        raise ValueError("Number of raw files and dmt files must match")
    for i, raw in enumerate(rawlist):
        dmt = dmtlist[i]
        rawimp = ImpFac.create_importer(raw, silent=True)
        dmtimp = ImpFac.create_importer(dmt)
        dmtdat = dmtimp.get_cdms_data()
        print("DMT Data Length:", len(dmtdat), start_scan, start_time)
        dmttimes = np.array([rawimp.get_scan_time(scan) for scan in dmtdat[:,2]])
        dmttimes += start_time
        # Add dmttimes column
        dmtdat[:,4] = dmttimes
        dmtdat[:,2] += start_scan
        outdata.append(dmtdat)
        start_scan = np.amax(rawimp.scans) + start_scan
        start_time = np.amax(rawimp.times) + start_time

    # Concatenate all data
    outdata = np.vstack(outdata)

    # Print final scan and time
    print("Final Scan:", start_scan)
    print("Final Time:", start_time)

    # Write a _conf.dat text file that has
    HTanalysistime = start_time
    HTmaxscans = start_scan

    # Create unidecfiles folder. The conf (.dat) should live inside the
    # unidecfiles folder. The .npz should be placed in the parent_dir (i.e.
    # the 'k injection' folder) if provided; otherwise it will be placed
    # inside the unidecfiles folder (legacy behavior).
    if parent_dir:
        unidec_folder = os.path.join(parent_dir, outfile + "_unidecfiles")
    else:
        unidec_folder = outfile + "_unidecfiles"

    os.makedirs(unidec_folder, exist_ok=True)

    confout = os.path.join(unidec_folder, outfile + "_conf.dat")
    with open(confout, "w") as f:
        f.write(f"HTanalysistime {HTanalysistime}\n")
        f.write(f"HTmaxscans {HTmaxscans}\n")

    # Determine where to save the merged .npz: in parent_dir if provided
    if parent_dir:
        npz_path = os.path.join(parent_dir, outfile + ".npz")
    else:
        npz_path = os.path.join(unidec_folder, outfile + ".npz")

    np.savez_compressed(npz_path, data=np.array(outdata))
    print("Saved merged data to:", npz_path)
    print("Saved conf to:", confout)


if __name__ == "__main__":
    import sys

    # Default folder (edit as needed). You can override by passing a folder
    # path as the first command-line argument.
    # Use a raw string for Windows path to avoid escape-sequence warnings
    folder = r"Z:\Group Share\BHT\Temp\Paper\DMT\Amgen IgM-like Fusion Protein 21 Injections Set 1"
    if len(sys.argv) > 1:
        folder = sys.argv[1]

    os.chdir(folder)

    #strip_dmt_files(folder)

    filedict = get_file_list(folder)
    print(merge_looper(filedict))



