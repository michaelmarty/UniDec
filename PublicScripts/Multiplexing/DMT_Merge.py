import os
import re
from unidec.UniDecImporter.ImporterFactory import ImporterFactory as ImpFac
import numpy as np
import matplotlib.pyplot as plt

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
    rawfiles = []
    for file in os.listdir(folder):
        if re.match(r'.*\.raw$', file, re.IGNORECASE):
            rawfiles.append(file)
    rawfiles.sort()

    # Split each raw file by the last "_" and sort by the numeric suffix if present
    splitfiles = [sort_key(f) for f in rawfiles]
    # Get unique base names
    unique_bases = {}
    for base, num in splitfiles:
        if num == -1:
            continue
        if base not in unique_bases:
            unique_bases[base] = []
        unique_bases[base].append((base, num))
    # Flatten the list while maintaining order
    return unique_bases

def merge_looper(file_dict):
    outfiles = []
    for base in file_dict:
        nums = [num for _, num in file_dict[base] if num != -1]
        files = [f"{base}_{num}.raw" for num in sorted(nums)]
        dmtfiles = [f"{base}_{num}.dmt" for num in sorted(nums)]
        outfile = f"{base}_merged"
        merge_files(files, dmtfiles, outfile)
        outfiles.append(outfile)
        print("Merged files into:", outfile)
    return outfiles

def merge_files(rawlist, dmtlist, outfile):
    start_scan = 0
    start_time = 0
    outdata = []
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

    # Save to Numpy compressed file
    np.savez_compressed(outfile, data=np.array(outdata))


if __name__ == "__main__":

    folder = r"C:\Users\bht442\Desktop\New folder\18 Injections"
    os.chdir(folder)

    #strip_dmt_files(folder)

    filedict = get_file_list(folder)
    print(merge_looper(filedict))



