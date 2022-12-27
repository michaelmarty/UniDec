__author__ = 'Michael.Marty'

from unidec.modules.unidectools import get_importer
from unidec.modules.mzMLimporter import *
import os
import fnmatch
import h5py


def parse(path, times, timestep, volts, outputheader, directory, output="txt"):
    if os.path.isfile(path) or os.path.isdir(path):
        if output == "hdf5":
            outfile = outputheader + ".hdf5"
            outpath = os.path.join(directory, outfile)
            hdf = h5py.File(outpath, "a")
            try:
                del hdf["ms_dataset"]
            except:
                pass
            msdataset = hdf.require_group("ms_dataset")

            if volts is not None:
                msdataset.attrs["v1name"] = "Collision Voltage"
            else:
                msdataset.attrs["v1name"] = "timestart"
            msdataset.attrs["timestep"] = timestep
            config = hdf.require_group("config")
            config.attrs["metamode"] = -1

        num = 0
        d = get_importer(path)
        for v, time in enumerate(times):
            data = d.get_data(time_range=(time, time + timestep))
            if not ud.isempty(data):
                if output == "txt":
                    if volts is not None:
                        outfile = outputheader + "_" + str(int(volts[v])) + ".txt"
                    else:
                        outfile = outputheader + "_" + str(v) + ".txt"
                    if not os.path.isdir(directory):
                        os.mkdir(directory)
                    outpath = os.path.join(directory, outfile)
                    np.savetxt(outpath, data)
                    print("Saved:", outpath)
                elif output == "hdf5":
                    group = msdataset.require_group(str(v))
                    replace_dataset(group, "raw_data", data=data)
                    # group=msdataset.require_group(str(v))
                    if volts is not None:
                        group.attrs["Collision Voltage"] = volts[v]
                    group.attrs["timestart"] = time
                    group.attrs["timeend"] = time + timestep
                    group.attrs["timemid"] = time + timestep * 0.5
                    num += 1
                    pass
        if output == "hdf5":
            msdataset.attrs["num"] = num
            hdf.close()
    else:
        print("File not found:", path)
    return outpath


def parse_multiple(paths, timestep, newdir, starttp, endtp, voltsarr=None, outputname=None):
    outfile = outputname + ".hdf5"
    outpath = os.path.join(newdir, outfile)
    hdf = h5py.File(outpath, "a")
    try:
        del hdf["ms_dataset"]
    except:
        pass
    msdataset = hdf.require_group("ms_dataset")

    if voltsarr is not None:
        msdataset.attrs["v1name"] = "Collision Voltage"
    else:
        msdataset.attrs["v1name"] = "timestart"
    msdataset.attrs["timestep"] = timestep
    msdataset.attrs["v2name"] = "Original File"
    config = hdf.require_group("config")
    config.attrs["metamode"] = -1

    num = 0
    v = 0
    print(starttp, endtp, timestep)
    for path in paths:
        if os.path.isfile(path) or os.path.isdir(path):
            d = get_importer(path)
            for t in np.arange(starttp, endtp, timestep):
                data = d.get_data(time_range=(t, t + timestep))
                if not ud.isempty(data):
                    group = msdataset.require_group(str(num))
                    replace_dataset(group, "raw_data", data=data)
                    # group=msdataset.require_group(str(v))
                    try:
                        if voltsarr[v] is not None:
                            temp = int(t / timestep)
                            group.attrs["Collision Voltage"] = voltsarr[v][temp]
                    except:
                        pass
                    group.attrs["timestart"] = t
                    group.attrs["timeend"] = t + timestep
                    group.attrs["timemid"] = t + timestep * 0.5
                    splits = path.split(sep="\\")
                    group.attrs["Original File"] = splits[len(splits) - 1]
                    num += 1
                    pass
            v += 1
        else:
            print("File not found: ", path)
    msdataset.attrs["num"] = num
    hdf.close()
    return outpath


def extract(file, directory, timestep=1.0, output="txt"):
    print(file)
    path = os.path.join(directory, file)
    name = os.path.splitext(file)[0]
    newdir = os.path.join(directory, name)
    if output == "hdf5":
        newdir = directory
    # name=ud.smartdecode(name)
    splits = name.split(sep="_")
    try:
        for i, s in enumerate(splits):
            if s.lower() == "ramp":
                start = float(splits[i + 1])
                end = float(splits[i + 2])
                step = float(splits[i + 3])
                volts = np.arange(start, end + step, step)
                times = np.arange(0.0, len(volts) * timestep, timestep)
                print("Ramp:", times, volts)
                outpath = parse(path, times, timestep, volts, name, newdir, output=output)
                return outpath
    except:
        print("Error parsing ramp keyword. Ignoring.")

    d = get_importer(path)
    maxtime = d.get_max_time()

    times = np.arange(0, maxtime, timestep)
    outpath = parse(path, times, timestep, None, name, newdir, output=output)
    return outpath


def extract_scans(file, directory, scanbins=1, output="txt"):
    print(file)
    scanbins = int(float(scanbins))
    path = os.path.join(directory, file)

    if os.path.isfile(path) or os.path.isdir(path):
        d = get_importer(path)

        name = os.path.splitext(file)[0]
        newdir = os.path.join(directory, name)
        if output == "hdf5":
            newdir = directory
        if os.path.splitext(file)[1] == ".mzML":
            maxtime = d.get_max_time()
            maxscans = d.get_max_scans()
        else:
            maxtime = d.get_max_time()
            maxscans = d.get_max_scans() - 1
        scans = np.arange(0, maxscans, scanbins)

        if output == "hdf5":
            outfile = name + ".hdf5"
            outpath = os.path.join(newdir, outfile)
            hdf = h5py.File(outpath, "a")
            try:
                del hdf["ms_dataset"]
            except:
                pass
            msdataset = hdf.require_group("ms_dataset")
            msdataset.attrs["v1name"] = "timemid"
            config = hdf.require_group("config")
            config.attrs["metamode"] = -1
        num = 0

        for v, scan in enumerate(scans):
            data = d.get_data(scan_range=(scan, scan + scanbins))

            if not ud.isempty(data):
                if output == "txt":
                    outfile = name + "_" + str(v) + ".txt"
                    if not os.path.isdir(newdir):
                        os.mkdir(newdir)
                    outpath = os.path.join(newdir, outfile)
                    np.savetxt(outpath, data)
                    print("Saved:", outpath)
                elif output == "hdf5":
                    group = msdataset.require_group(str(v))
                    replace_dataset(group, "raw_data", data=data)
                    times = d.get_times_from_scans([scan, scan + scanbins])
                    group.attrs["timestart"] = times[0]
                    group.attrs["timeend"] = times[2]
                    group.attrs["timemid"] = times[1]
                    group.attrs["scanstart"] = scan
                    group.attrs["scanend"] = scan + scanbins
                    num += 1
                    pass
        if output == "hdf5":
            msdataset.attrs["num"] = num
            hdf.close()
            print(outpath)
            return outpath
    else:
        print("File not found:", path)
        return None


def extract_timepoints(files, directories, starttp=None, endtp=None, timestep=1.0, outputname="Combined"):
    paths = []
    names = []
    newdir = 0
    for f, d in zip(files, directories):
        path = os.path.join(d, f)
        name = os.path.splitext(f)[0]
        newdir = d
        paths.append(path)
        names.append(name)
    voltsarr = []
    for n in names:
        splits = n.split(sep="_")
        found = 0
        for i, s in enumerate(splits):
            if s.lower() == "ramp":
                found = 1
                try:
                    start = float(splits[i + 1])
                    end = float(splits[i + 2])
                    step = float(splits[i + 3])
                    volts = np.arange(start, end + step, step)
                    voltsarr.append(volts)
                except:
                    found = 0
                break
        if not found:
            print("File didn't have ramp")
    return parse_multiple(paths, timestep, newdir, starttp, endtp, voltsarr, outputname)


def extract_scans_multiple_files(files, dirs, startscan=1.0, endscan=1.0, outputname="Combined", existing_path=None,
                                 vars=None, keys=None):
    paths = []
    names = []
    startscan = int(float(startscan))
    endscan = int(float(endscan))
    newdir = 0
    for f, d in zip(files, dirs):
        path = os.path.join(d, f)
        name = os.path.splitext(f)[0]
        newdir = d
        paths.append(path)
        names.append(name)
    outfile = outputname + ".hdf5"
    if existing_path is None:
        outpath = os.path.join(newdir, outfile)
    else:
        outpath = existing_path
    print("Creating HDF5 File:", outpath)
    hdf = h5py.File(outpath, "a")
    try:
        if existing_path is None:
            del hdf["ms_dataset"]
    except:
        pass
    msdataset = hdf.require_group("ms_dataset")
    if vars is not None and keys is not None:
        msdataset.attrs["v1name"] = "Position"
    else:
        msdataset.attrs["v1name"] = "timemid"
    msdataset.attrs["v2name"] = "Original File"
    config = hdf.require_group("config")
    config.attrs["metamode"] = -1
    num = 0
    for path in paths:
        print("Opening", path)
        if os.path.isfile(path):
            importer = get_importer(path)
            if endscan == -1:
                endscan = importer.get_max_scans()
            data = importer.get_data(scan_range=(startscan, endscan))
            if not ud.isempty(data):
                group = msdataset.require_group(str(num))
                replace_dataset(group, "raw_data", data=data)
                times = importer.get_times_from_scans([startscan, endscan])
                group.attrs["timestart"] = times[0]
                group.attrs["timeend"] = times[2]
                group.attrs["timemid"] = times[1]
                group.attrs["scanstart"] = startscan
                group.attrs["scanend"] = endscan
                splits = path.split(sep="\\")
                group.attrs["Original File"] = splits[len(splits) - 1]
                if vars is not None and keys is not None:
                    v = vars[num]
                    for k in keys:
                        group.attrs[str(k)] = str(v[keys[k]])
                num += 1
        else:
            print("File not found: ", path)
    msdataset.attrs["num"] = num
    hdf.close()
    return outpath


def get_files(directory, timestep=1.0, output="txt"):
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, '*.mzml'):
            extract(file, directory, timestep=timestep, output=output)


if __name__ == '__main__':
    directory = "z:\\mtmarty\\Test"
    directory = "C:\\Data\\"
    directory = "Z:\\wresager\\Bleaching\\"
    file = "20170721_WCR_RHO_bleaching.RAW"
    directory = "C:\\Python\\UniDec3\\TestSpectra\\"
    file = "test_ms.raw"
    timestep = 2.0
    # extract("test.mzML", directory, timestep, output="hdf5")
    extract_scans(file, directory, scanbins=100, output="hdf5")
    # get_files(directory, output="hdf5")
    exit()
