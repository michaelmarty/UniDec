__author__ = 'Michael.Marty'

from unidec_modules.mzMLimporter import *
from unidec_modules.unidectools import removeduplicates
import os
import fnmatch
import string
import h5py
from unidec_modules.hdf5_tools import *

try:
    from unidec_modules.data_reader import *
except:
    print "Could not import data reader: mzmlparse_auto"


def parse(path, times, timestep, volts, outputheader, directory, output="txt"):
    if os.path.isfile(path):
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

        for v, time in enumerate(times):
            if os.path.splitext(path)[1] == ".mzML":
                data = mzMLimporter(path).get_data(time_range=(time, time + timestep))
            else:
                data = DataImporter(path).get_data(time_range=(time, time + timestep))
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
                    print "Saved:", outpath
                elif output == "hdf5":
                    group = msdataset.require_group(str(v))
                    replace_dataset(group, "raw_data", data=data)
                    # group=msdataset.require_group(str(v))
                    if volts is not None:
                        group.attrs["Collision Voltage"] = volts[v]
                    group.attrs["timestart"] = time
                    group.attrs["timeend"] = time + timestep
                    num += 1
                    pass
        if output == "hdf5":
            msdataset.attrs["num"] = num
            hdf.close()
    else:
        print "File not found:", path


def extract(file, directory, timestep=1.0, output="txt"):
    print file
    path = os.path.join(directory, file)
    name = os.path.splitext(file)[0]
    newdir = os.path.join(directory, name)
    if output == "hdf5":
        newdir = directory
    splits = string.split(name, sep="_")
    for i, s in enumerate(splits):
        if s.lower() == "ramp":
            start = float(splits[i + 1])
            end = float(splits[i + 2])
            step = float(splits[i + 3])
            volts = np.arange(start, end + step, step)
            times = np.arange(0.0, len(volts) * timestep, timestep)
            print "Ramp:", times, volts
            parse(path, times, timestep, volts, name, newdir, output=output)
            return True
    if os.path.splitext(file)[1] == ".mzML":
        maxtime = mzMLimporter(path).get_max_time()
    else:
        maxtime = DataImporter(path).get_max_time()

    times = np.arange(0, maxtime, timestep)
    parse(path, times, timestep, None, name, newdir, output=output)


def extract_scans(file, directory, scanbins=1, output="txt"):
    print file
    scanbins = int(float(scanbins))
    path = os.path.join(directory, file)
    name = os.path.splitext(file)[0]
    newdir = os.path.join(directory, name)
    if output == "hdf5":
        newdir = directory
    if os.path.splitext(file)[1] == ".mzML":
        maxtime = mzMLimporter(path).get_max_time()
        maxscans = mzMLimporter(path).get_max_scans()
    else:
        maxtime = DataImporter(path).get_max_time()
        maxscans = DataImporter(path).get_max_scans()
    print maxscans
    scans = np.arange(0, maxscans, scanbins)
    if os.path.isfile(path):
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
            if os.path.splitext(path)[1] == ".mzML":
                d = mzMLimporter(path)
            else:
                d = DataImporter(path)

            data = d.get_data(scan_range=(scan, scan + scanbins))

            if not ud.isempty(data):
                if output == "txt":
                    outfile = name + "_" + str(v) + ".txt"
                    if not os.path.isdir(newdir):
                        os.mkdir(newdir)
                    outpath = os.path.join(newdir, outfile)
                    np.savetxt(outpath, data)
                    print "Saved:", outpath
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
    else:
        print "File not found:", path


def get_files(directory, timestep=1.0, output="txt"):
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, '*.mzml'):
            extract(file, directory, timestep=timestep, output=output)


if __name__ == '__main__':
    directory = "z:\\mtmarty\\Test"
    directory = "C:\\Data\\"
    directory = "Z:\\wresager\\Bleaching\\"
    file = "20170721_WCR_RHO_bleaching.RAW"
    timestep = 2.0
    # extract("test.mzML", directory, timestep, output="hdf5")
    extract_scans(file, directory, scanbins=100, output="hdf5")
    # get_files(directory, output="hdf5")
    exit()
