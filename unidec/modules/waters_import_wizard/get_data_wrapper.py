import os, sys

try:
    from unidec.UniDecImporter.Waters.Waters import WatersDataImporter as WDI
except Exception as e:
    print("Error importing Waters Importer, get_data_wrapper.py", e)
import unidec.tools as ud
import numpy as np


def get_imms_data(start, end, im_bin_size, raw_file, pusher, function_no, scan_start, scan_end, directory=None):
    print(directory)
    if directory is None:
        cdcreader_path = os.getcwd()
    else:
        cdcreader_path = directory

    cdcexe = os.path.join(cdcreader_path, 'CDCReader.exe')
    print(cdcreader_path, cdcexe)
    if not os.path.isfile(cdcexe):
        cdcexe = os.path.join(cdcreader_path, 'bin', 'CDCReader.exe')
        if not os.path.isfile(cdcexe):
            print(u"Unable to find CDCReader.exe")
            print(cdcexe)
            sys.exit()

    if os.path.isfile(cdcexe):
        immsfile = os.path.splitext(raw_file)[0] + '_' + str(function_no) + '_imms.txt'
        msfile = os.path.splitext(raw_file)[0] + '_' + str(function_no) + '_ms.txt'
        print("IMMS Processing: ", raw_file)
        print("Outputs:\n\t", immsfile, "\n\t", msfile)

        call_params = [cdcexe,
                       "-r", raw_file,
                       "-m", msfile,
                       "-i", immsfile,
                       "--im_bin=" + str(im_bin_size),
                       "--fn=" + function_no,
                       "--ms_smooth_window=0",
                       "--ms_number_smooth=0",
                       "--ms_bin=" + str(im_bin_size)]

        if scan_start > 0.:
            call_params.append("--scan_start=" + str(scan_start))

        if scan_end is not None and scan_end > 0.:
            call_params.append("--scan_end=" + str(scan_end))
        try:
            call_params.append("--mass_start=" + str(int(start)))
            call_params.append("--mass_end=" + str(int(end)))
        except ValueError:
            pass

        p = ud.exe_call(call_params)
        # p = subprocess.call(call_params, shell=False, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)  # , stdout=devnull, shell=False)

        if p != 0:
            print("CONVERSION ERROR! Call Parameters:", call_params)
            print("Std out", p)

        # arrival_times = numpy.array(range(200)) * (float(pusher) / 1000.)
        return None, None, None, None, None


def new_get_data_MS(start, end, bin_size, raw_file, function_no, scan_start=None, scan_end=None, directory=None,
                    time_start=None, time_end=None):
    msfile = os.path.splitext(raw_file)[0] + '_' + str(function_no) + '_ms.txt'
    print("MS Processing: ", raw_file, "to", msfile)

    reader = WDI(raw_file, function=function_no)

    if time_start is not None:
        if time_start >= 0. and time_end is not None and time_end > 0.:
            time_range = [time_start, time_end]
        else:
            time_range = None
        print("Time Range:", time_range)
        data = reader.get_avg_scan(time_range=time_range, mzbins=bin_size)

    else:
        if scan_start > 0. and scan_end is not None and scan_end > 0.:
            scan_range = [scan_start, scan_end]
        else:
            scan_range = None

        data = reader.get_avg_scan(scan_range=scan_range, mzbins=bin_size)

    try:
        start = float(start)
        end = float(end)
    except:
        start = None
        end = None

    if start is not None and end is not None:
        data = ud.datachop(data, start, end)

    np.savetxt(msfile, data)
