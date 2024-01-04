import os, subprocess, sys

try:
    from unidec.modules.waters_importer.WatersImporter import WatersDataImporter as WDI
except Exception as e:
    print("Error importing Waters Importer, get_data_wrapper.py", e)
import unidec.tools as ud
import numpy as np


def new_get_data(start, end, im_bin_size, raw_file, pusher, function_no, scan_start, scan_end, directory=None):
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

        p = subprocess.call(call_params, shell=False, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)  # , stdout=devnull, shell=False)

        if p != 0:
            print("CONVERSION ERROR! Call Parameters:", call_params)
            print("Std out", p)

        try:
            '''
            # This allows the data to be imported if you wanted to use it
            # msfile=raw_file + '_' + str(function_no) + '_ms.txt'
            ms_data = numpy.loadtxt(msfile, dtype='float')
            A = ms_data[:, 0]
            B = ms_data[:, 1]

            # immsfile=raw_file + '_' + str(function_no) + '_imms.txt'
            imms_data = numpy.loadtxt(immsfile, dtype='float')
            mymap2 = dict((mz_value, numpy.zeros(200)) for mz_value in numpy.unique(imms_data[:, 0]))

            for row in imms_data:
                mymap2[row[0]][row[1]] = row[2]

            mz_values = mymap2.keys()
            mz_values.sort()
            if pusher is not None:
                arrival_times = numpy.array(range(200)) * (float(pusher) / 1000.)
            else:
                arrival_times = numpy.array(range(200))
                print "No Pusher Provided, Using Bin Number"

            X, Y, C = [], [], []
            for mz in mz_values:
                X.append(numpy.ones(200) * mz)
                Y.append(arrival_times)
                C.append(mymap2[mz])

            X = numpy.array(X)
            Y = numpy.array(Y)
            C = numpy.array(C)

            return A, B, X, Y, C
            '''
            return None, None, None, None, None
        except Exception as e:
            print("ERROR")
            print(e)
            return None, None, None, None, None


'''
def new_get_data_MS_old(start, end, bin_size, raw_file, function_no, scan_start, scan_end, dir=None):
    if dir is None:
        reader_path = os.getcwd()
    else:
        reader_path = dir

    rawexe = os.path.join(reader_path, 'rawreadertim.exe')

    if not os.path.isfile(rawexe):
        rawexe = os.path.join(reader_path, "bin", 'rawreadertim.exe')
        if not os.path.isfile(rawexe):
            print("Unable to find rawreadertim.exe")
            print(rawexe)
            sys.exit()

    if os.path.isfile(rawexe):
        msfile = os.path.splitext(raw_file)[0] + '_' + str(function_no) + '_ms.txt'
        print("MS Processing: ", raw_file, "to", msfile)

        call_params = [rawexe,
                       "-r", raw_file,
                       "-m", msfile,
                       "--fn=" + function_no,
                       "--ms_smooth_window=0",
                       "--ms_number_smooth=0",
                       "--ms_bin=" + str(bin_size)]

        if scan_start > 0.:
            call_params.append("--scan_start=" + str(scan_start))

        if scan_end is not None and scan_end > 0.:
            call_params.append("--scan_end=" + str(scan_end))

        try:
            call_params.append("--mass_start=" + str(int(start)))
            call_params.append("--mass_end=" + str(int(end)))
        except ValueError:
            pass

        p = subprocess.call(call_params, shell=False)

        if p != 0:
            print("CONVERSION ERROR! Call Parameters:", call_params)
            print("Std out", p)
'''


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
        data = reader.get_data(time_range=time_range, mzbins=bin_size)

    else:
        if scan_start > 0. and scan_end is not None and scan_end > 0.:
            scan_range = [scan_start, scan_end]
        else:
            scan_range = None

        data = reader.get_data(scan_range=scan_range, mzbins=bin_size)

    try:
        start = float(start)
        end = float(end)
    except:
        start = None
        end = None

    if start is not None and end is not None:
        data = ud.datachop(data, start, end)

    np.savetxt(msfile, data)
