#!/usr/bin/env python

use_mp_import = False

import os
from math import floor
import subprocess
from multiprocessing import Process, Queue, cpu_count
import sys
import platform
from pubsub import pub

from unidec_modules.tims_import_wizard import TagTypes as tt
from unidec_modules.tims_import_wizard import get_data_wrapper
try:
    from unidec_modules.waters_importer.Importer import WatersDataImporter as WDI
except Exception as e:
    print("Error importing Waters Importer, data_importer.py, ", e)

param_file_cache = {}


def auto_from_wizard(lines, exedir):
    '''
    Auto load from the import wizard
    '''
    # make list to convert items for importing
    conv = {'Filename': tt.FILE_NAME,
            'Full Path': tt.FILE_PATH,
            # 'Sample': tt.SAMPLE,
            # 'Description': tt.DESCRIPTION,
            'Pusher': tt.PUSHER,
            'Drift Voltage': tt.DRIFT_V,
            'Collision Voltage': tt.COLLISION_V,
            'Drift Pressure': tt.DRIFT_PRESSURE,
            'Drift Temperature': tt.TEMP,
            'Atom': tt.ATOM,
            'm/z Start': tt.START,
            'm/z End': tt.END,
            'm/z Bin Size': tt.BIN,
            'Calibration 1': tt.TCAL1,
            'Calibration 2': tt.TCAL2,
            'EDC': tt.EDC,
            'Function': tt.FUNCTION,
            'Scan Start': tt.SCAN_START,
            'Scan End': tt.SCAN_END}

    # Make empty job_list
    job_list = []

    # ADD EACH ITEM HERE
    for c, line in enumerate(lines):
        l = split_and_strip(line)

        # set mode type
        if c == 0:  # first line is header
            mode = l[1]

        elif c == 1:  # second line is headers too
            # map column headers
            mapper = {}
            for index, value in enumerate(l):
                mapper[index] = value.strip()

        # 3rd line and onwards are data line
        # process the line
        elif c > 1:
            parse = {'2D': None,
                     '1D': None,
                     tt.RANGE: None,
                     tt.BIN: "1",
                     tt.FILE_NAME: None,
                     tt.FILE_PATH: None,
                     tt.ATOM: None,
                     # tt.DESCRIPTION: None,
                     tt.COLLISION_V: None,
                     tt.TYPE: mode.lower(),  # 'ims mode'
                     tt.FUNCTION: '0',
                     'exedir': exedir}

            if mode.lower() == 'linear':
                parse[tt.DRIFT_PRESSURE] = None
                parse[tt.TEMP] = None
                parse[tt.DRIFT_V] = None
            elif mode.lower() == "ms":
                pass
            else:  # t-wave
                parse[tt.TCAL1] = None
                parse[tt.TCAL2] = None
                parse[tt.EDC] = None

            for index in mapper:
                if l[index] == "None":
                    l[index] = None
                if mapper[index] != '':
                    try:
                        parse[conv[mapper[index]]] = l[index]
                    except Exception as e:
                        print("Unable to import: ", index, mapper[index], l[index], e)

            # fix range
            parse[tt.RANGE] = (parse[tt.START], parse[tt.END])

            # ensure CV is float rather than int string
            if parse[tt.COLLISION_V] is not None:
                parse[tt.COLLISION_V] = float(parse[tt.COLLISION_V])

            if mode.lower() == 't-wave':
                if parse[tt.TCAL1] is not None:
                    parse[tt.TCAL1] = float(parse[tt.TCAL1])

                if parse[tt.TCAL2] is not None:
                    parse[tt.TCAL2] = float(parse[tt.TCAL2])

                if parse[tt.EDC] is not None:
                    parse[tt.EDC] = float(parse[tt.EDC])

            # alter atom from string to integar
            if not mode.lower() == "ms":
                atoms = {'He': 4, 'N2': 28}
                parse[tt.ATOM] = atoms[parse[tt.ATOM]]

            if use_mp_import:
                job_list.append(parse)
            else:
                process_from_wizard(**parse)

    if use_mp_import:
        # new style send for parallel processing
        mp_process_from_wizard(job_list)


def run_get_data(job_kwargs):
    try:
        start = float(job_kwargs[tt.START])
        end = float(job_kwargs[tt.END])
    except (ValueError, KeyError, AttributeError, TypeError):
        start = None
        end = None

    try:
        scan_start = int(float(job_kwargs[tt.SCAN_START]))
    except (ValueError, KeyError, AttributeError, TypeError):
        scan_start = 0

    try:
        scan_end = int(float(job_kwargs[tt.SCAN_END]))
    except (ValueError, KeyError, AttributeError, TypeError):
        scan_end = None
    # print job_kwargs[tt.FUNCTION]
    if int(job_kwargs[tt.FUNCTION]) != 1:
        job_kwargs[tt.FILE_NAME] = job_kwargs[tt.FILE_NAME] + '_%02d' % int(job_kwargs[tt.FUNCTION])
        # print job_kwargs[tt.FILE_NAME]
    if job_kwargs[tt.TYPE].lower() != 'ms':
        A, B, X, Y, C = get_data_wrapper.new_get_data(start,
                                                      end,
                                                      job_kwargs[tt.BIN],
                                                      job_kwargs[tt.FILE_PATH],
                                                      job_kwargs[tt.PUSHER],
                                                      job_kwargs[tt.FUNCTION],
                                                      scan_start,
                                                      scan_end,
                                                      dir=job_kwargs['exedir'])

        # A = [i[0] for i in X]
        # B = [i.sum() for i in C]

        # X is m/z
        # Y is drift time
        # C is height
        # A is unique m/z values
        # B is total intensity summed over drift tim
        # add point data here to kwarg dict
        job_kwargs['2D'] = (X, Y, C)
        job_kwargs['1D'] = (A, B)

        MakeUniDecConfig(job_kwargs)
    else:
        get_data_wrapper.new_get_data_MS(start,
                                         end,
                                         job_kwargs[tt.BIN],
                                         job_kwargs[tt.FILE_PATH],
                                         job_kwargs[tt.FUNCTION],
                                         scan_start,
                                         scan_end)
        # job_kwargs['1D'] = (A, B)


def MakeUniDecConfig(job_kwargs):
    raw_file = job_kwargs[tt.FILE_PATH]
    immsfile = os.path.splitext(raw_file)[0] + '_' + str(job_kwargs[tt.FUNCTION]) + '_imms.txt'
    dirnew = os.path.splitext(immsfile)[0] + "_unidecfiles"
    if not os.path.isdir(dirnew):
        os.mkdir(dirnew)

    header = os.path.splitext(job_kwargs[tt.FILE_NAME])[0] + '_' + str(job_kwargs[tt.FUNCTION]) + '_imms'
    filename = os.path.join(dirnew, header + "_conf.dat")

    if job_kwargs[tt.TYPE] == "linear":
        twaveflag = 0
    else:
        twaveflag = 1

    f = open(filename, 'w+')
    f.write("twaveflag " + str(twaveflag) + "\n")
    f.write("mzbins " + str(job_kwargs[tt.BIN]) + "\n")
    if job_kwargs[tt.PUSHER] is not None:
        f.write("pusher " + str(job_kwargs[tt.PUSHER]) + "\n")
    if job_kwargs[tt.ATOM] is not None:
        f.write("gasmass " + str(job_kwargs[tt.ATOM]) + "\n")
    if job_kwargs[tt.COLLISION_V] is not None:
        f.write("collision_v " + str(job_kwargs[tt.COLLISION_V]) + "\n")
    if twaveflag == 0:
        if job_kwargs[tt.TEMP] is not None:
            f.write("temp " + str(job_kwargs[tt.TEMP]) + "\n")
        if job_kwargs[tt.DRIFT_PRESSURE] is not None:
            f.write("pressure " + str(job_kwargs[tt.DRIFT_PRESSURE]) + "\n")
        if job_kwargs[tt.DRIFT_V] is not None:
            f.write("volt " + str(job_kwargs[tt.DRIFT_V]) + "\n")
    else:
        if job_kwargs[tt.EDC] is not None:
            f.write("edc " + str(job_kwargs[tt.EDC]) + "\n")
        if job_kwargs[tt.TCAL1] is not None:
            f.write("tcal1 " + str(job_kwargs[tt.TCAL1]) + "\n")
        if job_kwargs[tt.TCAL2] is not None:
            f.write("tcal2 " + str(job_kwargs[tt.TCAL2]) + "\n")

    f.close()


def process_from_wizard(**kwargs):
    '''
    Processes row in import wizard file
    and adds to data model.
    '''
    run_get_data(kwargs)
    pub.sendMessage('RAW DATA ADDED TO MODEL')


# uses multiprocessing module to import data across multiple cores
def mp_process_from_wizard(job_list):
    worker_queue = Queue()
    result_queue = Queue()

    for parse in job_list:
        worker_queue.put(parse)

    core_worker = 1
    try:
        core_worker = cpu_count() - 1
    except Exception as e:
        pass

    workers = [Process(target=data_import_worker, args=(worker_queue, result_queue)) for i in range(core_worker)]

    for each in workers:
        each.start()

    # we don't join yet as if data are large then worker process
    # does not buffer it all through to the queue
    # therefore have to wait before get (and therefore remove) from the queue
    # before joining

    result_count = 0
    # for job_kwargs in iter(result_queue.get_nowait, None):
    while True:
        if result_count == len(job_list):
            break
        try:
            job_kwargs = result_queue.get()
        except Exception as e:  # Queue.Empty:
            break
        pub.sendMessage('RAW DATA ADDED TO MODEL')
        result_count += 1

    for each in workers:
        each.join()


def data_import_worker(worker_queue, result_queue):
    while True:
        try:
            # false means it does not block waiting for more items
            job_kwargs = worker_queue.get(False)
        except Exception as e:  # Queue.Empty:
            break

        run_get_data(job_kwargs)
        result_queue.put(job_kwargs)
    '''
    try:
        #for job_kwargs in iter(worker_queue.get_nowait, None):
        while True:
            try:
                # false means it does not block waiting for more items
                job_kwargs = worker_queue.get(False)
            except Exception, e: # Queue.Empty:
                break

            run_get_data(job_kwargs)
            result_queue.put(job_kwargs)

    except Exception, e:
        print 'Failed Import'
        print e, e.message
        result_queue.put(job_kwargs)

    return
    '''


def split_and_strip(line, delimiter=','):
    tmp = []
    for l in line.split(delimiter):
        if '\n' in l and l != '\n':
            l = l.rsplit('\n')[0]

        tmp.append(l)

    return tmp


def parse_file(file_path, exp_type='linear', collision=None, debug=False, dir=None):
    file_name = os.path.split(file_path)[1]
    exedir = dir
    # check a data folder
    if file_name.endswith('.raw'):
        des = header_desc(file_path)

        # encapsulate
        out = {tt.FILE_PATH: file_path,
               tt.FILE_NAME: file_name,
               tt.DRIFT_PRESSURE: None,
               tt.PUSHER: None,
               tt.TEMP: None,
               tt.ATOM: 'He',
               tt.DRIFT_V: None,
               tt.TCAL1: None,
               tt.TCAL2: None,
               tt.EDC: None,
               tt.CONE: None,
               tt.EXTRACT: None,
               tt.COLLISION_V: None,
               tt.START: None,
               tt.END: None,
               tt.BIN: "1",
               tt.TYPE: exp_type,
               # tt.DESCRIPTION: None,
               # tt.SAMPLE: None,
               tt.FUNCTION: '1',
               tt.SCAN_START: None,
               tt.SCAN_END: None}

        # only for t-wave
        if out[tt.TYPE] == 't-wave':
            # wave_velocity = search_extern(file_path, 'IMS Wave Velocity (m/s)')
            # if wave_velocity is not None:
            #    out[tt.WAVE_VELOCITY] = float(wave_velocity)

            # wave_height = search_extern(file_path, 'IMS Wave Height (V)')
            # if wave_height is not None:
            #    out[tt.WAVE_HEIGHT] = float(wave_height)

            # out[tt.COLLISION_V] = search_extern(file_path, 'Trap Collision Energy')
            edc = search_extern(file_path, 'EDC Delay Coefficient')
            if edc is not None:
                out[tt.EDC] = float(edc)

            out[tt.ATOM] = 'N2'

        # only for linear mode
        if out[tt.TYPE] == 'linear':
            try:
                # print des.split()
                out[tt.DRIFT_PRESSURE] = str(float(des.split()[2]))
                out[tt.TEMP] = str(float(des.split()[4]))
            except Exception as e:
                # modified for jon's format
                des = des.replace("=", " ")
                weird_des = des.split()

                for c, item in enumerate(weird_des):
                    # print item
                    # get temp
                    if item.lower() == 'temp' or item.lower() == 't':
                        if ',' in weird_des[c + 1]:
                            x = float(weird_des[c + 1].split(',')[0])
                        else:
                            try:
                                x = float(weird_des[c + 1])
                            except (ValueError, IndexError):
                                x = ''
                        out[tt.TEMP] = x

                    # get drift pressure
                    if item.lower() == 'pres' or item.lower() == 'p':
                        if ',' in weird_des[c + 1]:
                            x = float(weird_des[c + 1].split(',')[0])
                        else:
                            x = float(weird_des[c + 1])
                        out[tt.DRIFT_PRESSURE] = x

                    # get collision voltage
                    if item.lower() == 'trap':
                        if ',' in weird_des[c + 1]:
                            x = float(weird_des[c + 1].split(',')[0])
                        else:
                            x = float(weird_des[c + 1])
                        out[tt.COLLISION_V] = x

                pass

            out[tt.DRIFT_V] = search_extern(file_path, 'Transfer Collision Energy')

        if out[tt.TYPE] != 'ms':
            # try to grab the pusher frequency from the stat code
            # pusher = get_stat_code(file_path, 76, dir=exedir)
            pusher = get_stat_name(file_path, "Transport RF")
            if pusher is not None:
                out[tt.PUSHER] = floor((1. / pusher) * 1000 * 1000)
        else:
            out[tt.FUNCTION] = "0"

        # try to grab the collision voltage from the stat code
        # cv = get_stat_code(file_path, 62, dir=exedir)
        cv = get_stat_name(file_path, "Collision Energy")
        if cv is not None:
            out[tt.COLLISION_V] = cv
        if collision != None:
            out[tt.COLLISION_V] = str(collision)

        # get some general instrument parameters
        out[tt.CONE] = search_extern(file_path, 'Sampling Cone')
        out[tt.EXTRACT] = search_extern(file_path, 'Extraction Cone')

        # make a good guess at the pusher frequency
        # out[tt.PUSHER] = ''
        # pusher_interval = search_extern(file_path, 'Pusher Interval', split_skip=True)
        # if pusher_interval is not None:
        #    out[tt.PUSHER] = pusher_interval

        # get range of m/z in file
        start, end = GetStartEndMass(file_path)[1:]
        out[tt.START] = start
        out[tt.END] = end

        if debug:
            for key in out:
                print(key, out[key])

        return out
    return None


def GetLines(InputFileName):
    # Read input file
    # Creates an array containing elements for each lines of the input file.
    if os.path.isfile(InputFileName):
        InputFile = open(InputFileName, 'r')  # Open file to read
        Lines = InputFile.read().splitlines()  # string of lines from pdb file
        InputFile.close()
    else:
        print("Could not find file: ", InputFileName)
        return []
    return Lines


def GetStartEndMass(RawFile):
    ce, s, e = None, None, None
    if (platform.uname()[0] == 'Windows'):
        ParamFile = RawFile + "\_extern.inf"
    else:
        ParamFile = RawFile + "/_extern.inf"

    if (os.path.exists(ParamFile) == 0):
        sys.stdout.write('Cannot find param file %s\n' % ParamFile)
        # sys.exit(100)

    for line in GetLines(ParamFile):
        if line.find('Start Mass') >= 0:
            s = float(line.split()[2])
        if line.find('End Mass') >= 0:
            try:
                e = float(line.split()[2])
            except ValueError:
                e = float(line.split()[3])
        if line.startswith('Trap Collision Energy'):
            ce = float(line.split()[3])
    return ce, s, e


'''
def get_stat_code_old(raw_file, stat_code, dir=None):
    if dir is None:
        path = os.getcwd()
    else:
        path = dir
    rawreader_path = os.path.join(path, 'rawreadertim.exe')

    if not os.path.isfile(rawreader_path):
        rawreader_path = os.path.join(path, "unidec_bin", 'rawreadertim.exe')
        if not os.path.isfile(rawreader_path):
            print("Unable to find rawreadertim.exe")
            print(rawreader_path)
            return None
    # default/hard coded to function 0, odd, 0 is for rawreader and 1 is for CDCReader
    process = subprocess.Popen([rawreader_path, raw_file, '--fn=0', '--stat=' + str(stat_code)],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    (stdout, stderr) = process.communicate()
    try:
        param = float(stdout)
        if param > 0.:
            return param
    except ValueError:
        return None
'''


def get_stat_code(raw_file, stat_code, dir=None):
    param = WDI(raw_file, do_import=False).get_stat_code(stat_code)
    # print(param)
    try:
        param = float(param)
        # print(param)
        return param
    except Exception:
        return None


def get_stat_name(raw_file, stat_name):
    param = WDI(raw_file, do_import=False).get_stat_name(stat_name)
    # print(param)
    try:
        param = float(param)
        # print(param)
        return param
    except Exception:
        return None


def header_desc(file_path):
    des_path = os.path.join(file_path, '_HEADER.TXT')
    if os.path.isfile(des_path):
        # get description from header
        f = open(des_path, 'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            if line.startswith('$$ Sample Description:'):
                return line.split('$$ Sample Description:')[1].strip()

    else:
        return None


def search_extern(file_path, search_string, split_skip=False):
    # note that split_skip will go +1 on what is returned to get around non-ASCII unit issue
    # e.g. Pusher Interval (uS) where u is Greek mu, search for "Pusher Interval" and then set split_skip = True
    '''
    Pass in each string left of target value.
    An example,
        target
        ------
        Trap Wave Velocity (m/s)	300
        search string='Trap Wave Velocity (m/s)'
        returns :: 300
    '''
    global param_file_cache

    param_path = os.path.join(file_path, '_extern.inf')
    if os.path.isfile(param_path):
        if param_path in param_file_cache:
            size, mtime, lines = param_file_cache[param_path]
            try:
                stat = os.stat(param_path)

            except os.error:
                del param_file_cache[param_path]

            if size != stat.st_size or mtime != stat.st_mtime:
                del param_file_cache[param_path]

        else:
            # load the file and add it to the cache
            try:
                stat = os.stat(param_path)

            except os.error:
                pass

            f = open(param_path, 'r')
            lines = f.readlines()
            f.close()
            param_file_cache[param_path] = (stat.st_size, stat.st_mtime, lines)

        for line in lines:
            if line.startswith(search_string):
                if split_skip:
                    index = 1

                else:
                    index = 0

                tmp = line.split(search_string)[1].split()[index]
                try:
                    tmp = str(float(tmp))

                except ValueError:
                    tmp = None

                return tmp


if __name__ == '__main__':
    pass
