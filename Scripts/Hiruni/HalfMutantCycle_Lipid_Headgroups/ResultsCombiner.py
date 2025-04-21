"""
Methods for parsing Agilent Chemstation files.

"""

import os
import struct
import numpy as np
from lxml import etree

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt


class DataFile:
    """
    Class representing a chromatogram data file.

    Args:
        path (str): Path of the file.
        detector (str): Detector for the file.
        xlabels (numpy.ndarray): 1D array with retention times (in minutes).
        ylabels (numpy.ndarray): 1D array with y-axis labels
            (e.g. mz, wavelength).
        data (numpy.ndarray): 2D array with data values (e.g. intensity).
        metadata (dict): Metadata for the file.

    Attributes:
        name (str): Name of the file.
        detector (str): Name of the detector. Options: UV, MS, FID, CAD, ELSD.
        xlabels (numpy.ndarray): 1D array with retention times (in minutes).
        ylabels (numpy.ndarray): 1D array with y-axis labels
            (e.g. mz, wavelength).
        data (numpy.ndarray): 2D array with data values (e.g. intensity).
            The rows correspond to the retention times and the columns
            correspond to the y-axis labels.
        metadata (dict): Depends on the vendor and file format.

    """

    def __init__(self, path, detector, xlabels, ylabels, data, metadata):

        if not isinstance(path, str) or \
                not detector in {'UV', 'MS', 'FID', 'CAD', 'ELSD', None} or \
                not isinstance(xlabels, np.ndarray) or xlabels.ndim != 1 or \
                not isinstance(ylabels, np.ndarray) or ylabels.ndim != 1 or \
                not isinstance(data, np.ndarray) or data.ndim != 2 or \
                not isinstance(metadata, dict):
            raise Exception("Wrong argument parameters for DataFile.")

        self.name = os.path.basename(path)
        self.detector = detector
        self.xlabels = xlabels
        self.ylabels = ylabels
        self.data = data
        self.metadata = metadata
        warnings.filterwarnings("ignore", category=FutureWarning)

    def __repr__(self):
        return f"{self.name}"

    def get_info(self):
        """
        Returns a string summary of the DataFile.

        """
        return f"\n{'-' * len(self.name)}\n" \
               f"{self.name}\n" \
               f"{'-' * len(self.name)}\n" \
               f"Detector: {self.detector}\n" \
               f"Xlabels: {self.xlabels}\n" \
               f"Ylabels: {self.ylabels}\n" \
               f"Data: {self.data}\n" \
               f"Metadata: {self.metadata}\n"

    def extract_traces(self, labels=None):
        """
        Extracts data corresponding to the specified :code:`labels`.

        Raises an exception if any :code:`labels` are invalid.

        Args:
            labels (int/float/list, optional): Ylabel(s) to extract.

        Returns:
            2D numpy array containing data for the specified ylabel(s).
            The rows correspond to the ylabels and the columns corrrespond \
                to the retention times.

        """

        if labels is None:
            return self.data.T

        if isinstance(labels, int) or \
                isinstance(labels, float) or \
                labels == '':
            labels = [labels]

        if not isinstance(labels, list):
            raise Exception("Invalid type for labels.")

        for label in np.array(labels):
            if label not in self.ylabels:
                raise Exception(f"Label {label} not in {self.name}.")

        indices = np.searchsorted(self.ylabels, labels)
        traces = self.data[:, indices].T

        return traces

    def export_csv(self, filename, labels=None, delim=','):
        """
        Outputs a CSV containing data for the specified :code:`labels`.

        Args:
            filename (str): Filename for the output CSV.
            labels (int/float/list, optional):   Ylabel(s) to export.
            delim (str, optional): Delimiter used in the output CSV.

        """
        f = open(filename, 'w+')
        f.write(self.to_csvstr(labels, delim))
        f.close()

    def to_csvstr(self, labels=None, delim=','):
        """
        Returns a string representation of a CSV containing data \
            for the specified :code:`labels`.

        Args:
            labels (int/float/list, optional): Ylabel(s) to return.
            delim (str, optional): Delimiter used in the CSV representation.

        """
        str_traces_tp = self.extract_traces(labels).T.astype(str)

        if isinstance(labels, int) or \
                isinstance(labels, float) or \
                labels == '':
            labels = [labels]

        if labels is None:
            labels = self.ylabels

        str_labels = [str(label) for label in labels]

        csvstr = f"RT (min){delim}{f'{delim}'.join(str_labels)}\n"
        for i in range(self.xlabels.size):
            csvstr += f"{self.xlabels[i]}{delim}" \
                      f"{f'{delim}'.join(str_traces_tp[i])}\n"
        return csvstr

    def plot(self, label, **kwargs):
        """
        Shows a basic matplotlib plot for the specified :code:`label`.

        Args:
            label (int/float): Ylabel to be plotted.
            **kwargs (optional): Keyword arguments for matplotlib.

        """
        import matplotlib.pyplot as plt
        plt.plot(self.xlabels, self.extract_traces(label).T, **kwargs)
        plt.show()

    def get_data(self):
        return np.transpose([self.xlabels, self.data[:,0]])

"""
MAIN PARSING METHODS

"""


def parse_allfiles(path, prec=0, requested_files=None):
    """
    Finds and parses Agilent Chemstation data files \
        with a .ch, .uv, or .ms extension from a .D directory.

    Args:
        path (str): Path to the .D directory.
        prec (int, optional): Number of decimals to round mz values.
        requested_files (list, optional): List of filenames to parse.

    Returns:
        List with a DataFile for each parsed data file.

    """
    datafiles = []
    for name in os.listdir(path):
        if requested_files and name.lower() not in requested_files:
            continue
        datafile = parse_file(os.path.join(path, name), prec)
        if datafile:
            datafiles.append(datafile)
    return datafiles


def parse_file(path, prec=0):
    """
    Parses an Agilent Chemstation data file.

    Supported extensions are .ch, .uv, and .ms.

    Args:
        path (str): Path to the data file.
        prec (int, optional): Number of decimals to round mz values.

    Returns:
        DataFile representing the file, if it can be parsed. Otherwise, None.

    """
    ext = os.path.splitext(path)[1].lower()
    if ext == '.ch':
        return parse_ch(path)
    elif ext == '.uv':
        return parse_uv(path)
    elif ext == '.ms':
        return parse_ms(path, prec)
    return None


"""
.ch PARSING METHODS

"""


def parse_ch(path):
    """
    Parses an Agilent .ch file.

    These files contain data from a FID, CAD, ELSD, or UV channel. \
    Files that contain FID data have a different format than other .ch files.

    This method calls the appropriate subroutine by file format.

    Args:
        path (str): Path to the .ch file.

    Returns:
        DataFile with data from a channel, if the file can be parsed. \
            Otherwise, None.

    """
    with open(path, 'rb') as f:
        head = read_string(f, offset=0, gap=1)
        if head == '179':
            return parse_ch_fid(path)
        elif head == '130' or head == '30':
            return parse_ch_other(path, head)
        return None


def parse_ch_fid(path):
    """
    Parses an Agilent .ch file with FID channel data.

    This method should not be called directly. Use :obj:`parse_ch` instead.

    Learn more about this file format :ref:`here <ch_fid>`.

    Args:
        path (str): Path to the .ch file with FID data.

    Returns:
        DataFile with FID data, if the file can be parsed. Otherwise, None.

    """
    data_offsets = {
        'num_times': 0x116,
        'scaling_factor': 0x127C,
        'data_start': 0x1800
    }

    f = open(path, 'rb')
    raw_bytes = f.read()

    # Extract the number of retention times.
    f.seek(data_offsets['num_times'])
    num_times = struct.unpack(">I", f.read(4))[0]
    if num_times == 0:
        return None

    # Compute retention times using the first and last times.
    start_time = struct.unpack(">f", f.read(4))[0]
    end_time = struct.unpack(">f", f.read(4))[0]
    delta_time = (end_time - start_time) / (num_times - 1)
    times = np.arange(start_time, end_time + 1e-3, delta_time)

    # Extract the raw data values.
    data = np.ndarray(
        num_times, '<d', raw_bytes, data_offsets['data_start'], 8)
    data = data.copy().reshape(-1, 1)

    # Convert times into minutes.
    times /= 60000

    # Scale the absorbances.
    f.seek(data_offsets['scaling_factor'])
    scaling_factor = struct.unpack('>d', f.read(8))[0]
    data *= scaling_factor

    # No ylabel for FID data.
    ylabels = np.array([''])

    # Extract metadata from file header.
    metadata_offsets = {
        'notebook': 0x35A,
        'date': 0x957,
        'method': 0xA0E,
        'instrument': 0xC11,
        'unit': 0x104C,
        'signal': 0x1075
    }
    metadata = read_header(f, metadata_offsets)
    f.close()

    return DataFile(path, 'FID', times, ylabels, data, metadata)


def parse_ch_other(path, head):
    """
    Parses an Agilent .ch file with CAD, ELSD, or UV channel data.

    This method should not be called directly. Use :obj:`parse_ch` instead.

    IMPORTANT: ELSD data may be mistakenly labeled as CAD on rare occasions. Users may need to make this distinction on their own when decoding Agilent CAD or ELSD data.

    Learn more about this file format :ref:`here <ch_other>`.

    Args:
        path (str): Path to the .ch file with UV, CAD, or ELSD data.

    Returns:
        DataFile with CAD, ELSD, or UV data, if parsable. Otherwise, None.

    """
    if head == '130':
        data_offsets = {
            'time_range': 0x11A,
            'scaling_factor': 0x127C,
            'data_start': 0x1800
        }
        metadata_offsets = {
            'notebook': 0x35A,
            'date': 0x957,
            'method': 0xA0E,
            'instrument': 0xC11,
            'unit': 0x104C,
            'signal': 0x1075
        }
        gap = 2
    elif head == '30':
        data_offsets = {
            'time_range': 0x11A,
            'scaling_factor': 0x284,
            'data_start': 0x400
        }
        metadata_offsets = {
            'notebook': 0x18,
            'date': 0xB2,
            'method': 0xE4,
            'instrument': 0xDA,
            'unit': 0x244,
            'signal': 0x254
        }
        gap = 1
    else:
        return None

    f = open(path, 'rb')
    byte_unpack = struct.Struct('>B').unpack
    short_unpack = struct.Struct('>h').unpack
    int_unpack = struct.Struct('>i').unpack

    # Extract the raw data values.
    # Count the total number of retention times.
    f.seek(data_offsets['data_start'])
    absorbances = []
    absorb_accum = 0
    while True:
        # If the segment header is invalid, stop reading.
        head = byte_unpack(f.read(1))[0]
        if head != 0x10:
            break
        num_times_seg = byte_unpack(f.read(1))[0]

        # If the next short is equal to -0x8000
        #     then the next absorbance value is the next integer.
        # Otherwise, the short is a delta from the last absorbance value.
        for _ in range(num_times_seg):
            check_int = short_unpack(f.read(2))[0]
            if check_int == -0x8000:
                absorb_accum = int_unpack(f.read(4))[0]
            else:
                absorb_accum += check_int
            absorbances.append(absorb_accum)

    # Process the extracted values.
    # If no values are extracted, this file is invalid.
    data = np.array(absorbances)
    num_times = data.size
    if num_times == 0:
        return None

    # Calculate retention times using the first and last times.
    f.seek(data_offsets['time_range'])
    start_time, end_time = struct.unpack('>ii', f.read(8))
    delta_time = (end_time - start_time) / (num_times - 1)
    times = np.arange(start_time, end_time + 1e-3, delta_time)

    # Convert time to minutes
    times /= 60000

    # Process the absorbances.
    f.seek(data_offsets['scaling_factor'])
    scaling_factor = struct.unpack('>d', f.read(8))[0]
    data = data.reshape(-1, 1) * scaling_factor

    # Read file metadata.

    metadata = read_header(f, metadata_offsets, gap=gap)
    f.close()

    # Determine the detector and ylabels using metadata.
    detector = None
    ylabel = ''
    signal = metadata['signal']
    if '=' in signal:
        ylabel = signal.split('=')[1].split(',')[0]
        detector = 'UV'
    elif 'ADC' in signal:
        detector = 'ELSD' if 'CHANNEL' in signal else 'CAD'
    ylabels = np.array([ylabel])

    return DataFile(path, detector, times, ylabels, data, metadata)


"""
.uv PARSING METHODS

"""


def decode_uv_delta(f, data_offsets, num_times, num_wavelengths):
    uint_unpack = struct.Struct('<I').unpack
    int_unpack = struct.Struct('<i').unpack
    short_unpack = struct.Struct('<h').unpack

    f.seek(data_offsets["data_start"])
    times = np.empty(num_times, dtype=np.uint32)
    data = np.empty((num_times, num_wavelengths), dtype=np.int64)
    for i in range(num_times):
        f.read(4)
        times[i] = uint_unpack(f.read(4))[0]
        f.read(14)
        # If the next short is equal to -0x8000
        #     then the next absorbance value is the next integer.
        # Otherwise, the short is a delta from the last absorbance value.
        absorb_accum = 0
        for j in range(num_wavelengths):
            check_int = short_unpack(f.read(2))[0]
            if check_int == -0x8000:
                absorb_accum = int_unpack(f.read(4))[0]
            else:
                absorb_accum += check_int
            data[i, j] = absorb_accum

    return times, data


def decode_uv_array(f, data_offsets, num_times, num_wavelengths):
    uint_unpack = struct.Struct('<I').unpack

    f.seek(data_offsets["data_start"])
    times = np.empty(num_times, dtype=np.uint32)
    data = np.empty((num_times, num_wavelengths), dtype=np.float64)
    for i in range(num_times):
        f.read(4)
        times[i] = uint_unpack(f.read(4))[0]
        f.read(14)
        for j in range(num_wavelengths):
            data[i, j] = struct.unpack('<d', f.read(8))[0]

    return times, data


def parse_uv(path):
    """
    Parses an Agilent .uv file.

    These files contain UV spectra.

    Learn more about this file format :ref:`here <uv>`.

    Args:
        path (str): Path to the Agilent .uv file.

    Returns:
        DataFile with UV data, if the file can be parsed. Otherwise, None.

    """

    f = open(path, 'rb')
    uint_unpack = struct.Struct('<I').unpack
    int_unpack = struct.Struct('<i').unpack
    short_unpack = struct.Struct('<h').unpack

    # Validate file header.
    head = read_string(f, 0, gap=1)

    if head == '131':
        data_offsets = {
            'num_times': 0x116,
            'scaling_factor': 0xC0D,
            'data_start': 0x1000
        }
        metadata_offsets = {
            "notebook": 0x35A,
            "date": 0x957,
            "method": 0xA0E,
            "unit": 0xC15,
            "signal": 0xC40,
            "vialpos": 0xFD7
        }
        file_type = read_string(f, 347, gap=2)
        if file_type.startswith('LC'):
            decode = decode_uv_delta
        elif file_type.startswith('OL'):
            decode = decode_uv_array
        else:
            return None
        gap = 2
    elif head == '31':
        data_offsets = {
            'num_times': 0x116,
            'scaling_factor': 0x13E,
            'data_start': 0x200
        }
        metadata_offsets = {
            "notebook": 0x18,
            "date": 0xB2,
            "method": 0xE4,
            "unit": 0x146
        }
        decode = decode_uv_delta
        gap = 1
    else:
        f.close()
        return None

    # Extract the number of retention times.
    f.seek(data_offsets["num_times"])
    num_times = struct.unpack(">I", f.read(4))[0]
    # If there are none, the file may be a partial.
    if num_times == 0:
        f.close()
        return parse_uv_partial(path)

    # Compute the wavelengths by taking the range from
    #     the header of the first data segment
    f.seek(data_offsets["data_start"] + 0x8)
    start_wlen, end_wlen, delta_wlen = \
        tuple(num // 20 for num in struct.unpack("<HHH", f.read(6)))
    wavelengths = np.arange(start_wlen, end_wlen + 1, delta_wlen)
    num_wavelengths = wavelengths.size

    # Extract the retention times and absorbances from each data segment.
    times, data = decode(f, data_offsets, num_times, num_wavelengths)

    # Covert times to minutes.
    times = times / 60000

    # Scale the absorbances.
    f.seek(data_offsets['scaling_factor'])
    scaling_factor = struct.unpack('>d', f.read(8))[0]
    data = data * scaling_factor

    # Read file metadata.

    metadata = read_header(f, metadata_offsets, gap=gap)
    f.close()

    return DataFile(path, 'UV', times, wavelengths, data, metadata)


def parse_uv_partial(path):
    """
    Parses a partial Agilent .uv file.

    Learn more about this file format :ref:`here <uv>`.

    Args:
        path (str): Path to the partial .uv file.

    Returns:
        DataFile with UV data, if the file can be parsed. Otherwise, None.

    """
    data_offsets = {
        'num_times': 0x116,
        'scaling_factor': 0xC0D,
        'data_start': 0x1000
    }

    f = open(path, 'rb')
    uint_unpack = struct.Struct('<I').unpack
    int_unpack = struct.Struct('<i').unpack
    short_unpack = struct.Struct('<h').unpack

    # Compute the wavelengths by taking the range from
    #     the header of the first data segment.
    # If this process fails, then the file is not a partial.
    f.seek(data_offsets["data_start"] + 0x8)
    try:
        start_wlen, end_wlen, delta_wlen = \
            tuple(num // 20 for num in struct.unpack("<HHH", f.read(6)))
        wavelengths = np.arange(start_wlen, end_wlen + 1, delta_wlen)
    except Exception:
        return None

    # Extract the retention times and absorbances from each data segment.
    f.seek(data_offsets['data_start'])
    times = []
    absorbances = []
    while True:
        try:
            f.read(4)
            time = uint_unpack(f.read(4))[0]
            times.append(time)
            f.read(14)
            # If the next short is equal to -0x8000
            #     then the next absorbance value is the next integer.
            # Otherwise, the short is a delta from the last absorbance value.
            absorb_accum = 0
            for _ in range(wavelengths.size):
                check_int = short_unpack(f.read(2))[0]
                if check_int == -0x8000:
                    absorb_accum = int_unpack(f.read(4))[0]
                else:
                    absorb_accum += check_int
                absorbances.append(absorb_accum)
        except Exception:
            break

    # Process the extracted values.
    times = np.array(times) / 60000
    data = np.array(absorbances).reshape((times.size, wavelengths.size))

    # Scale the absorbances.
    f.seek(data_offsets['scaling_factor'])
    scaling_factor = struct.unpack('>d', f.read(8))[0]
    data = data * scaling_factor

    # Read file metadata.
    metadata_offsets = {
        "notebook": 0x35A,
        "date": 0x957,
        "method": 0xA0E,
        "unit": 0xC15,
        "signal": 0xC40,
        "vialpos": 0xFD7
    }
    metadata = read_header(f, metadata_offsets)
    f.close()

    return DataFile(path, 'UV', times, wavelengths, data, metadata)


"""
.ms PARSING METHODS

"""


def parse_ms(path, prec=0):
    """
    Parses an Agilent .ms file.

    These files contain MS spectra and SIM.

    Learn more about this file format :ref:`here <ms>`.

    Args:
        path (str): Path to Agilent .ms file.
        prec (int, optional): Number of decimals to round mz values.

    Returns:
        DataFile with MS data, if the file can be parsed. Otherwise, None.

    """
    data_offsets = {
        'type': 0x4,
        'data_start': 0x10A,
        'lc_num_times': 0x116,
        'gc_num_times': 0x142
    }

    f = open(path, 'rb')
    short_unpack = struct.Struct('>H').unpack
    int_unpack = struct.Struct('>I').unpack

    # Validate file header.
    # If invalid, the file may be a partial.
    head = int_unpack(f.read(4))[0]
    if head != 0x01320000:
        f.close()
        return parse_ms_partial(path, prec)

    # Determine the type of .ms file based on header.
    # Read the number of retention times from different offsets by type.
    type_ms = read_string(f, data_offsets['type'], 1)
    if type_ms == "MSD Spectral File":
        f.seek(data_offsets['lc_num_times'])
        num_times = int_unpack(f.read(4))[0]
    else:
        f.seek(data_offsets['gc_num_times'])
        num_times = struct.unpack('<I', f.read(4))[0]

    # Go to the data start offset.
    f.seek(data_offsets['data_start'])
    f.seek(short_unpack(f.read(2))[0] * 2 - 2)

    # Extract retention times and data pair counts for each time.
    # Store the bytes holding mz-intensity pairs.
    times = np.empty(num_times, dtype=np.uint32)
    pair_counts = np.zeros(num_times, dtype=np.uint16)
    pair_bytearr = bytearray()
    for i in range(num_times):
        f.read(2)
        times[i] = int_unpack(f.read(4))[0]
        f.read(6)
        pair_counts[i] = short_unpack(f.read(2))[0]
        f.read(4)
        pair_bytes = f.read(pair_counts[i] * 4)
        pair_bytearr.extend(pair_bytes)
        f.read(10)

    # Minor processing on the extracted data.
    raw_bytes = bytes(pair_bytearr)
    times = times / 60000
    total_paircount = np.sum(pair_counts)

    # Calculate the mz values.
    mzs = np.ndarray(total_paircount, '>H', raw_bytes, 0, 4)
    mzs = np.round(mzs / 20, prec)

    # Calculate the intensity values.
    int_encs = np.ndarray(total_paircount, '>H', raw_bytes, 2, 4)
    int_heads = int_encs >> 14
    int_tails = int_encs & 0x3fff
    int_values = np.multiply(8 ** int_heads, int_tails, dtype=np.uint32)
    del int_encs, int_heads, int_tails, raw_bytes

    # Make the array of `ylabels` with mz values.
    ylabels = np.unique(mzs)
    ylabels.sort()

    # Make the `data` array with intensities.
    mz_indices = np.searchsorted(ylabels, mzs)
    data = np.zeros((num_times, ylabels.size), dtype=np.uint32)
    cur_index = 0
    for i in range(num_times):
        stop_index = cur_index + pair_counts[i]
        np.add.at(
            data[i],
            mz_indices[cur_index:stop_index],
            int_values[cur_index:stop_index])
        cur_index = stop_index
    del mz_indices, mzs, int_values, pair_counts

    # Read file metadata.
    metadata_offsets = {
        'date': 0xB2,
        'method': 0xE4
    }
    metadata = read_header(f, metadata_offsets, 1)
    f.close()

    return DataFile(path, 'MS', times, ylabels, data, metadata)


def parse_ms_partial(path, prec=0):
    """
    Parses a partial Agilent .ms file.

    IMPORTANT: This method only supports LC .ms partials.

    Learn more about this file format :ref:`here <ms>`.

    Args:
        path (str): Path to the partial .ms file.
        prec (int, optional): Number of decimal to round mz values.

    Returns:
        DataFile with MS data, if the file can be parsed. Otherwise, None.

    """
    f = open(path, 'rb')
    short_unpack = struct.Struct('>H').unpack
    int_unpack = struct.Struct('>I').unpack

    # Partial .ms files do not store the start offset.
    # Shallow validation of filetype by checking that offset is null.
    f.seek(0x10A)
    if short_unpack(f.read(2))[0] != 0:
        f.close()
        return None

    # The start offset for data in .ms files is technically variable,
    #     but it has been constant for every .ms file we have tested.
    # Since the start offset is not stored in partials, this code uses that
    #     "constant" common starting offset. It may not work in all cases.
    f.seek(0x2F2)

    # Extract retention times and data pair counts for each time.
    # Store the bytes holding mz-intensity pairs.
    times = []
    pair_counts = []
    pair_bytearr = bytearray()
    while True:
        try:
            f.read(2)
            time = int_unpack(f.read(4))[0]
            f.read(6)
            pair_count = short_unpack(f.read(2))[0]
            f.read(4)
            pair_bytes = f.read(pair_count * 4)
            f.read(10)
            times.append(time)
            pair_counts.append(pair_count)
            pair_bytearr.extend(pair_bytes)
        except Exception:
            break

    # Minor processing on the extracted data.
    raw_bytes = bytes(pair_bytearr)
    times = np.array(times) / 60000
    pair_counts = np.array(pair_counts)
    num_times = times.size
    total_paircount = np.sum(pair_counts)

    # Calculate the mz values.
    mzs = np.ndarray(total_paircount, '>H', raw_bytes, 0, 4)
    mzs = np.round(mzs / 20, prec)

    # Calculate the intensity values.
    int_encs = np.ndarray(total_paircount, '>H', raw_bytes, 2, 4)
    int_heads = int_encs >> 14
    int_tails = int_encs & 0x3fff
    int_values = np.multiply(8 ** int_heads, int_tails, dtype=np.uint32)
    del int_encs, int_heads, int_tails, raw_bytes

    # Make the array of `ylabels` with mz values.
    ylabels = np.unique(mzs)
    ylabels.sort()

    # Fill the `data` array with intensities.
    mz_indices = np.searchsorted(ylabels, mzs)
    data = np.zeros((num_times, ylabels.size), dtype=np.uint32)
    cur_index = 0
    for i in range(num_times):
        stop_index = cur_index + pair_counts[i]
        np.add.at(
            data[i],
            mz_indices[cur_index:stop_index],
            int_values[cur_index:stop_index])
        cur_index = stop_index
    del mz_indices, mzs, int_values, pair_counts

    # Read file metadata.
    metadata_offsets = {
        'date': 0xB2,
        'method': 0xE4
    }
    metadata = read_header(f, metadata_offsets, 1)
    f.close()

    return DataFile(path, 'MS', times, ylabels, data, metadata)


""" 
FILE METADATA PARSING METHODS

"""


def read_header(f, offsets, gap=2):
    """
    Extracts metadata from the header of an Agilent data file.

    Args:
        f (_io.BufferedReader): File opened in 'rb' mode.
        offsets (dict): Dictionary mapping properties to file offsets.
        gap (int): Distance between two adjacent characters.

    Returns:
        Dictionary containing metadata as string key-value pairs.

    """
    metadata = {}
    for key, offset in offsets.items():
        string = read_string(f, offset, gap)
        if string:
            metadata[key] = string
    return metadata


def read_string(f, offset, gap=2):
    """
    Extracts a string from the specified offset.

    This method is primarily useful for retrieving metadata.

    Args:
        f (_io.BufferedReader): File opened in 'rb' mode.
        offset (int): Offset to begin reading from.
        gap (int): Distance between two adjacent characters.

    Returns:
        String at the specified offset in the file header.

    """
    f.seek(offset)
    str_len = struct.unpack("<B", f.read(1))[0] * gap
    try:
        return f.read(str_len)[::gap].decode().strip()
    except Exception:
        return ""


""" 
DIRECTORY METADATA PARSING METHODS 

"""


def parse_metadata(path, datafiles):
    """
    Parses Agilent metadata at the directory level.

    First, the DataFiles are checked for date and vial position metadata.

    Then, several files are scanned for the vial position. \
        This method can look inside the AcqData directory, which may be \
        misleading because this method resides in the Chemstation module.

    Args:
        path (str): Path to the directory.
        datafiles (list): List of DataFile objects.

    Returns:
        Dictionary containing directory metadata.

    """
    metadata = {}
    metadata['vendor'] = "Agilent"

    # Scan each DataFile for the date and vial position.
    # These may be stored in multiple files but the values are constant.
    # In MS files, the time may be saved in a different format.
    for datafile in datafiles:
        if 'date' not in metadata and 'date' in datafile.metadata:
            metadata['date'] = datafile.metadata['date']
        if 'vialpos' not in metadata and 'vialpos' in datafile.metadata:
            metadata['vialpos'] = datafile.metadata['vialpos']
    if 'date' in metadata and 'vialpos' in metadata:
        return metadata

    # Scan certain files for the vial position.
    dircontents = set(os.listdir(path))

    # sequence.acam_
    if "sequence.acam_" in dircontents:
        vialnum = get_xml_vialnum(os.path.join(path, "sequence.acam_"))
        if vialnum:
            metadata['vialpos'] = vialnum
            return metadata

    # sample.acaml
    if "sample.acaml" in dircontents:
        vialnum = get_xml_vialnum(os.path.join(path, "sample.acaml"))
        if vialnum:
            metadata['vialpos'] = vialnum
            return metadata

    # AcqData/sample_info.xml
    if "AcqData" in dircontents:
        acqdata_path = os.path.join(path, "AcqData")
        if "sample_info.xml" in os.listdir(acqdata_path):
            tree = etree.parse(os.path.join(acqdata_path, "sample_info.xml"))
            root = tree.getroot()
            for samplefield in root.xpath('//Field[Name="Sample Position"]'):
                vialnum = samplefield.find("Value")
                if vialnum is not None and len(vialnum.text.split()) == 1:
                    metadata['vialpos'] = vialnum.text
                    return metadata

    # runstart.txt
    if "runstart.txt" in dircontents:
        f = open(os.path.join(path, "runstart.txt"))
        lines = f.read().splitlines()
        f.close()
        for line in lines:
            stripped = line.strip()
            if "Alsbottle" not in stripped:
                continue
            vialnum = stripped.split()[-1]
            if int(vialnum):
                metadata['vialpos'] = vialnum
                return metadata

    # RUN.LOG
    if "RUN.LOG" in dircontents:
        f = open(os.path.join(path, "RUN.LOG"), 'rb')
        plaintext = f.read().decode('ascii', 'ignore').replace("\x00", "")
        f.close()
        for line in plaintext.splitlines():
            vialpos = None
            if "Method started" in line:
                split = line.split()
                vialpos = get_nextstr(split, "vial#")
                if not vialpos:
                    vialpos = get_nextstr(split, "location")
            elif "Instrument running sample" in line:
                split = line.split()
                vialpos = get_nextstr(split, "Vial")
                if not vialpos:
                    vialpos = get_nextstr(split, "location")
                if not vialpos:
                    vialpos = get_nextstr(split, "sample")
            if vialpos:
                metadata['vialpos'] = vialpos.replace("'", "")
                break

    return metadata


def get_xml_vialnum(path):
    """
    Returns the VialNumber from an XML document, if it exists.

    Args:
        path (str): Path to the XML document.

    """
    tree = etree.parse(path)
    root = tree.getroot()
    for vialnum in root.xpath("//*[local-name()='VialNumber']"):
        if vialnum.text:
            return vialnum.text
    return None


def get_nextstr(str_list, target_str):
    """
    Returns the string at the next index in :obj:`str_list`, if it exists.

    Args:
        str_list (str): List of strings.
        target_str (str): Initial string to find.

    """
    try:
        next_str = str_list[str_list.index(target_str) + 1]
        return next_str
    except Exception:
        return None

if __name__ == "__main__":
    # Example usage
    path = "Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-03-03 19-34-20\GFP.D"

    paths = [path]
    for path in paths:
        datafiles = parse_allfiles(path)
        #metadata = parse_metadata(path, datafiles)
        #print(metadata)
        infos = [d.get_info() for d in datafiles]
        print(infos)
        datas=np.array([d.get_data() for d in datafiles])

        plt.plot(datas[2,:,0],datas[2, :,1]/np.amax(datas[2,:,1])*100)
        #plt.plot(datas[3,:,0],datas[3, :,1]/np.amax(datas[3,:,1])*100, color="blue")
        #plt.show()
        #print(datafiles[0].get_info())
    exit()

    print(datafiles[0].get_info())
    print(datafiles[0].extract_traces())
    print(datafiles[0].extract_traces(1))
    print(datafiles[0].extract_traces(2))
    print(datafiles[0].extract_traces([1,2]))

    #save csv file in parent path directory
    datafiles[0].export_csv('Z:\\Group Share\\Seun\\LX-MS\\LC\\LC Data\\OF_direct inject 2024-03-03 19-34-20\GFPtest.csv')


    #print(datafiles[0].metadata)
