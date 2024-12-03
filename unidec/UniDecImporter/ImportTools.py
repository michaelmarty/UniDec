import re

def header_test(path, deletechars=None, delimiter=" |\t|,", strip_end_space=True):
    """
    A quick test of a file to remove any problems from the header.

    If scans through each line of an input file specificed by path and count the number of lines that cannot be
    completely converted to floats. It breaks the loop on the first hit that is entirely floats, so any header lines
    blow that will not be counted.
    :param path: Path of file to be scanned
    :param deletechars: Characters to be removed from the file before scanning
    :param delimiter: Delimiter to be used in the file
    :param strip_end_space: Boolean to strip the end space from the line
    :return: Integer length of the number of header lines
    """
    header = 0
    try:
        with open(path, "r") as f:
            for line in f:
                # If the line is empty, skip it and add to the header
                if line == "\n":
                    header += 1
                    continue
                # Remove any characters that are defined in deletechars
                if deletechars is not None:
                    for c in deletechars:
                        line = line.replace(c, "")

                # Strip the end space if strip_end_space is True
                if strip_end_space:
                    if line[-1] == "\n" and len(line) > 1:
                        line = line[:-1]
                    if line[-1] == " " and len(line) > 1:
                        line = line[:-1]

                # Split the line by the delimiter and try to convert each element to a float
                for sin in re.split(delimiter, line):  # line.split(delimiter):
                    try:
                        float(sin)
                    except ValueError:
                        # print(sin, line)
                        header += 1
                        break
        # If the header is greater than 0, print the header length
        if header > 0:
            print("Header Length:", header)
    except (ImportError, OSError, AttributeError, IOError) as e:
        print("Failed header test", e)
        header = 0
    return int(header)

#
# def waters_convert2(path, config=None, outfile=None, time_range=None):
#     data = WDI(path).get_data(time_range=time_range)
#
#     if outfile is None:
#         outfile = os.path.join(path, "converted_rawdata.txt")
#     np.savetxt(outfile, data)
#     return data
