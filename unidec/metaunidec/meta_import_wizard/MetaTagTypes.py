# the available tag types
# usage:
# import TagTypes as tt
# tt.SAMPLE

FILE_PATH = 'file_path'
FILE_NAME = 'file_name'
TIME_START = 'start'
TIME_END = 'end'
SCAN_START = 'scan_start'
SCAN_END = 'scan_end'
V1="v1"
V2="v2"
TYPE="type"

tag_type_formatting = {FILE_NAME: 'File Name',
                       FILE_PATH: 'Full Path',
                       TIME_START: 'Time Start',
                       TIME_END: 'Time End',
                       SCAN_START: 'Scan Start',
                       SCAN_END: 'Scan End',
                       V1:"Variable 1",
                       V2:"Variable 2"
                       }


def format_tag_name(tag_type):
    global tag_type_formatting
    try:
        return tag_type_formatting[tag_type]
    except KeyError:
        print('Unknown tag_type: ' + str(tag_type))
