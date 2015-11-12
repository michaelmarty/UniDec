# the available tag types
# usage:
# import TagTypes as tt
# tt.SAMPLE

SAMPLE = 'sample'
DESCRIPTION = 'description'
DRIFT_V = 'dv'
COLLISION_V = 'cv'
FILE_PATH = 'file_path'
FILE_NAME = 'file_name'
DRIFT_PRESSURE = 'drift_pressure'
TEMP = 'temp'
ATOM = 'atom'
WAVE_VELOCITY = 'wave_velocity'
WAVE_HEIGHT = 'wave_height'
EDC = 'edc'
CONE = 'cone'
EXTRACT = 'extract'
START = 'start'
END = 'end'
BIN = 'bin'
RANGE = 'range'
TYPE = 'type'
PUSHER = 'pusher'
FUNCTION = 'function'
SCAN_START = 'scan_start'
SCAN_END = 'scan_end'
TCAL1="tcal1"
TCAL2="tcal2"

tag_type_formatting = {SAMPLE: 'Sample',
                       DESCRIPTION: 'Description',
                       DRIFT_V: 'Drift Voltage',
                       COLLISION_V: 'Collision Voltage',
                       FILE_PATH: 'File Path',
                       FILE_NAME: 'File Name',
                       DRIFT_PRESSURE: 'Drift Pressure',
                       TEMP: 'Temperature',
                       ATOM: 'Atom',
                       WAVE_VELOCITY: 'Wave Velocity',
                       WAVE_HEIGHT: 'Wave Height',
                       EDC: 'EDC',
                       CONE: 'Cone Voltage',
                       EXTRACT: 'Extract',
                       START: 'Start',
                       END: 'End',
                       BIN: 'm/z Bin Size',
                       RANGE: 'm/z Range',
                       TYPE: 'Experiment Type',
                       PUSHER: 'Pusher',
                       FUNCTION: 'Function',
                       SCAN_START: 'Scan Start',
                       SCAN_END: 'Scan End',
                       TCAL1:"Calibration 1",
                       TCAL2:"Calibration 2"
                       }


def format_tag_name(tag_type):
    global tag_type_formatting
    try:
        return tag_type_formatting[tag_type]
    except KeyError:
        print 'Unknown tag_type: ' + str(tag_type)
