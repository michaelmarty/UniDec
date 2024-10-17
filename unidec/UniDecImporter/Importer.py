import os

import numpy as np



class Importer:

    def __init__(self, file_path, **kwargs):
        self._file_path = file_path
        self._params = kwargs

    def grab_data(self):
        pass

    def get_max_time(self):
        pass

    def get_max_scan(self):
        pass

    def get_scans_from_times(self, time_range):
        pass

    def get_times_from_scans(self, scan_range):
        pass

    def get_tic(self):
        pass

    def get_polarity(self):
        pass

