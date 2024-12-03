import os
# import numpy as np


class Importer:
    def __init__(self, file_path, **kwargs):
        self._file_path = file_path
        self._params = kwargs
        self.ext = os.path.splitext(file_path)[1]
        self.scans = None
        self.times = None
        self.centroided = False
        self.polarity = "Positive"

    def grab_data(self):
        pass

    def grab_scan_data(self, scan):
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

    def get_scan_time(self, scan):
        pass

    def get_ms_order(self, scan):
        pass

    def check_centroided(self):
        return self.centroided

