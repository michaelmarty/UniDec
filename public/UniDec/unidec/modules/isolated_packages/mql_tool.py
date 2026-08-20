try:
    from massql import msql_engine
    from massql import msql_fileloading

    found = True
except ImportError:
    print("massql not found")
    found = False
    msql_engine = None
    msql_fileloading = None

import numpy as np
import os


class MQL_TOOL:
    def __init__(self, file):
        self.file = file
        if not found:
            print('Error, MassQL module not found')
            raise ImportError
        if not os.path.isfile(self.file):
            print("MassQL file not found:", self.file)
            print("Try selecting peaks first")
            raise FileNotFoundError

        self.ms1_df, self.ms2_df = msql_fileloading.load_data(self.file)

    def query(self, input_query, pks):
        self.input_query = input_query
        self.results_df = msql_engine.process_query(self.input_query, self.file, ms1_df=self.ms1_df, ms2_df=self.ms2_df)
        mz = self.results_df["comment"].to_numpy(dtype=float)
        for p in pks.peaks:
            if np.any(p.mass == mz):
                p.ignore = False
                print(p.mass, "Matched Query")
            else:
                p.ignore = True
                print(p.mass, "No Match")
        return pks
