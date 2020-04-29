import numpy as np
from pymsfilereader import MSFileReader
import os
import time

def get_raw_metadata(path):
    rawfile = MSFileReader(path)
    #print(rawfile.Version())
    #for i in range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1):
    #    print(i, rawfile.GetCollisionEnergyForScanNum(i, MSOrder=1))

    #print(rawfile.GetComment1())

    # print('GetFileName', rawfile.GetFileName())
    # This one below seems wrong... 1970?
    print('GetCreationDate', rawfile.GetCreationDate())

    # Injection would be useful to addto our reports
    print('GetSeqRowInjectionVolume', rawfile.GetSeqRowInjectionVolume())

    # Instrument:
    print('GetInstModel', rawfile.GetInstModel())

    # Tried to find a way extract certain lines from this, but never got that
    print(rawfile.GetInstMethod())

    # These are all columns in the Xcalibur sequence editor with differnet names. They could be repurposed.
    print('GetSeqRowUserLabel', rawfile.GetSeqRowUserLabel(index=0))  # "Study" - could be solution info
    print('GetSeqRowUserLabel', rawfile.GetSeqRowUserLabel(
        index=1))  # "Client" - could be the collaborator's name. Might be able to insert this into title of the report.
    print('GetSeqRowUserLabel',
          rawfile.GetSeqRowUserLabel(index=2))  # "Laboratory" - one of these could be used as "shipment date"
    print('GetSeqRowUserLabel', rawfile.GetSeqRowUserLabel(index=3))  # "Company"
    print('GetSeqRowUserLabel', rawfile.GetSeqRowUserLabel(index=4))  # "Phone"

    # Not sure which comment is the one in the Xcalibur sequence editor
    print('GetSeqRowComment', rawfile.GetSeqRowComment())
    print('GetComment1()', rawfile.GetComment1())

    # Andrew - The date modified could be used as the "Data recorded" date as in a raw file, this date is the end of the run.
    # Python program to explain os.path.getmtime() method

    # Get the time of last modifation of the specified path since the epoch
    modification_time = os.path.getmtime(path)
    # print("Last modification time since the epoch:", modification_time)
    # convert the time in
    # seconds since epoch
    # to local time
    local_time = time.ctime(modification_time)
    print("Last modification time(Local time):", local_time)

    output_string="Parameter 1: " + str(1) +";"
    return output_string


if __name__ == "__main__":
    dir = "C:\\Python\\UniDec3\\TestSpectra"
    filename = "test.raw"
    path = os.path.join(dir, filename)
    get_raw_metadata(path)
