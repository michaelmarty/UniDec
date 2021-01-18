import numpy as np
from pymsfilereader import MSFileReader
import os
import time

def get_raw_metadata(path):
    rawfile = MSFileReader(path)
# Most are columns in the Xcalibur sequence editor the names are Xcalibur default.
# "Study" - could be solution info - OBE mM AmAc
    header1 = rawfile.GetSeqRowUserLabel(index=0)
    print(header1)
    study1 = rawfile.GetSeqRowUserText(index=0)
    print(study1)
# "Client" - collaborator's name, PI - contact.
    header2 = rawfile.GetSeqRowUserLabel(index=1)
    print(header2)
    client2 = rawfile.GetSeqRowUserText(index=1)
    print(client2)
# "Laboratory" - collaborator PI's name - contact name
    header3 = rawfile.GetSeqRowUserLabel(index=2)
    print(header3)
    lab3 = rawfile.GetSeqRowUserText(index=2)
    print(lab3)
# "Company" - Research Group name (previous) "Sample Conditions" (new)
    header4 = rawfile.GetSeqRowUserLabel(index=3)
    print(header4)
    company4 = rawfile.GetSeqRowUserText(index=3)
    print(company4)
# "Phone" - This will be "Shipment date"
    header5 = rawfile.GetSeqRowUserLabel(index=4)
    print(header5)
    phone5 = rawfile.GetSeqRowUserText(index=4)
    #shipment = phone.replace('Phone', '')
    print(phone5)
# Andrew - The date modified = "Data recorded" - this date is the end of the run.
    # Get the time of last modifation of the specified path since the epoch
    #modification_time = os.path.getmtime(path)
    # to local time
    local_time = time.strftime("%a, %d %b %Y, %H:%M", time.localtime(os.path.getmtime(path)))
    print(local_time)
# Instrument:
    inst = rawfile.GetInstModel()
    print(inst)
# Injection would be useful to add to our reports
    inject = rawfile.GetSeqRowInjectionVolume()
    print(inject)
# "comment" - for misc info
    GetSeqRowComment = rawfile.GetSeqRowComment()
    print(GetSeqRowComment)
    

    output_string= str(header3) + ": " + str(lab3) + ";" + str(header2) + ": " + str(client2) + ";" + str(header5) + ": " + str(phone5) + ";" \
                   + "Data Recorded: " + str(local_time) + ";" + "Instrument: " + str(inst) + ";" + str(header4) + ": " \
                  + str(company4) + ";" + str(header1) + ": " + str(study1) + ";" + "Injection (uL): " + str(inject) + ";" + "Comment: " + str(GetSeqRowComment)
    #print(output_string)
    return output_string


if __name__ == "__main__":
    dir = "C:\\Users\\norri\\OneDrive\\Desktop\\UniDec_report tests"
    #dir = "D:\\Data"
    filename = "20201022_Ryan_RK_hetcrn-n-s-2_3xdilu_6k_HCD120.raw"
    path = os.path.join(dir, filename)
    get_raw_metadata(path)

def get_raw_samplename(path):
    rawfile = MSFileReader(path)
# "sample name" - will be in the title of the report
    rawsamplename = rawfile.GetSeqRowSampleName()
    print(rawsamplename)
    output_string= str(rawsamplename)
    return output_string

if __name__ == "__main__":
    dir = "C:\\Users\\norri\\OneDrive\\Desktop\\UniDec_report tests"
    #dir = "D:\\Data"
    filename = "20201022_Ryan_RK_hetcrn-n-s-2_3xdilu_6k_HCD120.raw"
    path = os.path.join(dir, filename)
    #get_raw_companyname(path)
    
    #for i in range(rawfile.FirstSpectrumNumber, rawfile.LastSpectrumNumber + 1):
    #    print(i, rawfile.GetCollisionEnergyForScanNum(i, MSOrder=1))

    #print(rawfile.GetComment1())

    # print('GetFileName', rawfile.GetFileName())
    # This one below seems wrong... 1970?
    #print('GetCreationDate', rawfile.GetCreationDate())

# Tried to find a way extract certain lines from this, but never got that
    #print(rawfile.GetInstMethod())

# These are all columns in the Xcalibur sequence editor with differnet names. They could be repurposed.
    #print('GetSeqRowUserLabel',rawfile.GetSeqRowUserLabel(index=2))  # "Laboratory" - one of these could be used as "shipment date"
    #print('GetSeqRowUserLabel', rawfile.GetSeqRowUserLabel(index=3))  # "Company"
# Not sure which comment is the one in the Xcalibur sequence editor
    #print('GetSeqRowComment', rawfile.GetSeqRowComment())
    #print('GetComment1()', rawfile.GetComment1())
    #print('GetSeqRowSampleName', rawfile.GetSeqRowSampleName())
    #print('GetSeqRowRawFileName', rawfile.GetSeqRowRawFileName())
    #print('GetSeqRowLevelName', rawfile.GetSeqRowLevelName())

