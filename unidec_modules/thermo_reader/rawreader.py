import numpy as np
from pymsfilereader import MSFileReader
import os
import time

def get_raw_metadata(path):
    rawfile = MSFileReader(path)
# Most are columns in the Xcalibur sequence editor
# "Company" - Research Group name
    company = rawfile.GetSeqRowUserText(index=3)
    print(company)
# "Laboratory" - collaborator PI's name
    lab = rawfile.GetSeqRowUserText(index=2)
    print(lab)
# "Client" - collaborator's name.
    client = rawfile.GetSeqRowUserText(index=1)
    print(client)
# "Phone" - This will be "shipment date"
    phone = rawfile.GetSeqRowUserText(index=4)
    shipment = phone.replace('Phone', '')
    print(shipment)
# Andrew - The date modified = "Data recorded" - this date is the end of the run.
    # Get the time of last modifation of the specified path since the epoch
    #modification_time = os.path.getmtime(path)
    # to local time
    local_time = time.strftime("%a, %d %b %Y, %H:%M", time.localtime(os.path.getmtime(path)))
    print(local_time)
# Instrument:
    inst = rawfile.GetInstModel()
    print(inst)
# "Study" - could be solution info - OBE mM AmAc
    study = rawfile.GetSeqRowUserText(index=0)
    print(study)
# Injection would be useful to add to our reports
    inject = rawfile.GetSeqRowInjectionVolume()
    print(inject)
    

    output_string="Lab: " + str(company) + ";" +"Client: " + str(lab)+ " - " + str(client) + ";" + "Shipment Date: " + str(shipment) + ";" + "Data Recorded: " + str(local_time) + ";" + "Instrument: " + str(inst) + ";" + "Study: " + str(study) + ";" + "Injection (uL): " + str(inject)
    #print(output_string)
    return output_string


if __name__ == "__main__":
    dir = "C:\\Users\\norri\\OneDrive\\Desktop\\UniDec_report tests"
    #dir = "D:\\Data"
    filename = "20200603_Ivan_D3-5A2_3k_HCD005.raw"
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
    filename = "20200603_Ivan_D3-5A2_3k_HCD005.raw"
    path = os.path.join(dir, filename)
    get_raw_companyname(path)
    
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

