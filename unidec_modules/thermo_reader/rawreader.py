import numpy as np
from pymsfilereader import MSFileReader
import os
import time

def get_raw_metadata(path):
    rawfile = MSFileReader(path)
# Most are columns in the Xcalibur sequence editor
# "Laboratory" - collaborator PI's name 
    lab = rawfile.GetSeqRowUserLabel(index=2)
    lab1 = lab.replace('Laboratory', '')
    print(lab1)

# "Client" - collaborator's name.
    client = rawfile.GetSeqRowUserLabel(index=1)
    client1 = client.replace('Client', '')
    print(client1)    
# "Phone" - This will be "shipment date"
    phone = rawfile.GetSeqRowUserLabel(index=4)
    shipment = phone.replace('Phone', 'Shipment Date: ')
    print(shipment)  
# Andrew - The date modified = "Data recorded" - this date is the end of the run.
    # Get the time of last modifation of the specified path since the epoch
    #modification_time = os.path.getmtime(path)
    # to local time
    local_time = time.strftime("%a, %d %b %Y, %H:%M", time.gmtime(os.path.getmtime(path)))
    print(local_time)
# Instrument:
    inst = rawfile.GetInstModel()
    print(inst)
# "Study" - could be solution info - OBE mM AmAc
    study = rawfile.GetSeqRowUserLabel(index=0)
    study1 = study.replace('Study', 'Study: ')
    print(study1)
# Injection would be useful to add to our reports
    inject = rawfile.GetSeqRowInjectionVolume()
    print(inject)

    

    output_string="Client: " + str(lab1)+ "- " + str(client1) + ";" + str(shipment) + ";" + "Data Recorded: " + str(local_time) + ";" + "Instrument: " + str(inst) + ";" + str(study1) + ";" + "Injection (uL): " + str(inject)
    print(output_string)
    return output_string


if __name__ == "__main__":
    dir = "C:\\Users\\norri\\OneDrive\\Desktop"
    filename = "20200223_UHMR1_OBE_Yang_YH_c5-Hfuse_21.raw"
    path = os.path.join(dir, filename)
    get_raw_metadata(path)

    
    
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

